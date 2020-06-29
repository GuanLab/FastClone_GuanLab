"""Subclonal composition inference"""

import pickle

import logbook
import numpy
import pandas
import scipy.special
import scipy.optimize
import scipy.signal
import scipy.stats
from itertools import chain


_LOG = logbook.Logger(__name__)


_MINIMAL_MUTATIONS = 50

_SMALLEST_CLONE = 10

_MINIMAL_SUBCLONE = 1e-7


def infer_single_sample(mutations, purity):
    """Infer subclones for a single sample.

    Parameters
    ----------
    mutations : pandas.DataFrame
        The mutation table.
    purity : float | str | NoneType
        The estimated tumour purity.

    Returns
    -------
    numpy.ndarray[float]
        The subclone table.
    numpy.ndarray[float]
        The mutation assignment to subclones.
    """
    if isinstance(purity, str):
        peaks = numpy.asarray([float(x) for x in purity.split(':')])
        return _assign_mutations_to_subclones(mutations, peaks)
    filter_ = ((mutations['major_copy1'] == mutations['minor_copy1']) &
               (mutations['state1'] == 1.0))
    if filter_.sum() < _MINIMAL_MUTATIONS:
        filter_ = mutations['state1'] == 1.0
    if filter_.sum() > _MINIMAL_MUTATIONS:
        peaks = _infer_without_2state_cna(mutations[filter_], purity)
    else:
        peaks = _infer_with_2state_cna(mutations, purity)
    return _assign_mutations_to_subclones(mutations, peaks)


def _infer_without_2state_cna(mutations, purity):
    """Infer subclones for a single sample.

    Parameters
    ----------
    mutations : pandas.DataFrame
        The mutation table.
    purity : float | NoneType
        The estimated tumour purity.

    Returns
    -------
    numpy.ndarray[float]
        The subclone cell prevalance.
    """
    r = 1.0 if purity is None else purity
    f = mutations['allelic_count'] / mutations['total_count']
    ratios = f * ((mutations['major_copy1'] + mutations['minor_copy1']) * r
                  + mutations['normal_copy'] * (1 - r))
    peaks = _get_density_peaks(ratios)
    if purity is not None:
        peaks = _adjust_subclone_on_purity(peaks, purity)
    else:
        if peaks[-1] >= 1.1:
            peaks = [p for p in peaks if p < 1.1]
        if peaks[-1] >= 1:
            peaks = [p for p in peaks if p < 1.0] + [1.0]
        r = peaks[-1]
        ratios = f * ((mutations['major_copy1'] + mutations['minor_copy1']) * r
                      + mutations['normal_copy'] * (1 - r))
        peaks = _get_density_peaks(ratios)
        peaks = _adjust_subclone_on_purity(peaks, r)
    return numpy.asarray(peaks)


def _infer_with_2state_cna(mutations, purity):
    """Infer subclones for a single sample.

    Parameters
    ----------
    mutations : pandas.DataFrame
        The mutation table.
    purity : float | NoneType
        The estimated tumour purity.

    Returns
    -------
    numpy.ndarray[float]
        The subclone cell prevalance.
    """
    r = 1.0 if purity is None else purity
    f = mutations['allelic_count'] / mutations['total_count']
    ratios = []
    for i in range(len(mutations)):
        tmp_ratios = [f[i] * ((mj_copy + mutations['minor_copy1'].iloc[i]) * r + mutations['normal_copy'].iloc[i] * (1 - r))
                      for mj_copy in range(int(mutations['major_copy1'].iloc[i]))]
        ratios.append(tmp_ratios)
    ratios = list(chain(*ratios))
    ratios = pandas.Series(ratios)
    peaks = _get_density_peaks(ratios)
    if purity is not None:
        peaks = _adjust_subclone_on_purity(peaks, purity)
    else:
        if peaks[-1] >= 1.1:
            peaks = [p for p in peaks if p < 1.1]
        if peaks[-1] >= 1:
            peaks = [p for p in peaks if p < 1.0] + [1.0]
        r = peaks[-1]
        ratios = f * ((mutations['major_copy1'] + mutations['minor_copy1']) * r
                      + mutations['normal_copy'] * (1 - r))
        peaks = _get_density_peaks(ratios)
        peaks = _adjust_subclone_on_purity(peaks, r)
    return numpy.asarray(peaks)


def _get_density_peaks(samples, grid_count=200):
    """Identify density peaks

    Parameters
    ----------
    samples : numpy.ndarray[float]
        The samples for estimating density.
    grid_count : int
        The number of data point for scanning maxima. By default, it is 200.

    Returns
    -------
    list[float]
        The subclone cell prevalance.
    """
    kde = scipy.stats.gaussian_kde(samples.values.flatten())
    delta = 0.5 / grid_count
    x = numpy.linspace(delta, 1.0 - delta, grid_count)
    logpdf = kde.logpdf(x)
    maxidx = list(scipy.signal.argrelmax(logpdf)[0])
    if logpdf[0] > logpdf[1]:
        maxidx.append(0)
    if logpdf[-1] > logpdf[-2]:
        maxidx.append(len(logpdf) - 1)
    peaks = [float(scipy.optimize.fminbound(lambda x: -kde.logpdf(x),
                                            x[idx] - delta, x[idx] + delta,
                                            xtol=1e-6, full_output=False,
                                            disp=False))
             for idx in maxidx]
    peaks.sort()
    _LOG.info('Found peaks: {}', peaks)
    return peaks


def _adjust_subclone_on_purity(peaks, purity):
    """Adjust peaks based on tumour purity.

    Parameters
    ----------
    peaks : list[float]
        The identified density peaks of subclones.
    purity : float | NoneType
        The estimated tumour purity.

    Returns
    -------
    list[float]
        The adjusted density peaks of subclones
    """
    if purity is not None:
        max_frac = purity
        diff = numpy.abs(numpy.asarray(peaks) - max_frac)
        max_peak_idx = diff.argmin()
        if diff[max_peak_idx] > 0.1:
            peaks = peaks + [max_frac]
        else:
            peaks = peaks[:max_peak_idx + 1]
    if peaks[-1] >= 1.1:
        peaks = [p for p in peaks if p < 1.1]
    if peaks[-1] >= 1:
        peaks = [p for p in peaks if p < 1.0] + [1.0]
    _LOG.info('Adjusted peaks to: {}', peaks)
    return peaks


def _assign_mutations_to_subclones(mutations, peaks):
    """Assign mutations to subclones.

    Parameters
    ----------
    mutations : pandas.DataFrame
        The mutations for estimating subclone weights.
    peaks : list[float]
        The identified density peaks of subclones.

    Returns
    -------
    numpy.ndarray[int]
        The adjusted subclones for all mutations.
    numpy.ndarray[float]
        The probability of mutations associated to subclones.
    """
    _LOG.debug('Found raw subclone prevalences: {}', peaks)
    while True:
        weights, scores = _estimate_subclone_weights(mutations, peaks)
        assignment = scores.values.argmax(axis=1)
        idx, size = numpy.unique(assignment, return_counts=True)
        idx = idx[size >= _SMALLEST_CLONE]
        if idx.size == peaks.size or len(mutations) < _MINIMAL_MUTATIONS * 2:
            return peaks, scores
        peaks = peaks[idx]


def _estimate_subclone_weights(mutations, peaks):
    """Estimate subclone weights.

    Parameters
    ----------
    mutations : pandas.DataFrame
        The mutations for estimating subclone weights.
    peaks : list[float]
        The identified density peaks of subclones.

    Returns
    -------
    numpy.ndarray[float]
        The weights of subclones.
    numpy.ndarray[int]
        The assignment of mutations to subclones.
    """
    peaks = peaks[None, :, None]
    purity = peaks.max()
    weights = numpy.repeat(1 / peaks.size, peaks.size).astype('f4')
    major = mutations['major_copy1'].values[:, None, None]
    minor = mutations['minor_copy1'].values[:, None, None]
    state1 = mutations['state1'].values[:, None, None]
    denominator = state1 * purity * (major + minor) + (1 - state1 * purity) * 2
    numerator = numpy.dstack([
        numpy.where(peaks < state1, 0, (major * state1 * purity +
                                        peaks - state1 * purity)),
        numpy.where(peaks < state1, 0, (minor * state1 * purity +
                                        peaks - state1 * purity)),
        numpy.repeat(numpy.arange(1, (numpy.nanmax(major) + 1)) * peaks, len(mutations), axis=0)
    ])
    copy_number_category = numpy.unique(major)
    for copy_number in copy_number_category:
        numerator[numpy.argwhere(major == copy_number)[0][0]][0][(copy_number + 2):] = 0
    p = numerator / denominator
    x = mutations['allelic_count'].values[:, None, None]
    k = mutations['total_count'].values[:, None, None]
    l = numpy.nanmax(scipy.stats.binom.logpmf(x, k, p), axis=2) + numpy.log(weights)
    l_sample = scipy.special.logsumexp(l, axis=1)
    l_sum = numpy.nansum(l_sample)
    for t in range(5000):
        old_l_sum = l_sum
        scores = numpy.exp(l - l_sample[:, None])
        weights = numpy.nanmean(scores, axis=0)
        weights /= weights.sum()
        l = numpy.nanmax(scipy.stats.binom.logpmf(x, k, p), axis=2) + numpy.log(weights)
        l_sample = scipy.special.logsumexp(l, axis=1)
        l_sum = numpy.nansum(l_sample)
        if l_sum > old_l_sum and l_sum < old_l_sum + 0.001:
            _LOG.info('EM is fully converged.')
            break
    else:
        _LOG.warning('Maximum iteration reached but the EM did not converged.')
    scores = pandas.DataFrame(numpy.exp(l - l_sample[:, None]),
                              index=mutations.index, dtype='f4')
    return weights, scores
