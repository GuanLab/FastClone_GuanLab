"""Subclonal phylogeny inference"""

import logbook
import numpy
import pandas
import scipy.misc
import scipy.optimize
import scipy.signal
import scipy.stats


_LOG = logbook.Logger(__name__)


def infer(subclones):
    subclones = subclones[::-1]
    _LOG.info("Phylogeny inference model: Beta distribution")
    parents = numpy.zeros(subclones.size, dtype="i1")
    parents[0] = 0
    parents[1] = 0
    if parents.size < 3:
        _LOG.debug("Phylogeny: {}", parents)
        return (parents.size - parents - 1)[::-1]
    fractions = numpy.zeros(parents.shape, dtype="f4")
    leftover = numpy.copy(subclones)
    fractions[0] = leftover[0]
    fractions[1] = leftover[1] / leftover[0]
    leftover[0] -= leftover[1]
    parents[2] = _assign_3rd_subclone(leftover, fractions)
    for j in range(3, parents.size):
        parents[j] = _assign(j, leftover, fractions)
    _LOG.debug("Phylogeny: {}", parents + 1)
    return (parents.size - parents - 1)[::-1]

def _assign_3rd_subclone(leftover, fractions):
    if leftover[0] > leftover[2] + 0.1:
        fractions[2] = leftover[2] / leftover[0]
        leftover[0] -= leftover[2]
        return 0
    else:
        fractions[2] = leftover[2] / leftover[1]
        leftover[1] -= leftover[2]
        return 1

def _assign(idx, leftover, fractions):
    check = leftover[idx] < leftover[:idx]
    new_frac = leftover[idx] / leftover[check]
    fracs = numpy.hstack((
        numpy.repeat(fractions[:idx, None], new_frac.size, axis=1),
        new_frac,
    ))
    likelihood = numpy.asarray([
        scipy.stats.beta(
            *scipy.stats.beta.fit(fracs[:, j], f0=1, floc=0, fscale=1)
        ).logpdf(fracs[:, j]).sum()
        for j in range(new_frac.size)
    ])
    parent = numpy.where(check)[0][likelihood.argmax()]
    fractions[idx] = leftover[idx] / leftover[parent]
    leftover[parent] -= leftover[idx]
    return parent
