import pathlib
import sys

import fire
import logbook
import pandas

from . import loaders
from . import phylogeny
from . import subclone


_LOG = logbook.Logger(__name__)


class Entrypoint:
    """FastClone: a tumour subclone inference tool

    Check subcommands for sample loading.
    Chain multiple subcommands load multiple samples and infer.
    Use '-- --help' to show help messages.

    Examples
    --------
    fastclone load-pyclone t1 t1.tsv None load-pyclone t2 t2.tsv 0.8 solve

    Notes
    -----
    If MuTect VCF and PyClone samples are provided, note that MuTect
    mutations are labelled as 'Chromosome:Coordinate:AltBase', such as
    'Y:15989697:G'. Make sure PyClone ID uses the same ID.
    """

    def __init__(self):
        self._samples = dict()

    def load_mutect_battenberg(self, name, vcf_path, vcf_sample,
                               battenberg_path):
        """Load a sample from MuTect and Battenberg files.

        Parameters
        ----------
        name : str
            The name of the sample.
        vcf_path : str
            The MuTect VCF file.
        vcf_sample : str
            The column of sample in the VCF file to read.
        battenberg_path : str
            The Battenberg file.
        """
        if name in self._samples:
            raise ValueError('duplicate sample: ' + name)
        _LOG.info('Loading {} ...', name)
        _LOG.info('VCF file {} (column: {})', vcf_path, vcf_sample)
        _LOG.info('Battenberg file {}', battenberg_path)
        self._samples[name] = loaders.load_mutect_and_battenberg(
            vcf_path, vcf_sample, battenberg_path,
        )
        return self

    def load_pyclone(self, name, path, purity):
        """Load a sample from a PyClone file.

        Parameters
        ----------
        path : str
            A path or a file handle of the VCF file.
        purity : float
            An estimated tumour turity. If unknown, use None.
        """
        if name in self._samples:
            raise ValueError('duplicate sample: ' + name)
        _LOG.info('Loading {} ...', name)
        _LOG.info('PyClone file {}', path)
        if purity is None:
            _LOG.info('Infer the purity later.')
        else:
            _LOG.info('Purity: {}', purity)
        self._samples[name] = (loaders.load_pyclone(path), purity)
        return self

    def load_pyclone_truth(self, name, path, truth):
        """Load a sample from a PyClone file with pickled data.

        Parameters
        ----------
        path : str
            A path or a file handle of the VCF file.
        truth : str
            A truth subclonal frequency file.
        """
        if name in self._samples:
            raise ValueError('duplicate sample: ' + name)
        _LOG.info('Loading {} ...', name)
        _LOG.info('PyClone file {}', path)
        _LOG.info('Truth file: {}', truth)
        self._samples[name] = (loaders.load_pyclone(path), truth)
        return self

    def solve(self, output):
        """Infer the subclonal composition.

        Parameters
        ----------
        output : str
            The output path.
        """
        if not self._samples:
            raise RuntimeError('inference needs at least one sample')
        if len(self._samples) > 1:
            raise NotImplementedError('mult-sample inference is not available')
        _LOG.info('Solving subclonal composition ...')
        pathlib.Path(output).mkdir(parents=True)
        samples = list(self._samples.items())
        subclones = None
        for name, sample in samples:
            sample_subclones, score = subclone.infer_single_sample(*sample)
            sample_subclones = pandas.DataFrame(sample_subclones,
                                                columns=[name])
            if subclones is None:
                subclones = sample_subclones
        phylogeny.infer(subclones, score, output)
        pandas.DataFrame(subclones).to_csv(output + '/subclones.csv')
        pandas.DataFrame(score).to_csv(output + '/scores.csv')


def main():
    with logbook.StderrHandler(level="INFO").applicationbound():
        fire.Fire(Entrypoint)


if __name__ == '__main__':
    sys.exit(main())
