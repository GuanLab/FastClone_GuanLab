"""Various input loaders."""


import pandas


__all__ = (
    'load_mutect_and_battenberg',
    'load_pyclone',
)


def load_mutect_and_battenberg(vcf_path, vcf_sample, battenberg_path):
    """Load SNVs and CNAs from MuTect VCF and Battenberg files.

    The file format is used in the DREAM SMC-Het Challenge.

    Parameters
    ----------
    vcf_path : file-like object, pathlib.Path or str
        A path or a file handle of the VCF file.
    vcf_sample : str
        The column of sample in the VCF file to read.
    battenberg_path : file-like object, pathlib.Path or str
        A path or a file handle of the Battenberg file.

    Returns
    -------
    pandas.DataFrame
        The mutation data.
    float | NoneType
        The estimated purity.
    """
    vcf = pandas.read_csv(vcf_path, comment='#', sep='\t', header=None)
    vcf.columns = _read_mutect_vcf_columns(vcf_path)
    data = pandas.DataFrame()
    data['chromosome'] = vcf['CHROM'].astype(str)
    data['coordinate'] = vcf['POS']
    data['base'] = vcf['ALT']
    sample = vcf[vcf_sample].str.split(':', expand=True)
    read_counts = sample[1].str.split(',', expand=True).astype('i4')
    data['total_count'] = read_counts[0] + read_counts[1]
    data['allelic_count'] = read_counts[1]
    data['normal_copy'] = 2
    data['major_copy1'] = 1
    data['minor_copy1'] = 1
    data['major_copy2'] = 1
    data['minor_copy2'] = 1
    data['state1'] = 1.0
    if (data['chromosome'] == 'Y').any():
        data.loc[data['chromosome'] == 'X', 'normal_copy'] = 1
        data.loc[data['chromosome'] == 'Y', 'normal_copy'] = 1
        data.loc[data['chromosome'] == 'X', 'minor_copy1'] = 0
        data.loc[data['chromosome'] == 'Y', 'minor_copy1'] = 0
        data.loc[data['chromosome'] == 'X', 'minor_copy2'] = 0
        data.loc[data['chromosome'] == 'Y', 'minor_copy2'] = 0
    purity = _load_battenberg_cnas(data, battenberg_path)
    data.index = (data['chromosome'] + ':' + data['coordinate'].astype(str)
                  + ':' + data['base'])
    data.drop(['chromosome', 'coordinate', 'base'], axis=1, inplace=True)
    data = data[data.allelic_count/data.total_count > 0.01]
    return data, purity


def _read_mutect_vcf_columns(vcf_path, add_smc_het_column=False):
    """Read MuTect VCF column names.

    Parameters
    ----------
    vcf_path : file-like object, pathlib.Path or str
        A path or a file handle of the VCF file.
    add_smc_het_column : bool
        Whether to add a column at the end for the somatic mutation
        calling truth. This is only for SMC-Het challenge truth files.
        It is disabled by default.

    Returns
    -------
    list[str]
        List of VCF column names.
    """
    for line in open(vcf_path):
        if line[1] != '#' and line[0] == '#':
            columns = line[1:].strip().split()
            if add_smc_het_column:
                columns.append(None)
            return columns
    raise RuntimeError('VCF header line not found')


def _load_battenberg_cnas(data, path):
    """Load CNA information from the Battenberg file.

    Update the mutation data with new copy number data. Also, estimate
    the purity of the tumour tissue.

    Parameters
    ----------
    data : pandas.DataFrame
        The mutation data.
    path : file-like object, pathlib.Path or str
        A path or a file handle of the Battenberg file.

    Returns
    -------
    float
        The estimated purity of the tumour tissue.
    """
    table = pandas.read_csv(path, sep='\t')
    table['chr'] = table['chr'].astype(str)
    if 'Y' in table['chr']:
        data.loc[data['chromosome'] == 'X', 'normal_copy'] = 1
        data.loc[data['chromosome'] == 'Y', 'normal_copy'] = 1
    purities = []
    for __, row in table.iterrows():
        purity = _estimate_purity_from_battenberg(row)
        if purity is not None:
            purities.append(purity)
        filter_ = ((data['chromosome'] == row['chr'])
                   & (data['coordinate'] > row['startpos'])
                   & (data['coordinate'] <= row['endpos']))
        if row['frac1_A'] == 1.0:
            data.loc[filter_, 'major_copy1'] = int(row['nMaj1_A'])
            data.loc[filter_, 'minor_copy1'] = int(row['nMin1_A'])
        else:
            data.loc[filter_, 'major_copy1'] = int(row['nMaj1_A'])
            data.loc[filter_, 'minor_copy1'] = int(row['nMin1_A'])
            data.loc[filter_, 'state1'] = row['frac1_A']
            data.loc[filter_, 'major_copy2'] = int(row['nMaj2_A'])
            data.loc[filter_, 'minor_copy2'] = int(row['nMin2_A'])
    if purities:
        return sum(purities) / len(purities)
    else:
        return None


def _estimate_purity_from_battenberg(entry):
    """Estimate the tumour purity from the Battenberg entry.

    Parameters
    ----------
    entry : pandas.Series
        The Battenberg copy number entry.

    Returns
    -------
    float | NoneType
        The estimated purity of the tumour tissue.
    """
    if ((entry['frac1_A'] == 1.0 and entry['nMaj1_A'] == entry['nMin1_A'])
            or (entry['frac1_A'] != 1.0
                and entry['nMaj1_A'] == entry['nMin1_A']
                and entry['nMaj2_A'] == entry['nMin2_A'])):
        return None
    if entry['frac1_A'] == 1.0:
        return ((1 - 2 * entry["BAF"])
                / (1 + entry["BAF"] * (entry["nMaj1_A"] + entry["nMin1_A"])
                   - 2 * entry["BAF"] - entry["nMaj1_A"]))
    return ((1 - 2 * entry["BAF"])
            / (1 + entry["BAF"] * (entry["frac1_A"] * (entry["nMaj1_A"]
                                                       + entry["nMin1_A"])
                                   + entry["frac2_A"] * (entry["nMaj2_A"]
                                                         + entry["nMin2_A"]))
               - 2 * entry["BAF"] - entry["nMaj1_A"] * entry["frac1_A"]
               - entry["nMaj2_A"] * entry["frac2_A"]))


def load_pyclone(path):
    """Load SNVs and CNAs from the PyClone file.

    Parameters
    ----------
    path : file-like object, pathlib.Path or str
        A path or a file handle of the PyClone input file.
    battenberg_path : file-like object, pathlib.Path or str
        A path or a file handle of the Battenberg file.

    Returns
    -------
    pandas.DataFrame
        The mutation data.
    """
    data = pandas.read_csv(path, sep='\t')[[
        'mutation_id', 'ref_counts', 'var_counts', 'normal_cn',
        'major_cn', 'minor_cn',
    ]]
    data.set_index('mutation_id', inplace=True)
    data.index.name = None
    data['ref_counts'] += data['var_counts']
    data.columns = ['total_count', 'allelic_count', 'normal_copy',
                    'major_copy1', 'minor_copy1']
    data['major_copy2'] = 1
    data['minor_copy2'] = 1
    data['state1'] = 1.0
    data = data[data.allelic_count/data.total_count > 0.01]
    return data
