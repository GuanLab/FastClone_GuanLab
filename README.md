# FastClone

FastClone is a fast algorithm to infer tumour heterogeneity. Given somatic
mutation frequencies and copy number data, FastClone infers subclonal
composition and phylogeny. The algorithm won the first place in DREAM Somatic
Mutation Calling -- Heterogeneity Challenge in 2016.

## Installation

FastClone needs Python 3.5 or higher versions. It needs logbook, python-fire,
scikit-learn, and pandas. To install the package using Pip,

```
git clone https://github.com/GuanLab/FastClone.git
pip install FastClone
```

## Usage

FastClone accepts either MuTect VCF + Battenberg format (specified in the DREAM
SMC-Het Challenge) or PyClone format.

An example to load samples and infer:
```
fastclone load-pyclone t1 t1.tsv 0.8 solve [output_pathway]
```

Run `fastclone` for more help information.

If MuTect VCF and PyClone samples are provided, note that MuTect
mutations are labelled as 'Chromosome:Coordinate:AltBase', such as
'Y:15989697:G'. Make sure PyClone ID uses the same ID.
