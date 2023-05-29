
# getRPercId: a tool to extract references from BAM files and get read statistics


[![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/genomewalker/get-read-percid?include_prereleases&label=version)](https://github.com/genomewalker/get-read-percid/releases) [![get-read-percid](https://github.com/genomewalker/get-read-percid/workflows/getRPercId_ci/badge.svg)](https://github.com/genomewalker/get-read-percid/actions) [![PyPI](https://img.shields.io/pypi/v/get-read-percid)](https://pypi.org/project/get-read-percid/) [![Conda](https://img.shields.io/conda/v/genomewalker/get-read-percid)](https://anaconda.org/genomewalker/get-read-percid)

A simple tool to extract references from BAM files and get read statistics

# Installation

We recommend having [**conda**](https://docs.conda.io/en/latest/) installed to manage the virtual environments

### Using pip

First, we create a conda virtual environment with:

```bash
wget https://raw.githubusercontent.com/genomewalker/get-read-percid/master/environment.yml
conda env create -f environment.yml
```

Then we proceed to install using pip:

```bash
pip install get-read-percid
```

### Using mamba

```bash
mamba install -c conda-forge -c bioconda -c genomewalker get-read-percid
```

### Install from source to use the development version

Using pip

```bash
pip install git+https://git@github.com/genomewalker/get-read-percid.git
```

By cloning in a dedicated conda environment

```bash
git clone https://github.com:genomewalker/get-read-percid.git
cd get-read-percid
conda env create -f environment.yml
conda activate get-read-percid
pip install -e .
```


# Usage

getRPercId will take a BAM file and/or a list in the references in a BAM file and calculate the identity to the reference. In addition, if the set of references has been listed it will be extracted in an individual BAM file. The reference file can have an extra column in case the reference name has to be renamed, i.e., we want to profile it in [anvi'o](https://anvio.org/).

For a complete list of options:

```bash
$ getRPercId --help

usage: getRPercId [-h] -b FILE [--reference-list FILE] [--read-list FILE] [--prefix STR]
                  [--sort-memory STR] [--threads INT] [--chunk-size INT] [--debug] [--version]

A simple tool to extract references from BAM files and get read statistics

options:
  -h, --help            show this help message and exit

required arguments:
  -b FILE, --bam FILE   The BAM file used to generate the metaDMG results (default: None)
  --reference-list FILE
                        A list of references that we want to extract from the BAM file
                        (default: None)
  --read-list FILE      A list of reads that we want to extract from the BAM file (default:
                        None)

optional arguments:
  --prefix STR          Prefix used for the output files (default: None)
  --sort-memory STR     Set maximum memory per thread for sorting; suffix K/M/G recognized
                        (default: 1G)
  --threads INT         Number of threads (default: 1)
  --chunk-size INT      Chunk size for parallel processing (default: None)
  --debug               Print debug messages (default: False)
  --version             Print program version
```

One would run `getRPercId` as:

```bash
getRPercId --bam 3e3ce8e9f7.94.dmg-reads.dedup.filtered.bam --reference-list ref-list.csv --prefix 3e3ce8e9f7
```

The reference list file is a CSV file with the following format:

```bash
$ cat ref-list.csv
3300025461_7,c_000000000001
```

The program will output a CSV file with the following format when using a reference list:

```bash
$ gzcat 3e3ce8e9f7___c_000000000001--percid.csv.gz | head -n 2
M_A00706:51:HK27GDSXX:1:1342:26467:11694,99.03846153846155
M_A00706:51:HK27GDSXX:3:2465:29243:32690,97.45222929936305
```

And a CSV file with the following format when not using a reference list:

```bash
$ gzcat 3e3ce8e9f7.percid.csv.gz | head -n 2
M_A00706:51:HK27GDSXX:1:2134:16495:25582,3300000854_1,95.1219512195122
M_A00706:51:HK27GDSXX:3:1640:27950:10410,3300000854_1,97.5609756097561
```

where the first column is the read name, the second column is the reference name, and the third column is the identity to the reference.

One can also can provide a list of reads to extract from the BAM file. The read list file is a CSV file with the following format:

```txt
3300025461_7,M_A00706:51:HK27GDSXX:1:1101:10031:18865,c_000000000001
3300025461_7,M_A00706:51:HK27GDSXX:1:1101:10113:5384,c_000000000001
3300025461_7,M_A00706:51:HK27GDSXX:1:1101:10321:28573,c_000000000001
3300025461_7,M_A00706:51:HK27GDSXX:1:1101:10484:13510,c_000000000001
```

Where the first column is the reference, the second column is the read name, and the third column is the new reference name we want to give. This last column is optional.