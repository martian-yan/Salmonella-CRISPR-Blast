# Salmonella-CRISPR-Blast

A Python pipeline for fast identification of Salmonella lineages from CRISPR arrays.

## Dependencies

- blast+
- bedtools
- Python >= 3.9
- Python modules: (will be installed by `pip` automatically)
    - biopython
    - pandas

## Installation

### Install from source

Prepare a proper environment with the necessary dependencies, clone this project and install:

```
git clone https://github.com/martian-yan/Salmonella-WGCT-pipeline.git
cd Salmonella-WGCT-pipeline
python -m build
pip install dist/SWGCT-1.0.0-py3-none-any.whl
```

If anything goes wrong, check if you have the latest version of Python `setuptools` and `build`:

```
pip install --upgrade setuptools
pip install --upgrade build
```

## Usage

`Salmonella-CRISPR-Blast` includes 2 subcommands

1. Use `swgct run` to identify the CRISPR sequences and spacer arrays from each *Salmonella* genome
2. Then use `swgct summary` to summarise the CRISPR types from a batch of results.

## Paraments

### swgct run

```
$ swgct run -h
usage: swgct run [-h] [-o outdir] [-p prefixed] [-d db_path] [--algorithm {blast,minced}] [-k] fasta

optional arguments:
  -h, --help            show this help message and exit

Input/Output:
  fasta                 Input FASTA file (required)
  -o outdir, --outdir outdir
                        Output directory. By default it's a folder with the same name as input
  -p prefixed, --prefix prefixed
                        Prefixed output file base names. By default it's the same with '--outdir'

Other:
  -d db_path, --database_path db_path
                        BLAST databases of Salmonella CRISPR DR and spacers
  --algorithm {blast,minced}
                        The algorithm used for searching CRISPR alleles. They may generate different results. Default="blast". "minced" is an option when the "blast" does not work.
  -k, --keep            Keep the tmp files, which will be cleaned by default
```

### swgct summary

```
$ swgct summary -h
usage: swgct summary [-h] [-d dir]

Summarize the result from SWGCT, expected to be run under a directory containing all the results from 'swgct run'.

optional arguments:
  -h, --help         show this help message and exit
  -d dir, --dir dir  Directory to work with, default is current directory.
```