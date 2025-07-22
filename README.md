# panNrefNkmer

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

# Overview

A snakemake workflow to analyze short reads with a either a pangenome reference, linear reference, or no reference

# Setup

1. Make sure you have [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) installed. 

2. Install snakemake and slurm executer using the provided yaml file

```
mamba env create --name snakemake --file setup.yaml
```

3. Grab repository from github

```
git clone https://github.com/milesroberts-123/panNrefNkmer.git
```

4. Modify snakemake profile. The default profile runs snakemake on a slurm cluster (`workflow/profiles/default/config.yaml`), but you should still change the slurm account, slurm partition, and default resources to match your system.

# Inputs

## 1. Configuration

All variables are collected in `config/config.yaml`. Variable defintions are given in more detial in `config/config.schema.yaml`.

## 2. Reference genome sequences and annotations

For every genome specified in `config/config.yaml`, add a fasta file to `config/linear_genomes/sequence/` and a gff file to `config/linear_genomes/annotation/`.

For example, if the `linrefs` variable in your config file looks like this:

```
linrefs:
 - "ah7"
 - "arb0"
 - "belmonte494"
```

Then your file paths should look like this:

```
config/linear_genomes/
├── annotation
│   ├── ah7.gff
│   ├── arb0.gff
│   ├── belmonte494.gff
└── sequence
    ├── ah7.fa
    ├── arb0.fa
    ├── belmonte494.fa
```

# Usage

## Run whole workflow with conda envs on slurm cluster

`snakemake --sdm conda --rerun-incomplete --rerun-triggers mtime --scheduler greedy --retries 1 --keep-going`
