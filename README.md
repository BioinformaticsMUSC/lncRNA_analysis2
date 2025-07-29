# lncRNA_analysis2

A Snakemake pipeline for lncRNA analysis. This pipeline is an updated version of the TGAC lncRNA_analysis pipeline here: https://github.com/TGAC/lncRNA-analysis/(https://github.com/TGAC/lncRNA-analysis/)

---

## Features

- Trims and performs FastQC on fastq files
- Aligns to transcriptome using STAR aligner  
- Finds Fasta sequences from bam alignments
- Runs CPAT(https://cpat.readthedocs.io/en/latest/) to find non-coding probabilities of transcripts
- Runs CPC(https://cpc.gao-lab.org/) to calculate coding potentials of transcripts
- Finds concordant non-coding transcripts between CPC and CPAT
- Generates GTFs for quantification and lists of non-coding genes

---

## Installation

SnakeLord requires Python 3.9+.

Using pip:

```bash
git clone https://github.com/BioinformaticsMUSC/lncRNA_analysis2.git
cd lncRNA_analysis2
```

lncRNA_analysis2 will create conda environments as needed, but to run snakemake, use the following conda environment:
```bash
conda create -n snakemake_env python=3.10 snakemake
```

## Config

Before running, edit the config.yaml and samples.tsv files.


## Usage

To start the pipeline, navigate to the `lncRNA_analysis2` folder and run:
```bash
snakemake --cores 4 --use-conda
```

For sharing conda environments, you can run:
```bash
snakemake --cores 4 --use-conda --conda-prefix /path/to/shared/envs/folder
```

