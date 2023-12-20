# 1 About

CutOneStrand is a short pipeline designed to search for genomic positions that can be used
for targeting a specific strand with a cas9. It produces as an output a
list of SNPs that can be targeted on purpose after patient genotyping.
It has been designed in the context of [this
study](https://insermbiblio.inist.fr) and is especially usefull in the
context of [autosomal dominant disorders](https://www.genome.gov/genetics-glossary/Autosomal-Recessive-Disorder). For the moment it is only compatible with
spcas9 targeting NGG pam site. But further developments will lead to
other cas and pam to be added.

If you found this pipeline useful for you work. Please cite [Name of the
study](https://insermbiblio.inist.fr)

# 2 Installation

## 2.1 Requirements

The pipeline requires a Linux sytem with conda properly installed. If not
already set up, please refer to the [miniconda
documentation](https://docs.conda.io/en/latest/miniconda.html).

## 2.2 Get the code

Then clone this repository

```bash
    git clone git@github.com:clbenoit/CutOneStrand.git
```

## 2.3 Setup configuration

open **main.sh** file in **scripts** folder and replace **MYCONDAPATH** with your conda path installation. If you don't know where to find it execute :

```bash
    echo `which conda`
```

# 3 How to use

Launch the pipeline : `bash scripts/main.sh [args]`

```bash
############################################ HELP #####################################################

                 CutOneStrand version = 1.0.0

Usage: scripts/main.sh [args...]
Available arguments :

   -g, --gene           gene to scan for positions to cut on one strand only, ex : RYR1
   -c, --cas            cas9 you want to use to cut your gene (Only spcas9 available on v1.0)
   -f, --frequence      Minimal variant frequency in gnomAD v.3 population
   -o, --output,        Output file name to store results in
   -h, --help,          Show this help section

This pipeline was developped at the CHU Grenoble Alpes
Feel free to address any issue at : benoitclement.sand@gmail.com

#######################################################################################################
```

