# EGA - RefgenDetector

RefgenDetector is a bioinformatics tool that **infers the reference genome assembly** used to create aligment files (BAM/CRAM/header) and VCFs. 

## Aligment Files

It identifies major genome releases and derived assemblies across humans and multiple other species by analyzing contig names and lengths **from the header**. Benchmarking against 94 synthetic datasets achieved a 100% accuracy rate, while large-scale testing on 918,404 real-world files demonstrated 97.13% correctness, failing only when files’ headers are incomplete.

### Description

RefgenDetector is able to infer the following reference genomes:

**Primates**

👤 Homo sapiens

- hg16
- hg17
- hg18
- GRCh37
- GRCh38
- T2T

🐒 Pan troglodytes

- pantro3_0
- Pan_troglodytes-2.1

🐵 Macaca mulatta

- Mmul10
- rheMac8
- rheMac3

**Rodents**

🐭 Mus musculus

- mm7
- mm8
- mm9
- mm10
- mm39

🐀 Rattus norvegicus

- mRatBN7_2
- Rnor_6_0

**Other Mammals**

🐷 Sus scrofa

- Sscrofa10_2
- Sscrofa11_1

**Vertebrates (Non-Mammalian)**

🐟 Danio Rerio

- danRer10
- danRer11

**Invertebrates**

🪰 Drosophila Melanogaster

- dm5
- dm6

🐛 Caenorhabditis elegans

- WBcel215
- WBcel235

**Microorganisms & Plants**

🧫 Escherichia coli

- ASM886v2
- ASM584v2

🌱 Arabidopsis thaliana

- TAIR

🍺 Saccharomyces cerevisiae

- R64

## Variant Calling Files (VCFs)

From VCF files only 4 human assemblies can be inferred: 

- Hg18
- GRCh37
- GRCh38
- T2T

Two different sources of information are used to infer the reference genome from variant calling files

* **Header**
  
In the VCF specification it is recommended, but **not mandatory** that the VCF header includes tags describing the reference and contigs backing the data contained in the file. When present, the tool will analyze this information and output the reference genome version based on the contig lengths, following the same logic of the aligment files inference.

* **Variants**

To infer the reference genome from a VCF the tool will read the VCF file in chunks of 100.000 variants, avoiding to load the complete file in memory. The `POS` and `REF` columns will be extracted and compared to the pkl files.

The pkl files were created comparing the nucleotides in each position for hg18, GRCh37, GRCh38 and T2T. Each file contains a list of the positions where each reference had a different nucleotide (distinguishing positions). 

By getting the number of matches between these distinguishing positions and the `REF` present in the VCF we infer the reference genome version used to call the variants. 

## Requirements

- Python 3.10.6

Depending on how you want to install the package:

- pip
- Docker

Download the PKL files for the inference with VCFs: 

1. [Download the pkl reference](https://crgcnag-my.sharepoint.com/:u:/g/personal/mimarin_crg_es/IQByWKuqxkRMR7IfzwTy6CFiAUdAcXVkFH7IYQMj8wwrYTs?e=rt8p19)

2. Move the pkls to the correct path:

```
mv pkls.zip /refgenDetector/src/refgenDetector/
unzip /refgenDetector/src/refgenDetector/pkls.zip
```

## Installation

### Cloning this repository

1. Clone this repository

2. ``` $ cd PATH_WHERE_YOU_CLONED_THE_REPOSITORY/src/refgenDetector ```

3. ``$ python3 refgenDetector_main.py -h ``

### From pypi

``$ pip install refgenDetector``

## Usage

You can get the help menu by running:

```
$ refgenDetector -h
```

```
usage: INFERRING THE REFERENCE GENOME USED TO ALIGN BAM OR CRAM FILE [-h] -f FILE -t {BAM/CRAM,Header,VCF,BIM} [--md5] [-a] [-v MAX_N_VAR] [-m MATCHES] [-r]

optional arguments:
  -h, --help            show this help message and exit
  -f FILE, --file FILE  Input file path
  -t {BAM/CRAM,Header,VCF,BIM}, --type {BAM/CRAM,Header,VCF,BIM}
                        Type of files to analyze.
  --md5                 Print md5 values if present in header.
  -a, --assembly        Print assembly if present in header.
  -v MAX_N_VAR, --max_n_var MAX_N_VAR
                        Maximum number of variants to read before stopping inference. The file is processed in chunks of 100,000 variants, so this value must be a multiple of 100,000 (e.g. 100000,
                        200000, 300000, ...).
  -m MATCHES, --matches MATCHES
                        Number of matches required before stopping. [DEFAULT:5000]
  -r, --resources       When set, print execution time, CPU, memory, and disk I/O usage
```

## Test RefgenDetector

In the folder **examples** you can find headers, BAM and CRAMs to test the working of RefgenDetector.

*All this files belong to the [synthetics data cohort](https://ega-archive.org/synthetic-data) from the European
Genome-Phenome Archive ([EGA](https://ega-archive.org/)).*

### Test with headers in a TXT

In the folder TEST_HEADERS there are four headers obtained from synthetic BAM an CRAMs stored in the EGA. Each one of
them belongs to a different synthetic study:

- Test Study for EGA using data from 1000 Genomes Project - Phase
  3 [EGAS00001005042](https://ega-archive.org/studies/EGAS00001005042).
- Synthetic data - Genome in a Bottle - [EGAS00001005591](https://ega-archive.org/studies/EGAS00001005591).
- Human genomic and phenotypic synthetic data for the study of rare
  diseases - [EGAS00001005702](https://ega-archive.org/studies/EGAS00001005702).
- CINECA synthetic data.Please note: This study contains synthetic data (with cohort “participants” / ”subjects” marked
  with FAKE) has no identifiable data and cannot be used to make any inference about cohort data or
  results - [EGAS00001002472](https://ega-archive.org/studies/EGAS00001002472).

Further information about them can be found in the file *where_to_find_this_files.txt*, saved in the same folder.

To run RefgenDetector with the files:

1. Modify the txt *path_to_headers* so the paths match those in your computer.
2. Run:

   ``` $ refgenDetector -p /PATH_WHERE_YOU_CLONED_THE_REPOSITORY/refgenDetector/examples/path_to_headers -t Headers```

### Test with BAM and CRAMs

In the folder TEST_BAM_CRAM there are a BAM and a CRAM obtained from synthetic BAM an CRAMs stored in the EGA. They
belong to the synthetic study - Test Study for EGA using data from 1000 Genomes Project - Phase
3 [EGAS00001005042](https://ega-archive.org/studies/EGAS00001005042).

Further information about them can be found in the file *where_to_find_this_files.txt*, saved in the same folder.

To run RefgenDetector with the files:

1. Modify the txt *path_to_bam_cram* so the paths match those in your computer.

2. Run:

   ``` $ refgenDetector -p /PATH_WHERE_YOU_CLONED_THE_REPOSITORY/refgenDetector_pip-master/examples/path_to_bam_cram -t BAM/CRAM```



## Licence and funding

RefgenDetector is released under GNU General Public License v3.0.

It was funded by ELIXIR, the research infrastructure for life-science data (ELIXIR Beacon Implementation Studies
2019-2021 and 2022-2023).
