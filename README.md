# EGA - RefgenDetector

RefgenDetector is a python tool that infers the reference genome assembly used during the read alignment for BAM and
CRAM files. The proposed tool is designed to facilitate the analysis of genomic data with incomplete metadata
annotation. RefgenDetector can differentiate between major human reference genome releases, as well as commonly used
flavors, by utilizing the LN and SN mandatory fields in the BAM and CRAM headers. The tool includes dictionaries with
information on contig names and lengths, enabling it to accurately identify unique contigs and differentiate between
different flavors of the reference genome.

## Description

RefgenDetector is able to infer the following reference genomes:

- NCBI35/hg17
- NCBI36.1/hg18
- GRCh37
- hg19
- b37
- hs37d5
- GRCh38
- Verily's GRCh38
- hs38DH_extra
- T2T

## Requirements

- Python 3.10.6

Depending on how you want to install the package:

- pip 23.1.2
- Docker version 24.0.2

## Installation

### Cloning this repository

1. Clone this repository

2. ``` $ cd PATH_WHERE_YOU_CLONED_THE_REPOSITORY/src/refgenDetector ```

3. ``$ python3 refgenDetector.py -h ``

### From pypi

All the information can be found [here](https://pypi.org/project/refgenDetector/)

``$ pip install refgenDetector``

### Docker

All the instructions to run refgenDetector in a Docker container can be found [here](https://hub.docker.com/r/beacon2ri/refgendetector). 

## Usage

You can get the help menu by running:

```
$ refgenDetector -h
```

```
usage: INFERRING THE REFERENCE GENOME USED TO ALIGN BAM OR CRAM FILE
       [-h] -p PATH -t {BAM/CRAM,Headers} [-m] [-a]

options:
  -h, --help            show this help message and exit
  -p PATH, --path PATH  Path to main txt. It will consist of the paths to the
                        files to be analyzed (one path per line)
  -t {BAM/CRAM,Headers}, --type {BAM/CRAM,Headers}
                        All the files in the txt provided in --path must be
                        BAM/CRAMs or headers in a txt. Choose -tdepending on
                        the type of files you are going to analyze
  -m, --md5             [OPTIONAL] If you want to obtain the md5 of the
                        contigs present in the header, add --md5 to your
                        command. This will print the md5 values if the field
                        M5 was present in your header
  -a, --assembly        [OPTIONAL] If you want to obtain the assembly declared
                        in the header add --assembly to your command. This
                        will print the assembly if the field AS was present in
                        your header

```

In the main file (```-p argument```) you should add the paths to all the files you want to analyze. RefgenDetector
works with complete BAM and CRAMs and with txt files containing only the headers. The txt can be uncompressed, gzip
compressed, and with encodings utf-8 and iso-8859-1.

All the files included in this argument must be the same type, meaning, you should run RefgenDetector to analyze only
BAM/CRAMs or only headers.

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

   Check your installation has been successful by checking the test results are correct:

   PATH_TO_YOUR_COMPUTER_SETUP/refgenDetector_pip-master/examples/TEST_HEADERS/EGAF00001753746, b37
   PATH_TO_YOUR_COMPUTER_SETUP/refgenDetector_pip-master/examples/TEST_HEADERS/EGAF00005469864, hg19
   PATH_TO_YOUR_COMPUTER_SETUP/refgenDetector_pip-master/examples/TEST_HEADERS/EGAF00005572695.gz, hs37d5
   PATH_TO_YOUR_COMPUTER_SETUP/refgenDetector_pip-master/examples/TEST_HEADERS/EGAF00007462306, hs38DH_extra

### Test with BAM and CRAMs

In the folder TEST_BAM_CRAM there are a BAM and a CRAM obtained from synthetic BAM an CRAMs stored in the EGA. They
belong to the synthetic study - Test Study for EGA using data from 1000 Genomes Project - Phase
3 [EGAS00001005042](https://ega-archive.org/studies/EGAS00001005042).

Further information about them can be found in the file *where_to_find_this_files.txt*, saved in the same folder.

To run RefgenDetector with the files:

1. Modify the txt *path_to_bam_cram* so the paths match those in your computer.

2. Run:
`` $ refgenDetector -p /PATH_WHERE_YOU_CLONED_THE_REPOSITORY/refgenDetector_pip-master/examples/path_to_bam_cram-t BAM/CRAM ``

Check your installation has been successful by checking the test results are correct:

PATH_TO_YOUR_COMPUTER_SETUP/refgenDetector_pip-master/examples/TEST_BAM_CRAM/HG00096.GRCh38DH__1097r__10.10000-10100__21.5000000-5050000.bam, hs38DH_extra
PATH_TO_YOUR_COMPUTER_SETUP/refgenDetector_pip-master/examples/TEST_BAM_CRAM/HG00096.GRCh38DH__1097r__10.10000-10100__21.5000000-5050000.cram, hs38DH_extra

## Licence and funding

RefgenDetector is released under GNU General Public License v3.0.

It was funded by ELIXIR, the research infrastructure for life-science data (ELIXIR Beacon Implementation Studies
2019-2021 and 2022-2023).
