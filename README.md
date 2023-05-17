# EGA - RefgenDetector

RefgenDetector is a python tool that inferes the reference genome assembly used during the read alignment for BAM and CRAM files. The proposed tool is designed to facilitate the analysis of genomic data with incomplete metadata annotation. RefgenDetector can differentiate between major human reference genome releases, as well as commonly used flavors, by utilizing the LN and SN mandatory fields in the BAM and CRAM headers. The tool includes dictionaries with information on contig names and lengths, enabling it to accurately identify unique contigs and differentiate between different flavors of the reference genome.

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

## Installation 

TODO 

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

In the file (```-p argument```) path you should add the paths to all the files you want to analyze. RefgenDetector works with complete BAM and CRAMs and with txt files containing only the headers. The txt can be uncompressed, gzip compressed, and with encodings utf-8 and iso-8859-1. 

All the files included in this argument must be the same type, meaning, you should run RefgenDetector to analyze only BAM/CRAMs or only headers. 

## Licence and funding

RefgenDetector is released under GNU General Public License v3.0. 
It was funded by ELIXIR, the research infrastructure for life-science data (ELIXIR Beacon Implementation Studies 2019-2021 and 2022-2023).
