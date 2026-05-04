# Examples to test RefgenDetector

*All this files belong to the [synthetics data cohort](https://ega-archive.org/synthetic-data) from the European
Genome-Phenome Archive ([EGA](https://ega-archive.org/)).*

These links describe the process to download EGA data: 

- [Request Data Access](https://ega-archive.org/access/request-data/how-to-request-data/)

- [Download method 1: pyEGA](https://ega-archive.org/access/download/files/pyega3/)

- [Download method 2: Live Distribution](https://ega-archive.org/access/download/files/live-outbox/)

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

## Test with headers

In the folder TEST_HEADERS there are four headers obtained from synthetic BAM an CRAMs stored in the European 
Genome-Phenome Archive (EGA). 

Further information about them can be found in the file *where_to_find_this_files.txt*, saved in the same folder.

Example query: 

```
refgenDetector -t Header -f examples/TEST_HEADERS/EGAF00007462306
```

This is the results you should expect:

```
* Running refgenDetector 3.0.0 *
---
File: examples/TEST_HEADERS/EGAF00001753746 
Specie detected: Homo sapiens 
Reference genome version: b37
---
File: examples/TEST_HEADERS/EGAF00005469864 
Specie detected: Homo sapiens 
Reference genome version: hg19
---
File: examples/TEST_HEADERS/EGAF00005572695.gz 
Specie detected: Homo sapiens 
Reference genome version: hs37d5
---
File: examples/TEST_HEADERS/EGAF00007462306 
Specie detected: Homo sapiens 
Reference genome version: hs38DH_extra
---
File: examples/TEST_HEADERS/header_ASM58v2 
Specie detected: Escherichia Coli 
Reference genome version: ASM584v2
---
File: examples/TEST_HEADERS/header_danrer11 
Specie detected: Danio Rerio 
Reference genome version: danRer11
---
File: examples/TEST_HEADERS/header_dm6 
Specie detected: Drosophila Melanogaster 
Reference genome version: dm6
---
File: examples/TEST_HEADERS/header_MFA1912RKS 
Specie detected: Homo sapiens 
Reference genome version: hg16
---
File: examples/TEST_HEADERS/header_mm10 
Specie detected: Mus musculus 
Reference genome version: mm10
---
File: examples/TEST_HEADERS/header_pantro2.1.4 
Specie detected: Pan troglodytes 
Reference genome version: Pan_troglodytes-2.1
---
File: examples/TEST_HEADERS/header_R64 
Specie detected: Saccharomyces cerevisiae 
Reference genome version: R64
---
File: examples/TEST_HEADERS/header_rhemac3 
Specie detected: Macaca mulatta 
Reference genome version: rheMac3
---
File: examples/TEST_HEADERS/header_rnor_6 
Specie detected: Rattus norvegicus 
Reference genome version: Rnor_6_0
---
File: examples/TEST_HEADERS/header_tair 
Specie detected: Arabidopsis thaliana 
Reference genome version: TAIR
---
File: examples/TEST_HEADERS/header_WBcel235 
Specie detected: Caenorhabditis elegans 
Reference genome version: WBcel235
---


```

## Test with BAM and CRAMs

In the folder TEST_BAM_CRAM there are a BAM and a CRAM obtained from the synthetic data stored in the 
European 
Genome-Phenome Archive (EGA). 

Further information about them can be found in the file *where_to_find_this_files.txt*, saved in the same folder.

Example query: 

```
refgenDetector -t Header -f examples/TEST_HEADERS/EGAF00007462306
```

This is the results you should expect:


```
* Running refgenDetector 3.0.0 *
---
File: /home/mireia/Bioinfo/BioQC/RefgenDetector/refgenDetector/examples/TEST_BAM_CRAM/HG00096.GRCh38DH__1097r__10.10000-10100__21.5000000-5050000.bam 
Specie detected: Homo sapiens 
Reference genome version: hs38DH_extra
---
File: /home/mireia/Bioinfo/BioQC/RefgenDetector/refgenDetector/examples/TEST_BAM_CRAM/HG00096.GRCh38DH__1097r__10.10000-10100__21.5000000-5050000.cram 
Specie detected: Homo sapiens 
Reference genome version: hs38DH_extra
---
```

## Test with VCFs

Create an EGA account and download [EGAF00007553559](https://ega-archive.org/datasets/EGAD00001003338)

Query: 

```
refgenDetector -t VCF -f EGAF00007553559/EE-2564.ALL_2504_Samples.wgs.phase3.v5.20130502.genotypes_only_chr1and2.vcf.gz
```

This is the output you should get:

```
* Running refgenDetector v.3.0.4 *
---
++ INFORMATION INFERRED BY THE HEADER ++

File: EGAF00007553559/EE-2564.ALL_2504_Samples.wgs.phase3.v5.20130502.genotypes_only_chr1and2.vcf.gz
File type: VCF
Species detected: Homo sapiens 
Reference genome version  : hs37d5
Number of samples: 2504

++ INFORMATION INFERRED BY THE REF COLUMN ++ 
 
Variants being mapped from: chr1
Matches: 
{'hg18': 0, 'GRCh37': 9057, 'GRCh38': 0, 'T2T': 0}
Inferred Reference genome: GRCh37
---
```