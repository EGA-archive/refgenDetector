# Examples to test RefgenDetector

## Test with headers

In the folder Test_headers there are four headers obtained from synthetic BAM an CRAMs stored in the European 
Genome-Phenome Archive (EGA). 

Further information about them can be found in the file *where_to_find_this_files.txt*, saved in the same folder.

To run RefgenDetector with the files: 

1. Modify the txt *path_to_headers* so the paths match those in your computer. 
2. Run:
´´$ refgenDetector -p /PATH_WHERE_YOU_CLONED_THE_REPOSITORY/refgenDetector/examples/path_to_headers -t Headers´´

Check your installation has been successful by checking the test results are correct:

PATH_TO_YOUR_COMPUTER_SETUP/refgenDetector_pip-master/examples/TEST_HEADERS/EGAF00001753746, b37 PATH_TO_YOUR_COMPUTER_SETUP/refgenDetector_pip-master/examples/TEST_HEADERS/EGAF00005469864, hg19 PATH_TO_YOUR_COMPUTER_SETUP/refgenDetector_pip-master/examples/TEST_HEADERS/EGAF00005572695.gz, hs37d5 PATH_TO_YOUR_COMPUTER_SETUP/refgenDetector_pip-master/examples/TEST_HEADERS/EGAF00007462306, hs38DH_extra


## Test with BAM and CRAMs

In the folder Test_bam_cram there are a BAM and a CRAM obtained from the synthetic data stored in the 
European 
Genome-Phenome Archive (EGA). 

Further information about them can be found in the file *where_to_find_this_files.txt*, saved in the same folder.

To run RefgenDetector with the files: 

1. Modify the txt *path_to_bam_cram* so the paths match those in your computer. 
2. Run:
``
TODO
``
