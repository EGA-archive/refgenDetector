#!/usr/bin/env python

""" refgenDetector.py: Script to infer the reference genome used to create a BAM or CRAM"""

__author__ = "Mireia Marin Ginestar"
__version__ = "0.1"
__maintainer__ = "Mireia Marin Ginestar"
__email__ = "mireia.marin@crg.eu"
__status__ = "Developement"

from reference_genome_dictionaries import GRCh38, GRCh37, hs37d5, hg16, hg17, hg18, hg19, \
    b37, verily_difGRCh38, T2T
import argparse
import logging
import sys
import csv
import gzip
import pysam

logger = logging.getLogger("reference_genome_logger")

def intersection_targetfile_referencerepo(dict_SN_LN, reference_genome):
    """
    Find the matches between the target file and the repository of unique contigs per reference genome
    """
    if len(set(dict_SN_LN.values()).intersection(reference_genome[0].values())) !=0:
        return reference_genome

def comparison(dict_SN_LN, target_file):
    """
    First, it defines the major release to which the header belongs to. Then, checks if theres any match with the
    flavors.
    """
    major_releases=[[hg16, "hg16"], [hg17, "hg17"], [hg18,"hg18"], [GRCh37, "GRCh37"], [GRCh38,"GRCh38"], [T2T,
                                                                                                             "T2T"]]
    # list of the major releases that the script can identify
    flavors=[[hs37d5, "hs37d5"], [b37,"b37"], [hg19, "hg19"]]
    # list of the GRCh37 flavors that the script caon identify
    major_release_list = []
    for refgen in major_releases: #infer the major release family
        major_release_list.append(intersection_targetfile_referencerepo(dict_SN_LN,refgen))
    family = [ref for ref in major_release_list if ref is not None] # deletes all the Nones and gets the family of
    # the header
    if len(family)==0: #if there wasnt any positive results, no reference genome can be inferred
        print("The reference genome can't be inferred")
    elif family[0][1] == "GRCh37": #if the header belongs to GRCh37 family
        for flav in flavors: #checks GRCH37 flavors
            match_flavor = intersection_targetfile_referencerepo(dict_SN_LN, flav)
            if match_flavor:  # if a GRCh37 flavor is detected
                break  # the loop stops. If it continues all hs37d5 will also be identified as b37,
                # as both contain NC_007605
        if match_flavor:  # if a flavor is inferred, print it
            print(f"{target_file}, {match_flavor[1]}")
        else:  # if there arent matches with the flavors, prints the major release
            print(f"{target_file},{family[0][1]}")
    elif family[0][1] =="GRCh38": # if the header belongs to GRCh38 family
        if len({key for key in dict_SN_LN.keys() if "HLA-" in key}) != 0:  # check the SN field (contigs name)
            # match the HLA nomenclature
            print(f"{target_file}, hs38DH_extra")
        elif intersection_targetfile_referencerepo(dict_SN_LN, [verily_difGRCh38, "Verily's GRCh38"]):
            # checks any of the Verily's contigs is present
            print(f"{target_file}, Verily's GRCh38")
        else:  # if there arent unique contig, the major release is printed
            print(f"{target_file}, {family[1]}")
    else: #Prints the family (with no flavors) the headers belong to
        print(f"{target_file},{family[0][1]}")

def get_info_txt(header_txt, md5, assembly):
    """
    Second function of the txt module. Extracts the SQ (sequence dictionary) records in the header, creates a
    dictionary with the contigs names and lengths, and, if present and requested by the user (adding -m and -a in the
    argument) prints AS and M5
    """
    header_reader = csv.reader(header_txt, delimiter="\t")
    dict_SQ = [line for line in header_reader if "@SQ" in line]  # creates a list with the SQ header
    # lines
    try:
        dict_SN_LN = {line[1].replace("SN:", ""): int(line[2].replace("LN:", "")) for line in
           dict_SQ}  #the dictonary values must be int due to the structure of the collection of reference dictionaries
    except ValueError:
        print(f"Check the LN field of your header {header_txt.name} only contains numbers")
        return
    comparison(dict_SN_LN, header_txt.name)  # run the next function
    # if present and asked by the user prints AS
    if assembly:  # if the assembly argument is selected by the user
        dict_assembly = [l for line in dict_SQ for l in line if "AS" in l][:1]  # it saves the first AS field of the
        # header
        if dict_assembly:  # if AS is present in the header
            print(f"*** {header_txt.name}, {header_txt[0]}")  # prints the value
    # if present and asked prints md5
    if md5:  # if the md5 argument is selected by the user
        for i in dict_SQ[0]: # checks in the first line if the M5 field is present
            if "M5" in i: # if it is (i = M5 field)
                dict_M5 = {line[1].replace("SN:", ""): i.replace("M5:", "") for line in
                      dict_SQ}  # creates a dictionary with the name of the contig and the md5 values found in the
                # header
                print(f"*** {header_txt.name}, MD5: {dict_M5}")


def get_info_bamcram(header_bam_cram, target_file, md5, assembly):
    """
    Second function of the BAM/CRAM module. Loop over the SQ (sequence dictionary) records in the header, creates a
    dictionary with the contigs names and lengths, if present and requested by the user (adding -m and -a in the argument) prints AS and M5
    """

    dict_SN_LN = {sq_record["SN"]: sq_record["LN"] for sq_record in
          header_bam_cram.get("SQ", [])}  # creates a dictionary with the name of the contigs and their length
    if assembly:
        dict_assembly = set(sq_record["AS"] for sq_record in header_bam_cram.get("SQ", []) if
                 "AS" in sq_record)  # if the AS (Assembly sequence) field is present, it keeps record in a dictionary
        if dict_assembly:  # if AS was in the header
            print(f"*** {target_file}, AS:{dict_assembly.pop()}")
    if md5: #if the user chose -m
        dict_M5 = set(sq_record["M5"] for sq_record in header_bam_cram.get("SQ", []) if
                 "M5" in sq_record)  # if the AS (Assembly sequence) field is present, it keeps record in a dictionary
        if dict_M5:
            print(f"*** {target_file}, M5:{dict_M5}") #prints the AS field just once
    comparison(dict_SN_LN, target_file)  # calls comparison () with the length values as a set


def process_data_bamcram(target_file, md5, assembly):
    """
    First function of the BAM/CRAM module. It opens each BAM or CRAM provided by the user and extracts the header.
    """
    try:
        save = pysam.set_verbosity(0)  # https://github.com/pysam-developers/pysam/issues/939
        bam_cram = pysam.AlignmentFile(target_file, "rb")  # open bam/cram using pysam library
        pysam.set_verbosity(save)
    except Exception: # printed if the user chose -t BAM/CRAM but the paths in -p were pointing to txts
        print("The BAM/CRAMs can't be opened, please check your path or that you are using the correct --type")
        return #the bam and cram in --path will be analyzed and the incorrect format will be skipped
    header_bam_cram = bam_cram.header  # extract header object from AligmentFile class
    get_info_bamcram(header_bam_cram, target_file, md5, assembly)



def process_data_txt(target_file, md5, assembly):
    """
    First function of the txt module. It opens each header (saved in a txt) provided by the user. Its prepared to open a
    txt compressed with gzip or uncompressed. It can read both utf-8 and iso-8859-1.
    """
    try: #if the file is indeed a txt
        try: #tries to open an uncompressed txt
            try: #tries to open the file with utf-8 encoding
                with open(target_file,"r") as header_txt:
                    get_info_txt(header_txt, md5, assembly)
            except UnicodeError: #tries to open the file with iso-8859-1 encoding
                with open(target_file,"r", encoding="iso-8859-1") as header_txt:  # tries to open the file with utf-8
                    # encoding
                    get_info_txt(header_txt, md5, assembly)
        except: #tries to open a compressed txt
            try: #tries to open a compressed file with utf-8 encoding
                with gzip.open(target_file,"rt") as header_txt:
                    get_info_txt(header_txt, md5, assembly)
            except: #tries to open a compressed file with iso-8859-1 encoding
                with gzip.open(target_file,"rt", encoding="iso-8859-1") as header_txt:  # tries to open the file with
                    # utf-8 encoding
                    get_info_txt(header_txt, md5, assembly)
    except: #if the file is not a txt it breaks
        print("Please, check you are using the correct --type")
        return # the txts in --path will be analyzed and the incorrect formats will be skipped

def main():
    """
    Process the users inputs and chooses to run BAM/CRAM module or txt module, depending on the -t argument
    """
    try:
        parser = argparse.ArgumentParser(prog="INFERRING THE REFERENCE GENOME USED TO ALIGN BAM OR CRAM FILE")
        #MANDATORY ARGUMENTS
        parser.add_argument("-p", "--path", help="Path to main txt. It will consist of the paths to the files to be "
                                                 "analyzed (one path per line)",
                            required=True)
        parser.add_argument("-t", "--type", choices=["BAM/CRAM", "Headers"], help="All the files in the txt provided "
                                                                                  "in --path must be BAM/CRAMs or "
                                                                                  "headers in a txt. Choose -t"
                                                                                  "depending on the type of files you are going to "
                                                                                  "analyze",
                                                                              required=True)
        #OPTIONAL ARGUMENTS
        parser.add_argument("-m", "--md5", required=False, action="store_true",
                            help="[OPTIONAL] If you want to obtain the md5 of the contigs present in the header, "
                                 "add --md5 to "
                                 "your command. This will print the md5 values if the field M5 was present in "
                                 "your header")
        parser.add_argument("-a", "--assembly", required=False, action="store_true",
                            help="[OPTIONAL] If you want to obtain the assembly declared in the header add --assembly "
                                 "to "
                                 "your command. This will print the assembly if the field AS was present in "
                                 "your header")
        args = parser.parse_args()
        try: #try to open the main txt (-p)
            with open(args.path,"r") as txt:  # reads the txt with the paths to analyze
                if args.type == "Headers":
                    for target_file in txt:  # for each target file in the txt, it calls the function to open the
                        # headers saved in a txt and passes the arguments md5 and assembly.
                        process_data_txt(target_file.strip(), args.md5, args.assembly)
                else: # the target files will be BAMs or CRAMs
                    for target_file in txt:  # for each target file in the txt, it calls the function to get headers
                        # from BAM and CRAMs and passes the arguments md5 and assembly.
                        process_data_bamcram(target_file.strip(), args.md5, args.assembly)
        except OSError: #if the file provided in --path cant be opened
            print("The file provided in --path doesn't exist. Make sure to include the complete path to a txt file "
                  "formed by paths to headers saved in a txts or to BAM/CRAMs files (one per line)")
    except Exception as e:
        logger.error("Error: {}".format(e))
        sys.exit(-1)



if __name__ == "__main__":  # the first executed function will be main()
    main()

