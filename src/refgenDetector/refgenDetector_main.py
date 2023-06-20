#!/usr/bin/env python

""" refgenDetector.py: Script to infer the reference genome used to create a BAM or CRAM"""

__author__ = "Mireia Marin Ginestar"
__version__ = "0.1"
__maintainer__ = "Mireia Marin Ginestar"
__email__ = "mireia.marin@crg.eu"
__status__ = "Developement"

import os
import sys
# Add the parent directory to the Python path
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.insert(0, parent_dir)
from refgenDetector.reference_genome_dictionaries import GRCh38, GRCh37, hs37d5, hg16, hg17, hg18, hg19, \
    b37, verily_difGRCh38, T2T
import argparse
import logging
import csv
import gzip
import pysam

logger = logging.getLogger("reference_genome_logger")

def intersection_targetfile_referencerepo(dict_SN_LN, reference_genome):
    """
    Find the matches between the target file and the repository of unique contigs per reference genome
    """
    return set(dict_SN_LN.values()).intersection(reference_genome.values())

def comparison(dict_SN_LN, target_file):
    """
    First, it defines the major release to which the header belongs to. Then, checks if there's any match with the flavors.
    """
    major_releases = {"hg16": hg16,"hg17": hg17,"hg18": hg18,"GRCh37": GRCh37,"GRCh38": GRCh38,"T2T": T2T} #major
    # release that can be inferred
    flavors_GRCh37 = {"hs37d5": hs37d5,"b37": b37,"hg19": hg19} #GRCh37 flavors that can be inferred
    major_release_list = [major_releases[ref] for ref in major_releases if
                          intersection_targetfile_referencerepo(dict_SN_LN, major_releases[ref])] #gets the major
    # release to which the header belongs to
    if len(major_release_list) == 0: # if there wasnt any major release inferred print:
        print("The reference genome can't be inferred")
    else:
        if major_release_list[0] == major_releases["GRCh37"]: #check for GRCh37 flavors
            match_flavor = next(
                (flavors_GRCh37[flav] for flav in flavors_GRCh37 if intersection_targetfile_referencerepo(dict_SN_LN,
                                                                                             flavors_GRCh37[flav])),
                None) #infers the flavor of GRCh37 to which the header belongs to
            if match_flavor: #if some flavor was defined it prints it
                print(f"{target_file}, {[k for k, v in flavors_GRCh37.items() if v == match_flavor][0]}")
            else: #if there wasnt any flavor inferred, the major release it printed
                print(f"{target_file}, GRCh37")
        elif major_release_list[0] == major_releases["GRCh38"]: #checks for GRCh38 flavors
            if any("HLA-" in key for key in dict_SN_LN.keys()): #first checks if the contigs contain in their names HLA-
                print(f"{target_file}, hs38DH_extra") #if so, the reference genome used was hs38DH_extra
            elif intersection_targetfile_referencerepo(dict_SN_LN, verily_difGRCh38):#checks if the Verily's unique
                # lengths are present
                print(f"{target_file}, Verily's GRCh38")
            else: # if no GRCh38 flavor is inferred, the major release is printed
                print(f"{target_file}, GRCh38")
        else: #print the major releases with no considered flavors.
            print(f"{target_file}, {[k for k, v in major_releases.items() if v == major_release_list[0]][0]}")

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
        print("Please, check you are using the correct --type or that the path to the header is correct")
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
