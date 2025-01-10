#!/usr/bin/env python

""" refgenDetector.py: Script to infer the reference genome used to create a BAM or CRAM"""

__author__ = "Mireia Marin Ginestar"
__version__ = "2.0"
__maintainer__ = "Mireia Marin Ginestar"
__email__ = "mireia.marin@crg.eu"
__status__ = "Developement"

import os
import sys
# Add the parent directory to the Python path
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.insert(0, parent_dir)
from refgenDetector.reference_genome_dictionaries import *
import argparse
import csv
import gzip
import pysam
import psutil
import time
from rich.console import Console

console = Console()

def monitor_resources(func):
    def wrapper(*args, **kwargs):
        process = psutil.Process()
        start_time = time.time()
        start_cpu_time = process.cpu_times()

        # Run the function
        result = func(*args, **kwargs)

        end_time = time.time()
        end_cpu_time = process.cpu_times()
        duration = end_time - start_time

        # Calculate CPU time
        cpu_user = end_cpu_time.user - start_cpu_time.user
        cpu_system = end_cpu_time.system - start_cpu_time.system
        total_cpu_time = cpu_user + cpu_system

        # Get memory usage (current RSS in MB)
        memory_usage = process.memory_info().rss / (1024 * 1024)  # Convert to MB

        # Get disk I/O
        io_counters = process.io_counters()
        bytes_read = io_counters.read_bytes
        bytes_written = io_counters.write_bytes

        # Print results
        print(f"Execution time: {duration:.2f} seconds")
        print(f"CPU time used: {total_cpu_time:.2f} seconds")
        print(f"Memory usage (RSS): {memory_usage:.2f} MB")
        if bytes_read > (1024*1024):
            print(f"Disk I/O - Read: {bytes_read / (1024 * 1024):.2f} MB, Written: {bytes_written / (1024 * 1024):.2f} MB")
        elif bytes_read > 1024:
            print(f"Disk I/O - Read: {bytes_read / (1024):.2f} KB, Written: {bytes_written / (1024):.2f} KB")
        else:
            print(f"Disk I/O - Read: {bytes_read:.2f} Bytes, Written: {bytes_written:.2f} Bytes")

        return result

    return wrapper

def intersection_targetfile_referencerepo(dict_SN_LN, reference_genome):
    """
    Find the matches between the target file and the repository of unique contigs per reference genome
    """
    return (len(set(dict_SN_LN.values()).intersection(reference_genome["ref_gen"].values())), reference_genome[
        "build"], reference_genome["species"])

def comparison(dict_SN_LN, target_file):
    """
    First, it defines the major release to which the header belongs to. Then, checks if there's any match with the flavors.
    """

    # List of matches for each major release
    matches = [
        intersection_targetfile_referencerepo(dict_SN_LN, major_releases[ref])
        for ref in major_releases
    ]

    # Find the major release with the maximum matches
    match = max(matches, key=lambda ref_gen_w_macthes: ref_gen_w_macthes[0]) # The key parameter specifies a function that extracts a
    # value from each element in the iterable to be used for comparisons

    if match[1] == "GRCh37": #check for GRCh37 flavors

        matches_flavors = [
            intersection_targetfile_referencerepo(dict_SN_LN, flavors_GRCh37[ref])
            for ref in flavors_GRCh37
        ]
        match_flavors = max(matches_flavors, key=lambda x: x[0])
        if match_flavors: #if some flavor was defined it prints it
            console.print(f"[bold]File:[/bold] {target_file} \n[bold]Specie detected:[/bold] {match_flavors[2]} "
            f"[bold]\nReference genome version:[/bold] {match_flavors[1]}")
        else: #if there wasnt any flavor inferred, the major release it printed
            console.print(f"[bold]File:[/bold] {target_file} \n[bold]Specie detected:[/bold] Homo sapiens \n["
                          f"bold]Reference genome version:[/bold] GRCh37")

    elif match[1] == "GRCh38": #checks for GRCh38 flavors

        if any("HLA-" in key for key in dict_SN_LN.keys()):
            #first checks if the contigs contain in their names HLA-
            console.print(f"[bold]File:[/bold] {target_file} \n[bold]Specie detected:[/bold] Homo sapiens \n[bold]"
                          f"Reference genome version:[/bold] hs38DH_extra") # if so, the reference genome used was
            # hs38DH_extra
        elif set(dict_SN_LN.values()).intersection(verily_difGRCh38.values()):#checks if the Verily's unique
            # lengths are present
            console.print(f"[bold]File:[/bold] {target_file} \n[bold]Specie detected:[/bold] Homo sapiens \n[bold]"
                          f"Reference genome version:[/bold] GRCh38_no_alt_plus_hs38d1")
        else: # if no GRCh38 flavor is inferred, the major release is printed
            console.print(f"[bold]File:[/bold] {target_file} \n[bold]Specie detected:[/bold] Homo sapiens \n["
                          f"bold]Reference genome version:[/bold] GRCh38)")
    else: # print the major releases with no considered flavors.
        console.print(f"[bold]File:[/bold] {target_file} \n[bold]Specie detected:[/bold] {match[2]} "
              f"\n[bold]Reference genome version:[/bold] {match[1]}")

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
            console.print(f"[bold]AS field:[/bold] {dict_assembly[0].split(':')[1]}")  # prints the value
    # if present and asked prints md5
    if md5:  # if the md5 argument is selected by the user
        for i in dict_SQ[0]: # checks in the first line if the M5 field is present
            if "M5" in i: # if it is (i = M5 field)
                dict_M5 = {line[1].replace("SN:", ""): i.replace("M5:", "") for line in
                      dict_SQ}  # creates a dictionary with the name of the contig and the md5 values found in the
                # header
                console.print(f"[bold]MD5 fields:[/bold] {dict_M5}")


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
            console.print(f"[bold]AS field:[/bold] {dict_assembly.pop()}")
    if md5: #if the user chose -m
        dict_M5 = set(sq_record["M5"] for sq_record in header_bam_cram.get("SQ", []) if
                 "M5" in sq_record)  # if the AS (Assembly sequence) field is present, it keeps record in a dictionary
        if dict_M5:
            console.print(f"[bold]M5 fields:[/bold]{dict_M5}") #prints the AS field just once
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
    #try: #if the file is indeed a txt
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
    """except: #if the file is not a txt it breaks
        print("Please, check you are using the correct --type or that the path to the header is correct")
        return # the txts in --path will be analyzed and the incorrect formats will be skipped"""

@monitor_resources
def main():
    """
    Process the users inputs and chooses to run BAM/CRAM module or txt module, depending on the -t argument
    """

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
            print("* Running refgenDetector v2.0.0 *")
            if args.type == "Headers":
                for target_file in txt:  # for each target file in the txt, it calls the function to open the
                    # headers saved in a txt and passes the arguments md5 and assembly.
                    console.print("[bold]---[/bold]")
                    process_data_txt(target_file.strip(), args.md5, args.assembly)
                    console.print("[bold]---[/bold]")
            else: # the target files will be BAMs or CRAMs
                for target_file in txt:  # for each target file in the txt, it calls the function to get headers
                    # from BAM and CRAMs and passes the arguments md5 and assembly.
                    console.print("[bold]---[/bold]")
                    process_data_bamcram(target_file.strip(), args.md5, args.assembly)
                    console.print("[bold]---[/bold]")
    except OSError: #if the file provided in --path cant be opened
        print("The file provided in --path doesn't exist. Make sure to include the complete path to a txt file "
              "formed by paths to headers saved in a txts or to BAM/CRAMs files (one per line)")



if __name__ == "__main__":  # the first executed function will be main()
    main()
