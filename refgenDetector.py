from refgenDetector.reference_genome_dictionaries import GRCh38, GRCh37, t2t, set_hs37d5, set_hg16, set_hg17, set_hg18, set_hg19, \
    set_b37, verily_difGRCh38
import argparse
import logging
import sys
import csv
import gzip
import pysam


logger = logging.getLogger('reference_genome_logger')

def comparison(bam, file):
    """
    Compares the dictionaries and the length in the header. The dictionary that has more matches is the one used to align the bam
    """
    set_bamvalues = set(bam.values())

    if len(set_bamvalues.intersection(set_hs37d5)) != 0:
        print(f"{file},hs37d5")
    elif len(set_bamvalues.intersection(set_hg16)) != 0:
        print(f"{file},hg16")
    elif len(set_bamvalues.intersection(set_hg17)) != 0:
        print(f"{file},hg17")
    elif len(set_bamvalues.intersection(set_hg18)) != 0:
        print(f"{file},hg18")
    else:
        match = max(  # gets which reference genome has more matches with the file that's being analyzed
            (len(set_bamvalues.intersection(GRCh37.values())), 'GRCh37'),
            (len(set_bamvalues.intersection(GRCh38.values())), 'GRCh38/hg38'),
            (len(set_bamvalues.intersection(t2t.values())), 'T2T')
        )

        if match[1] == 'GRCh37':
            if len(set_bamvalues.intersection(
                    set_hg19)) != 0:  # if the length of the NC_001807.4 is present the reference genome used was hg19
                print(f"{file},hg19")
            elif len(set_bamvalues.intersection(set_b37)) != 0:
                print(f"{file},b37")
            else:
                print(f"{file},{match[1]}")
        elif match[1] == 'GRCh38/hg38':
            if len({key for key in bam.keys() if "HLA" in key}) != 0:
                print(f"{file}, hs38DH_extra")
            elif len(set_bamvalues.intersection(verily_difGRCh38.values())) != 0:
                print(f"{file}, Verily_GRCh38")
            else:
                print(f"{file},GRCh38/hg38")
        elif match == 0:
            if len({key for key in bam.keys() if "HLA" in key}) != 0: #check first HLA as it contains Verily's decoy contigs
                print(f"{file}, hs38DH_extra only with decoy contigs")
        else:
            print(f"{file}, no flavor can be inferred")


def get_info_txt(file_try, md5, assembly):
    """
    process_data takes the headers in the txt one by one and creates a dictionary with the SN and LN fields.
    If the AS and md5 field is present and asked by the user, it prints them.
    """
    header_reader = csv.reader(file_try, delimiter='\t')
    SQ = [line for line in header_reader if '@SQ' in line]  # creates a list with the SQ header lines
    bam = {line[1].replace('SN:', ''): int(line[2].replace('LN:', '')) for line in
           SQ}  # we convert the dictonary values to int, beacuse that's how we've constructed our collection of dictionaries
    comparison(bam, file_try.name)  # run the next function
    # if present and asked prints AS
    if assembly is not False:  # if the assembly argument is selected by the user
        ass = [l for line in SQ for l in line if 'AS' in l][:1]  # it saves the first AS field of the header
        if ass:  # if AS is present in the header
            print(f"*** {file_try.name}, {ass[0]}")  # prints the value
    # if present and asked prints md5
    if md5 is not False:  # if the md5 argument is selected by the user
        for i in SQ[0]:  # if the M4 field is in the header
            if 'M5' in i:
                m5 = {line[1].replace('SN:', ''): line[5].replace('M5:', '') for line in
                      SQ}  # creates a dictionary with the name of the contig and the md5 values found in the header
                print(f"*** {file_try.name}, MD5: {m5}")


def get_info_bamcram(header, file, md5, assembly):
    """Loop over the SQ (sequence dictionary) records in the header, creates a dictionary with the contigs names and lentghs,
    if present and asked for the user it prints AS and M5"""

    LN = {sq_record["SN"]: sq_record["LN"] for sq_record in
          header.get('SQ', [])}  # creates a dictionary with the name of the contigs and their length
    # assembly block
    if assembly is True:
        AS = set(sq_record["AS"] for sq_record in header.get('SQ', []) if
                 'AS' in sq_record)  # if the AS (Assembly sequence) field is present, it keeps record in a dictionary
        if AS:  # if AS was in the header
            print(f"*** {file}, AS:{AS.pop()}")
    # md5 block
    if md5 is True:
        if "M5" in header.get('SQ', [])[0]:
            M5 = {sq_record["SN"]: sq_record["M5"] for sq_record in header.get('SQ',
                                                                               [])}  # if the AS (Assembly sequence) field is present, it keeps record in a dictionary
            print(f"*** {file}, MD5: {M5}")
    comparison(LN, file)  # calls comparison () with the length values as a set


def process_data_bamcram(file, md5, assembly):
    if "bam" in file:  # if i's a BAM file
        save = pysam.set_verbosity(0)  # https://github.com/pysam-developers/pysam/issues/939
        f = pysam.AlignmentFile(file, "rb")  # open bam using pysam library
        pysam.set_verbosity(save)
    elif "cram" in file:  # if i's a CRAM file
        save = pysam.set_verbosity(0)  # https://github.com/pysam-developers/pysam/issues/939
        f = pysam.AlignmentFile(file, "rc")  # open cram using pysam library
        print(f"CRAM analyzed: {file}")
        pysam.set_verbosity(save)
    else:  # if the file is not a CRAM or a BAM
        print(f"The file provided, {file} is not a CRAM or a BAM and can't be analyzed\n")
        return
    header = f.header  # extract header object from AligmentFile class
    get_info_bamcram(header, file, md5, assembly)


def process_data_txt(file, md5, assembly):
    """
    try/excepts to make sure all the headers are read or that it skips a path when necessary
    """
    try:
        try:
            with gzip.open(file, 'rt') as file_try:  # try to open the file with utf-8 encoding
                get_info_txt(file_try, md5, assembly)
        except UnicodeError:
            with gzip.open(file, 'rt',
                           encoding='iso-8859-1') as file_except:  # try to open the file with iso-8859-1 enconding
                get_info_txt(file_except, md5, assembly)
    except UnicodeError:  # print ERROR if the encoding is unknown
        print("ERROR: The encoding of the file is not utf-8 or iso-8859-1")
        pass

    except OSError:
        try:
            with open(file, 'r') as file_try:  # if the txt with the header is not gz
                get_info_txt(file_try, md5, assembly)
        except UnicodeError:
            with open(file, 'r', encoding='iso-8859-1') as file_try:  # not ut8 encoding
                get_info_txt(file_try, md5, assembly)
        except UnicodeError:  # print ERROR if the encoding is unknown
            print("ERROR: The encoding of the file is not utf-8 or iso-8859-1")
            pass
        except OSError:  # if the path of the header is not correct
            print(f"{file} doesn't exist")


def main():
    """
    Parses the headers provided in a txt
    """
    try:
        parser = argparse.ArgumentParser(prog="INFERRING THE REFERENCE GENOME USED TO ALIGN BAM OR CRAM FILE")
        parser.add_argument("-p", "--path", help=" Path to txt with the files to analyze (one path per line)",
                            required=True)
        parser.add_argument("-t", "--type", choices=["BAM/CRAM", "Headers"], required=True)
        parser.add_argument("-m", "--md5", required=False, action='store_true',
                            help="Print contigs md5 value if present in the header")
        parser.add_argument("-a", "--assembly", required=False, action='store_true',
                            help="Print AS (assembly field) if present in the header")
        args = parser.parse_args()

        with open(args.path, 'r') as txt:  # reads the file provided
            for line in txt:  # for each file in the txt
                if args.type == "Headers":
                    process_data_txt(line.strip(), args.md5, args.assembly)
                else:
                    process_data_bamcram(line.strip(), args.md5, args.assembly)
    except Exception as e:
        logger.error("Error: {}".format(e))
        sys.exit(-1)



if __name__ == '__main__':  # the first executed function will be main()
    main()
