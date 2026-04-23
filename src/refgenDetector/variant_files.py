import gzip
import sys
import time
import pandas as pd
from dns.inet import inet_pton
from rich.console import Console
import os
import numpy as np 
import msgpack
try:
    # Works when installed as a pip package
    from .aligment_files import *
    from .chromosomes_dict import *
except ImportError:
    # Works when run directly as a script
    from aligment_files import comparison
    from chromosomes_dict import *



final_results = []
console = Console(highlight=False)
_msgpack_cache = {}
MSGPACK_DIR = "./msgpacks"

def gather_and_sum(lists):
    """
    It gathers and sums all the matches calculated in get_matches()
    Args:
        lists: list of lists with the matches calculated in get_matches() for each chunk.
    Returns:
        A dictionary with the total number of matches for each version, which is used to infer the reference genome version.
    """
    cumulative_sums = {}
    for lst in lists:
        for key, value in lst:
            if key in cumulative_sums:
                cumulative_sums[key] += value
            else:
                cumulative_sums[key] = value
    console.print(f"[bold]Matches: [/bold]", cumulative_sums)
    return cumulative_sums

def get_matches(snps_dict, chr_):
    """
    For each chromosome, it looks for the matches of the SNPs in the corresponding msgpack files, 
    which contain the reference genome information for each version. It returns a list with the number of matches for each version.
    Args:        
        snps_dict: dictionary with the position and the nucleotide of the SNPs in the chunk being processed. 
        chr_: chromosome of the chunk being processed, used to load the corresponding msgpack and get the matches.
    Returns:
        A list with the number of matches for each version, which is used to infer the reference genome version.
    """
    
    global _msgpack_cache
    start = time.time()

    genome_versions = ["hg18", "GRCh37", "GRCh38", "T2T"]
    matches = []

    # Extract arrays once
    positions = np.array(list(snps_dict.keys()), dtype=np.int64)
    nucleotides = np.array(list(snps_dict.values()))

    for version_name in genome_versions:
        cache_key = f"{version_name}-{chr_}"

        if cache_key not in _msgpack_cache:
            path = f"{MSGPACK_DIR}/{cache_key}.msgpack"
            if not os.path.exists(path):
                continue
            with open(path, "rb") as f:
                _msgpack_cache[cache_key] = msgpack.load(f, raw=False, strict_map_key=False)

        ref_dict = _msgpack_cache[cache_key]

        # Vectorized: look up all positions at once, compare arrays
        ref_nucs = np.array([ref_dict.get(p, None) for p in positions])
        match_count = int(np.sum(ref_nucs == nucleotides))

        matches.append([version_name, match_count])

    console.print("Getting matches. Took:", time.time() - start, "s")
    return matches


def trimming_indels(content, ref):
    """
    If a row is longer than one position it is deleted, deleting this way any indels
    Args:
        content: chunk of the file being read, containing the variants information. It must only contain SNPs, as indels are trimmed in this function.
        ref: column number of the reference allele in the file, which is different in vcfs and bim files. 
    Returns:
        A dataframe with only the SNPs in the chunk, and the gVCF status of the chunk.
    """
    try:
         # Keep only rows where REF is length 1 or <NON_REF>
        del_insertions = content[(content.iloc[:, 2].str.len() == 1) | (content.iloc[:, 2] == '<NON_REF>')]
        
        # Keep only rows where ALT is length 1 or <NON_REF>
        vcf_snps = del_insertions[
            (del_insertions.iloc[:, 3].str.len() == 1) | (del_insertions.iloc[:, 3] == '<NON_REF>')]
        
        # Select columns and make an explicit copy
        snps = vcf_snps.iloc[:, [1, ref]].copy()
        snps.columns = ['position', 'nucleotide']

        # Safe assignment using .loc
        snps.loc[:, 'nucleotide'] = snps['nucleotide'].str.upper()

        gVCF = False
        if content.iloc[:, 3].astype(str).str.contains('<NON_REF>').any():
            gVCF = True

        return snps, gVCF
    
    except Exception:
        console.print("Reference column is empty, please check your input file. Stopping scan.")
        sys.exit(1)


def call_trimming(content, chr):
    """ 
    The file must only contain SNPs. This function is necessary because the reference column has a different number in vcfs and in bim files 
    Args:   
        content: chunk of the file being read, containing the variants information. It must only contain SNPs, as indels are trimmed in the function.
        chr: chromosome of the chunk being processed, used to load the corresponding msgpack and get the matches.
    Returns:
        Calls get_matches() with the SNPs in the chunk, and adds the matches to the global variable final_results. It also returns the gVCF status of the chunk, which is True if any of the variants in the chunk has <NON_REF> in the ALT column, and False otherwise.
    """

    snps, gVCF = trimming_indels(content, 2)  # content : chr pos ref alt

    if len(snps) != 0:
        snps_dict = dict(zip(snps["position"].astype(int), snps["nucleotide"]))
        results = get_matches(snps_dict, chr)  # pass dict, not DataFrame
        final_results.append(results)
    else:
        console.print("There aren't FP SNPs in this chunk", style="bold red")
    
    return gVCF


def read_and_load(chunk):
    """
    Check the chromosome in the chunk and load only the msgpacks corresponding to that chromosome, 
    then call the function to trim indels and get the matches.
    Args:
        chunk: chunk of the file being read, containing the variants information. It must only contain SNPs, as indels are trimmed in the function.
    Returns:
        The gVCF status of the chunk, which is True if any of the variants in the chunk has <NON_REF> in the ALT column, and False otherwise.
    """
    for chromosome, group_content in chunk.groupby(chunk.columns[0]):
        chromosome_str = str(chromosome)
        if chromosome_str in chromosome_map:
            chr_key = chromosome_map[chromosome_str]  
            gVCF = call_trimming(group_content, chr_key)
            console.print("Variants being mapped from:", chr_key)
            
        else:
            console.print(f"Chromosome {chromosome_str} not found in chromosome map. Skipping variants from {chromosome_str}.")

    return gVCF


def read_chunks(complete_file, cols, n_matches=None, max_n_var=None):
    """
    Loads the file in batches to avoid loading it completely in memory.
    By default, all variants are read.

    If `n_matches` is provided (e.g. via -m in argparse), the file will be
    read in chunks until the total number of SNP matches reaches or exceeds
    `n_matches`. At that point the loop stops and the inference is done
    with the matches collected so far.

    Args:
        complete_file (str): path to the file to read
        cols (list[int]): column indices to load
        n_matches (int | None): total number of matches required before
                                  stopping early. If None, read all chunks.
    Returns:
        Chunks of 100.000 variants for the inference, until the stopping condition is met (if any).
    """
    chunk_counter = 0
    results = {}
    gVCF = False  # track if any chunk looks like gVCF
    console.print("[bold]\n++ INFORMATION INFERRED BY THE REF COLUMN ++ [/bold]\n ")
    try:

        for chunk in pd.read_csv(complete_file, sep="\t", comment="#", header=None, chunksize=100000, usecols=cols,):
            chunk_counter = chunk_counter+100000
            # Update global gVCF flag if any chunk reports True
            chunk_gvcf = read_and_load(chunk)
            gVCF = gVCF or chunk_gvcf

            # Update global results 
            results = gather_and_sum(final_results)

            # If user requested an early stop based on number of variants read
            if max_n_var is not None: 
                try:
                    if chunk_counter > max_n_var:
                        break
                except ValueError:
                        console.print("0 FP SPNs in this chunk", style="bold")

            # If user requested an early stop based on number of matches

            if n_matches is not None:
                results_matches = sum(results.values()) if results else 0
                if n_matches < results_matches:
                    # Enough evidence; stop reading more chunks
                    break
            

    finally:
        # Inference and messages
        try:
            if not results or max(results.values()) == 0:
                console.print("No SNPs found to infer the reference genome.", style="bold red")
            else:
                best_ref = max(results, key=results.get)
                best_matches = results[best_ref]
                total_matches = sum(results.values())

                if best_matches > total_matches / 2:
                    console.print(
                        f"[bold]Inferred Reference genome:[/bold] {best_ref}"
                    )
                else:
                    console.print(
                        "None of the versions has more than 50% of the total matches. "
                        "[bold]Reference genome version unknown.[/bold]",
                        style="red",
                    )

        except ValueError:
            # Covers cases like max([]) or similar
            console.print("No SNPs found to infer the reference genome.", style="bold red")

    if gVCF:
        console.print(f"[bold]gVCF by ALT column[/bold]")


def extract_columns(complete_file,  n_matches, max_n_var):
    """
    Extracts the columns of interest (chr, pos, ref, alt) and sends them to be processed. The inference is done with the matches collected until the stopping condition is met (if any).
    Args:        
        complete_file: opened file with the header
        n_matches: stop reading more chunks once the total number of matches reaches this value. By default, 5000 matches are required before stopping.
        max_n_var: stop reading more chunks once the total number of variants read exceeds this value. If None, read all chunks.    
    Returns:
        Calls the function to read the file in chunks and process them, with the columns of interest.
    """

    cols = [0, 1, 3, 4] # chr pos ref alt
    read_chunks(complete_file, cols,  n_matches, max_n_var)

def get_n_samples(header):
    """
    Counts the number of samples in the VCF file by looking at the last header line. 
    Args:        
        header: list of lines in the header of the VCF file
    Returns:
        The number of samples in the VCF file, inferred by the last header line. 
    """

    mandatory_columns = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] #fixed fields for variant information
    columns = header[-1].split("\t")

    try:
        for column in mandatory_columns:
            columns.remove(column)
        n_samples = len(columns)
        console.print(f"[bold]Number of samples:[/bold]", n_samples)
    except ValueError:
        console.print(f"[bold]Number of samples: [/bold]1")

def start_refgen_header(header, target_file):
    """
    Parses the information of the VCF header. If the contig information is present, it is used to infer the reference genome. 
    If not, a message is printed and the inference is done with the variants information. If the gVCF tag is present in the header, it is also reported.
    Args:        
        header: list of lines in the header of the VCF file
        target_file: path of the input file, used for the messages
    Returns:
        The reference genome inferred by the contig header information, if the file is gVCF according to the header
    """

    contig_list = [line for line in header if '##contig' in line and 'length' in line]
    if len(contig_list)!=0:
        contig_list2= [i.split(",") for i in contig_list]
        dict_contigs = {}

        for line in contig_list2:
            contig_id = None
            contig_length = None

            for part in line:
                if part.startswith('##contig=<ID='):
                    contig_id = part.replace('##contig=<ID=', '')
                elif part.startswith('length='):
                    contig_length = int(part.replace('length=', '').replace('>', ''))

            if contig_id is not None and contig_length is not None:
                dict_contigs[contig_id] = contig_length
        comparison(dict_contigs, target_file)  # run the next f

    else:
        console.print("[dark_orange]Contig information not in the header[/dark_orange] - [bold dark_orange]The reference genome can't be "
                      "inferred from the header information [/bold dark_orange]")
        
    
    gVCF= [line for line in header if '##ALT=<ID=NON_REF' in line]
    if gVCF:
        console.print(f"[bold]gVCF according to header[/bold]")
        

def extract_header(complete_file, input_file):
    """
    If present, extracts header and send it to match the refgenDetector database
    Args:        
        complete_file: opened file with the header
        Input_file: path of the input file, used for the messages
    Returns:
        The reference genome inferred by the contig header information 
    """
    header = []
    for line in complete_file:
        if line.startswith('#'):
            header.append(line.strip())
        else:
            break

    start_refgen_header(header, input_file)
    get_n_samples(header)

    

def open_vcf(input_file, n_matches, max_n_var):
    """
    Parse arguments and open the input VCF, compressed or not.
    Args:
         input_file: path of the input file
         n_matches (int | None): if provided, the function will stop reading more chunks once the total number of matches reaches. By default, 5000 matches are required before stopping. 
         max_n_var (int | None): if provided, the function will stop reading more chunks once the total number of variants read exceeds this value. If None, read all chunks.

    Returns:
        Calls the function to extract the header and the function to extract the columns of interest and infer the reference genome. The inference is done with the matches collected until the stopping condition is met (if any).
    """

    formats = ("vcf")
    compressed_formats = ("vcf.gz") ##TODO whem bgz read_chunks and possible read_and_load take a long time

    if input_file.endswith(compressed_formats):
        with gzip.open(input_file, "rt") as complete_file:
            extract_header(complete_file, input_file)
            extract_columns(complete_file, n_matches, max_n_var)

    elif input_file.endswith(formats):
        with open(input_file, "rt") as complete_file:
            extract_header(complete_file, input_file)
            extract_columns(complete_file, n_matches, max_n_var)
    else: 
        console.print(f"[bold][red] Only Formats Allowed: vcf and vcf.gz [/red][/bold]")




