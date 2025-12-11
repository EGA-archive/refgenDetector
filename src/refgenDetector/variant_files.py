import gzip
import sys
import time
import pandas as pd
from dns.inet import inet_pton
from aligment_files import comparison
from chromosomes_dict import *
from rich.console import Console
import json
import os

final_results = []
console = Console(highlight=False)
_pickle_cache = {}  # load pkl only once
PKL_DIR = os.path.join(os.path.dirname(__file__), "pkls")

def gather_and_sum(lists):
    """It gathers and sums all the matches calculated in get_matches()"""
    cumulative_sums = {}
    for lst in lists:
        for key, value in lst:
            if key in cumulative_sums:
                cumulative_sums[key] += value
            else:
                cumulative_sums[key] = value
    console.print(f"[bold]Matches: [/bold]", cumulative_sums)
    return cumulative_sums


def get_matches(snps, chr):
    """ Loads the pkl file depending on the chr that are present on the chunk. Compares the reference column from our
     input file and from the pkls, the number of matches to each version are returned to read_chunks()"""
    global _pickle_cache
    start = time.time()
    genome_versions = ["hg18", "GRCh37", "GRCh38", "T2T"]
    matches = []
    
    if not os.path.isdir(PKL_DIR):
        console.print(
            f"[bold red]Error:[/bold red] Required directory not found: {PKL_DIR}\n Execution aborted. \n"
            "If you haven't downloaded yet the pkl folder: https://crgcnag-my.sharepoint.com/:u:/g/personal/mimarin_crg_es/IQByWKuqxkRMR7IfzwTy6CFiAf5Lleoa9V8QnKi3ByJPzHU?e=fbcfUh"
        )
        sys.exit(1)   # Hard stop (cleaner than raising inside loops)
    
    # Load and check sequentially
    for version_name in genome_versions:
        cache_key = f"{version_name}-{chr}"
        if cache_key not in _pickle_cache:
            try:
                _pickle_cache[cache_key] = pd.read_pickle(f"./pkls/{cache_key}.pkl")
            except FileNotFoundError:
                continue
        
        version_df = _pickle_cache[cache_key]
        try:
            version_df = pd.read_pickle(f"./pkls/{version_name}-{chr}.pkl")
            merged_df = version_df.join(snps.set_index('position'), how='inner')
            matches_count = int((merged_df['nucleotide'] == merged_df[version_name]).sum())
            matches.append([version_name, matches_count])

        except FileNotFoundError as e:
            console.print(f"Warning: {version_name}-{chr}.pkl not found")
            continue

    #console.print("Getting matches. Took:", time.time() - start, "s")

    return matches


def trimming_indels(content, ref):
    """If a row is longer than one position it is deleted, deleting this way any indels"""
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
    """ The file must only contain SNPs. This function is necessary because the reference column has a different
    number in vcfs and in bim files """

    snps, gVCF = trimming_indels(content, 2)  # content : chr pos ref alt

    if len(snps) != 0:
        results = get_matches(snps, chr)
        final_results.append(results)
    else:
        console.print("There aren't FP SNPs in this chunk", style="bold red")
    
    return gVCF


def read_and_load(chunk):
    """To avoid loading a big pkl with information from all the chromosomes we first check which chr are there in
    the current chunk and then load only the necessary pkls
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

    If `min_matches` is provided (e.g. via -m in argparse), the file will be
    read in chunks until the total number of SNP matches reaches or exceeds
    `min_matches`. At that point the loop stops and the inference is done
    with the matches collected so far.

    Args:
        complete_file (str): path to the file to read
        cols (list[int]): column indices to load
        min_matches (int | None): total number of matches required before
                                  stopping early. If None, read all chunks.
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

def extract_columns(complete_file, n_matches, max_n_var):
    "Loads only the interesting columns"
    cols = [0, 1, 3, 4] # chr pos ref alt
    read_chunks(complete_file, cols, n_matches, max_n_var)

def get_n_samples(header):

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
    "Extracts header and send it to match the refgenDetector database"
    header = []
    for line in complete_file:
        if line.startswith('#'):
            header.append(line.strip())
        else:
            break

    start_refgen_header(header, input_file)
    get_n_samples(header)

    

def open_vcf(input_file, n_matches, max_n_var):

    formats = ("vcf")
    compressed_formats = ("vcf.gz")

    if input_file.endswith(compressed_formats):
        with gzip.open(input_file, "rt") as complete_file:
            extract_header(complete_file, input_file)
            extract_columns(complete_file, n_matches, max_n_var)

    elif input_file.endswith(formats):
        with open(input_file, "rt") as complete_file:
            extract_header(complete_file, input_file)
            extract_columns(complete_file, n_matches, max_n_var)




