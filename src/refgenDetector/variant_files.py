import gzip
import sys
import time
import pandas as pd
from aligment_files import comparison
from chromosomes_dict import *
from rich.console import Console
import argparse

final_results = []
console = Console(highlight=False)
_pickle_cache = {}  # Add this line here


def gather_and_sum(lists):
    """It gathers and sums all the matches calculated in get_matches()"""
    cumulative_sums = {}
    for lst in lists:
        for key, value in lst:
            if key in cumulative_sums:
                cumulative_sums[key] += value
            else:
                cumulative_sums[key] = value
    console.print("Matches:", cumulative_sums)
    return cumulative_sums


def get_matches(snps, chr):
    """ Loads the pkl file depending on the chr that are present on the chunk. Compares the reference column from our
     input file and from the pkls, the number of matches to each version are returned to read_chunks()"""
    global _pickle_cache
    start = time.time()
    genome_versions = ["hg18", "GRCh37", "GRCh38", "T2T"]
    matches = []
    
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
            matches_count = (merged_df['nucleotide'] == merged_df[version_name]).sum()
            matches.append([version_name, matches_count])

        except FileNotFoundError as e:
            console.print(f"Warning: {version_name}-{chr}.pkl not found")
            continue

    console.print("Getting matches. Took:", time.time() - start, "s")

    return matches


def trimming_indels(content, ref):
    """If a row is longer than one position it is deleted, deleting this way any indels"""
    try:
        # Masks for SNPs (handle missing values safely)
        m_ref = content.iloc[:, 2].str.len().eq(1)  # REF length == 1
        m_alt = content.iloc[:, 3].str.len().eq(1)  # ALT length == 1

        # Filter once, then take the two columns you need and make an explicit copy
        snps = content.loc[m_ref & m_alt, [content.columns[1], content.columns[ref]]].copy()

        # Rename and uppercase nucleotide column (vectorized)
        snps.columns = ["position", "nucleotide"]
        snps["nucleotide"] = snps["nucleotide"].astype("string").str.upper()

        return snps

    except Exception:
        console.print("Reference column is empty, please check your input file. Stopping scan.")
        sys.exit(1)

def call_trimming(content, file_type, chr):
    """ The file must only contain SNPs. This function is necessary because the reference column has a different
    number in vcfs and in bim files """

    start = time.time()
    if file_type == "VCF":
        snps = trimming_indels(content, 2)  # content : chr pos ref alt
    elif file_type == "BIM":
        snps = trimming_indels(content, 3)
    console.print("Trimming indels. Took:", time.time() - start, "s")

    if len(snps) != 0:
        results = get_matches(snps, chr)
        final_results.append(results)
    else:
        console.print("There aren't FP SNPs in this chunk", style="bold")


def read_and_load(chunk, input_file, file_type):
    """To avoid loading a big pkl with information from all the chromosomes we first check which chr are there in
    the current chunk and then load only the necessary pkls
    """
    for chromosome, group_content in chunk.groupby(chunk.columns[0]):
        chromosome_str = str(chromosome)
        if chromosome_str in chromosome_map:
            chr_key = chromosome_map[chromosome_str]
            console.print("Variants being mapped from:", chr_key)
            call_trimming(group_content, file_type, chr_key)
        else:
            console.print(f"Chromosome {chromosome_str} not found in chromosome map. Skipping variants from {chromosome_str}.")


def read_chunks(complete_file, input_file, file_type, cols,  chunks, matches):
    """Loads the file in batches to avoid loading completely on memory. If there are enough matches to define the
    version, the loop will stop. If it has loaded more chunks than desired, the loop breaks too. For bim files,
    the inferred version must have at least 50% of the matches"""

    results = {}

    try:
        for chunk in pd.read_csv(complete_file, sep="\t", comment='#', header=None, chunksize=100000, usecols=cols):
            read_and_load(chunk, input_file, file_type)
            results = gather_and_sum(final_results)

            try:
                if max(results.values()) > matches:
                    break

            except ValueError:
                    console.print("0 FP SPNs in this chunk", style="bold")
            
    finally:
        try:
            if max(results.values()) == 0:
                console.print("No SNPs found to infer the reference genome.", style="bold red")
            else:
                if max(results.values()) > sum(results.values())/2:
                    console.print("Inferred Reference genome:", max(results, key=results.get), style="bold dark_cyan")
                else:
                    console.print("Any of the versions have more than 50% of the total matches.  [bold] Reference "
                                  "genome version unknown. m", style="red")
        except ValueError as e:
            console.print("No SNPs found to infer the reference genome.", style="bold red")

def extract_columns(complete_file, input_file, file_type, chunks, matches):
    "Loads only the interesting columns"
    try:
        if file_type == "VCF":
            cols = [0, 1, 3, 4] # chr pos ref alt
            read_chunks(complete_file, input_file, file_type, cols, chunks, matches)
        elif file_type == "BIM":
            cols = [0, 3, 5, 4] # chr pos ref alt
            read_chunks(complete_file, input_file, file_type, cols, chunks, matches)
    except ValueError as e:
        console.print(f"ValueError: An error occurred while processing the columns for file {input_file}: {e}")
        sys.exit(1)
    except Exception as e:
        console.print(f"An unexpected error occurred while extracting columns from file {input_file}. Check your file "
              f"contains all mandatory columns. Error: {e}")
        sys.exit(1)


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
        
def extract_header(complete_file, input_file):
    "Extracts header and send it to match the refgenDetector database"
    header = []
    for line in complete_file:
        if line.startswith('#'):
            header.append(line.strip())
        else:
            break
    start_refgen_header(header, input_file)


def open_file_vcf_bim(input_file, file_type, chunks, matches):
    """
    Checks if the file_type is correct. Then, opens the file if it endswith vcf/bim or vcf.gz/bim.gz.
    """
    console.print(f"Starting pre-processing for [[light_slate_blue]{input_file}[/light_slate_blue]]")
    formats = ("vcf", "bim")
    compressed_formats = ("vcf.gz", "bim.gz")
    dict_formats = {
        "vcf.gz": "VCF",
        "bim.gz": "BIM",
        "vcf": "VCF",
        "bim": "BIM"
    }

    for ext, expected_type in dict_formats.items():
        if input_file.endswith(ext):
            if file_type != expected_type:
                raise ValueError(f"File type should be {expected_type} for .{ext} files.")
            break

    try:
        if input_file.endswith(compressed_formats):
            with gzip.open(input_file, "rt") as complete_file:
                extract_header(complete_file, input_file)
                extract_columns(complete_file, input_file, file_type, chunks, matches)
        elif input_file.endswith(formats):
            with open(input_file, "rt") as complete_file:
                extract_header(complete_file, input_file)
                extract_columns(complete_file, input_file, file_type, chunks, matches)
        else:
            console.print("File format not supported. Please state a vcf or bim file, gz compressed or uncompressed")
    except FileNotFoundError:
        console.print(f"Error: The file {input_file} was not found in the stated path.")
        sys.exit(1)
    except OSError as e:
        console.print(f"Error: OS error occurred while handling the file {input_file}: {e}")
        sys.exit(1)
    except Exception as e:
        console.print(f"An unexpected error occurred while processing the file {input_file}: {e}")
        sys.exit(1)