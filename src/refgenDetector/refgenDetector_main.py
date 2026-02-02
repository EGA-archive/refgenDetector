#!/usr/bin/env python

""" refgenDetector.py: Script to infer the reference genome used to create a BAM or CRAM"""

__author__ = "Mireia Marin Ginestar"
__version__ = "2.0"
__maintainer__ = "Mireia Marin Ginestar"
__email__ = "mireia.marin@crg.eu"
__status__ = "Developement"

version = "2.0.0"

import os
import sys
import argparse
import gzip
import pysam
import psutil
import time
from rich.console import Console

# Add the parent directory to the Python path
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.insert(0, parent_dir)

from reference_genome_dictionaries import *
from exceptions.NoFileException import *
from aligment_files import *
from variant_files import *

console = Console()


def monitor_resources(func):
    """Decorator to print resource usage (CPU, memory, I/O, runtime)."""
    def wrapper(*args, **kwargs):
        process = psutil.Process()
        start_time = time.time()
        start_cpu_time = process.cpu_times()

        result = func(*args, **kwargs)

        end_time = time.time()
        end_cpu_time = process.cpu_times()
        duration = end_time - start_time

        cpu_user = end_cpu_time.user - start_cpu_time.user
        cpu_system = end_cpu_time.system - start_cpu_time.system
        total_cpu_time = cpu_user + cpu_system
        memory_usage = process.memory_info().rss / (1024 * 1024)  # MB
        io_counters = process.io_counters()
        bytes_read = io_counters.read_bytes
        bytes_written = io_counters.write_bytes

        print(f"Execution time: {duration:.2f} seconds")
        print(f"CPU time used: {total_cpu_time:.2f} seconds")
        print(f"Memory usage (RSS): {memory_usage:.2f} MB")
        if bytes_read > (1024*1024):
            print(f"Disk I/O - Read: {bytes_read / (1024 * 1024):.2f} MB, Written: {bytes_written / (1024 * 1024):.2f} MB")
        elif bytes_read > 1024:
            print(f"Disk I/O - Read: {bytes_read / 1024:.2f} KB, Written: {bytes_written / 1024:.2f} KB")
        else:
            print(f"Disk I/O - Read: {bytes_read:.2f} Bytes, Written: {bytes_written:.2f} Bytes")

        return result

    return wrapper


def run_main(args):
    """Main logic of the tool (isolated from CLI parsing)."""
    console.print(f"[bold]* Running refgenDetector v.{version} *[/bold]")
    console.print(f"---")
    console.print(f"[bold]++ INFORMATION INFERRED BY THE HEADER ++[/bold]\n")
    console.print(f"[bold]File:[/bold] {args.file}")
    try:
        if args.type == "Header":  
            console.print("[bold]File type:[/bold] BAM/CRAM header")  
            process_data_txt(args.file, args.md5, args.assembly)
        elif args.type in ["VCF"]:
            console.print("[bold]File type:[/bold] VCF") 
            open_vcf(args.file, args.matches, args.max_n_var)
        else:
            console.print("[bold]File type:[/bold] BAM/CRAM") 
            process_data_bamcram(args.file, args.md5, args.assembly)
    except OSError:
        console.print(f"[red]The file {args.file} provided in --file can't be opened."
                      f"\nRun [bold]refgenDetector -h[/bold] to get more information about the usage of the tool.")
    console.print(f"---")

def main():
    parser = argparse.ArgumentParser(prog="INFERRING THE REFERENCE GENOME USED TO ALIGN BAM OR CRAM FILE")
    parser.add_argument("-f", "--file", help="Input file path", required=True)
    parser.add_argument("-t", "--type", choices=["BAM/CRAM", "Header", "VCF", "BIM"], required=True,
                        help="Type of files to analyze.")
    parser.add_argument("--md5", action="store_true", help="Print md5 values if present in header.")
    parser.add_argument("-a", "--assembly", action="store_true", help="Print assembly if present in header.")
    parser.add_argument("-v", "--max_n_var", type=int, help="Maximum number of variants to read before stopping inference. The file is processed in chunks of 100,000 variants, so this value must be a multiple of 100,000 (e.g. 100000, 200000, 300000, ...).") 
    parser.add_argument("-m", "--matches", type=int, default=5000, help="Number of matches required before stopping. [DEFAULT:5000]")
    parser.add_argument("-r", "--resources", action="store_true",
                        help="When set, print execution time, CPU, memory, and disk I/O usage.")
    args = parser.parse_args()

    # Conditional resource monitoring
    if args.resources:
        wrapped = monitor_resources(run_main)
        wrapped(args)
    else:
        run_main(args)


if __name__ == "__main__":
    main()
