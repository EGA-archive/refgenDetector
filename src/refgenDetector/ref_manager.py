import sys
import json
from pathlib import Path
try:
    # Works when installed as a pip package
    from .reference_genome_dictionaries import *
except ImportError:
    # Works when run directly as a script
    from reference_genome_dictionaries import *


DEFAULT_MAJOR_RELEASES = dict(major_releases)

CUSTOM_DB = Path("custom_references.json")

if CUSTOM_DB.exists():
    with open(CUSTOM_DB) as f:
        custom_refs = json.load(f)

    major_releases.update(custom_refs)


def load_from_fai(file_path):
    """
    Loads the contig information from a .fai file into a dictionary.
    Args:
        file_path (str): path to the .fai file to load

    Returns:
        dict: A dictionary where keys are contig names and values are their lengths.
    """

    contigs = {}

    with open(file_path) as f:
        for line in f:
            if line.strip() == "":
                continue

            parts = line.strip().split("\t")

            if len(parts) < 2:
                raise ValueError(f"Invalid .fai line: {line}")

            contigs[parts[0]] = int(parts[1])

    if not contigs:
        raise ValueError("Empty .fai file")

    return contigs


def load_custom_db():
    """
    Loads the custom reference database from a JSON file. If the file does not exist, returns an empty dictionary.

    Returns:
        dict: A dictionary containing the custom references, where keys are reference names and values are their data.
    """
    if CUSTOM_DB.exists():
        with open(CUSTOM_DB) as f:
            return json.load(f)
    return {}


def save_custom_db(db):
    """
    Saves the custom reference database to a JSON file.
    Args:
        db (dict): A dictionary containing the custom references, where keys are reference names and values are their data.

    """
    with open(CUSTOM_DB, "w") as f:
        json.dump(db, f, indent=4)


def get_all_references():
    """
    Merge default + custom references
    Returns:
        dict: A dictionary containing all references, where keys are reference names and values are their data.
    """
    custom = load_custom_db()
    merged = dict(major_releases)
    merged.update(custom)
    return merged


def compare_contigs(c1, c2):
    """
    Return True if contig dictionaries match exactly
    Args:
        c1 (dict): Contig dictionary of the new reference, where keys are contig names and values are their lengths.
        c2 (dict): Contig dictionary of an existing reference, where keys are contig names and values are their lengths.
    """
    return c1 == c2


def find_matching_reference(new_contigs, all_refs):
    """
    Check if contigs match an existing reference. Avoid adding duplicates to the database. 
    Args:
        new_contigs (dict): Contig dictionary of the new reference, where keys are contig names and values are their lengths.
        all_refs (dict): A dictionary containing all references, where keys are reference names and values are their data, including the contig dictionary under the key "
    Returns:        
        str | None: The name of the matching reference if a match is found, or None if not. 
    """

    for name, data in all_refs.items():
        if compare_contigs(new_contigs, data["ref_gen"]):
            return name

    return None


def add_reference(ref_name, species, fai_file):
    """
    Add a new reference to the custom database, after checking that it doesn't match an existing reference. The contig information is loaded from the provided .fai file.
    Args:
        ref_name (str): Name of the new reference to add.
        species (str): Species of the new reference to add. 
        fai_file (str): Path to the .fai file containing the contig information of the new reference. The .fai file must have the format: contig_name \t contig_length \t ... (other columns are ignored).
    Returns:    
        If a matching reference is found, it prints a message with the name of the matching reference and aborts the addition. 
        If no match is found, it adds the new reference to the custom database and prints a success message.
    """

    fai_file = Path(fai_file)

    if not fai_file.exists():
        raise FileNotFoundError(f"File not found: {fai_file}")

    if fai_file.suffix != ".fai":
        raise ValueError("Input must be a .fai file")

    print(f"Loading contigs from: {fai_file}")
    contigs = load_from_fai(fai_file)

    all_refs = get_all_references()

    # Duplicate detection
    match = find_matching_reference(contigs, all_refs)
    if match:
        print(f"This reference matches existing reference: {match}")
        print("Aborting to avoid duplication.")
        return

    db = load_custom_db()

    if ref_name in db:
        print(f"Warning: '{ref_name}' already exists and will be overwritten")

    db[ref_name] = {
        "ref_gen": contigs,
        "build": ref_name,
        "species": species
    }

    save_custom_db(db)

    print(f"Reference '{ref_name}' added successfully.")


def list_references():
    """
    List all available references, including both default and custom ones. 
    It indicates the origin of each reference (default or custom) and prints the total number of references.
    Args:     None
    Returns:
    """
    default = DEFAULT_MAJOR_RELEASES
    custom = load_custom_db()
    all_refs = get_all_references()

    print("\nAvailable references:\n")

    for name, data in all_refs.items():
        origin = "custom" if name in custom else "default"
        print(f"- {name} ({data['species']}) [{origin}]")

    print(f"\nTotal: {len(all_refs)} references")


def remove_reference(ref_name):
    """
    Remove a reference from the custom database by its name. 
    It checks if the reference exists in the custom database before attempting to remove it, and prints a message indicating whether the removal was successful or if the reference was not found.
    Args:      
        ref_name (str): Name of the reference to remove. It must be a reference that was added to the custom database, as default references cannot be removed.
    Returns:     
        If the reference is found and removed successfully, it prints a success message.
    """
    db = load_custom_db()

    if ref_name not in db:
        print(f"'{ref_name}' not found in custom references.")
        return

    del db[ref_name]
    save_custom_db(db)

    print(f"Removed reference '{ref_name}'")


def main():
    """
    Command-line interface for managing the reference database. It supports three commands:
        - add <name> <species> <fai>: Adds a new reference to the custom database with the specified name, species, and contig information loaded from the provided .fai file.
        - list: Lists all available references, including both default and custom ones, indicating their origin and the total number of references.
        - remove <name>: Removes a reference from the custom database by its name. Only references that were added to the custom database can be removed, default references cannot be removed.
    """

    if len(sys.argv) < 2:
        print(
            "Usage:\n"
            "  add <name> <species> <fai>\n"
            "  list\n"
            "  remove <name>"
        )
        sys.exit(1)

    command = sys.argv[1]

    if command == "add":
        if len(sys.argv) != 5:
            print("Usage: add <name> <species> <genome.fai>")
            sys.exit(1)

        add_reference(sys.argv[2], sys.argv[3], sys.argv[4])

    elif command == "list":
        list_references()

    elif command == "remove":
        if len(sys.argv) != 3:
            print("Usage: remove <name>")
            sys.exit(1)

        remove_reference(sys.argv[2])

    else:
        print(f"Unknown command: {command}")


if __name__ == "__main__":
    main()