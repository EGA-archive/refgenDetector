# src/refgenDetector/download_data.py
import os
import requests
import zipfile
from pathlib import Path

# URL of your reference data
REFGEN_URL = "https://crgcnag-my.sharepoint.com/:u:/g/personal/mimarin_crg_es/IQDa5CICZDAoRZmbfhBG3ZPEAVAnh4B4trHhGEHtBAIQxzI?e=bQZwX2"

# Directory of this file (same as refgenDetector_main.py)
PACKAGE_DIR = Path(__file__).parent

def get_refgen_data():

    zip_path = PACKAGE_DIR / "msgpacks.zip"
    extract_dir = PACKAGE_DIR / "msgpacks"
    # Only download if zip does not exist
    if not zip_path.exists():
        print("Downloading msgpacks.zip (~3GB)...")
        r = requests.get(REFGEN_URL, stream=True)
        r.raise_for_status()
        with open(zip_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
        print("Download complete!")

    # Only extract if folder does not exist
    if not extract_dir.exists():
        print("Extracting msgpacks.zip...")
        with zipfile.ZipFile(zip_path, "r") as zip_ref:
            zip_ref.extractall(extract_dir)
        print("Extraction complete!")

    return extract_dir
