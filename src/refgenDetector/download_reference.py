"""
Post-install script for refgenDetector.

Handles two installation scenarios:
  1. Cloned from GitHub: moves + decompresses .xz files from github_msgpacks/
  2. Installed via pip: downloads github_msgpacks/ from GitHub, then moves + decompresses.

Runs only once — skips everything if msgpacks/ already contains .msgpack files.
"""

import lzma
import shutil
import urllib.request
import urllib.error
import json
from pathlib import Path


# CONFIGURATION
GITHUB_USER    = "EGA-archive"       
GITHUB_REPO    = "refgenDetector"
GITHUB_BRANCH  = "main"

GITHUB_API_URL = (
    f"https://api.github.com/repos/{GITHUB_USER}/{GITHUB_REPO}"
    f"/contents/src/refgenDetector/github_msgpacks"
    f"?ref={GITHUB_BRANCH}"
)

GITHUB_RAW_BASE = (
    f"https://raw.githubusercontent.com/{GITHUB_USER}/{GITHUB_REPO}"
    f"/{GITHUB_BRANCH}/src/refgenDetector/github_msgpacks"
)

# PATH RESOLUTION

def _package_root() -> Path:
    """Return the installed package directory (contains __init__.py)."""
    # Works for both editable installs and regular pip installs
    try:
        import refgenDetector
        return Path(refgenDetector.__file__).parent
    except ImportError:
        # Fallback: relative to this script
        return Path(__file__).parent


def get_paths():
    pkg = _package_root()
    src_clone  = pkg / "github_msgpacks"   # present when repo is cloned
    src_pip    = pkg / "_downloaded_msgpacks"  # temp download dir for pip
    dst        = pkg / "msgpacks"
    return src_clone, src_pip, dst


# HELPERS

def is_already_setup(dst: Path) -> bool:
    """Return True if at least one .msgpack file exists in dst."""
    if not dst.exists():
        return False
    return any(dst.glob("*.msgpack"))


def decompress_xz(src_file: Path, dst_dir: Path) -> None:
    """Decompress a single .xz file into dst_dir, keeping the base name."""
    # e.g. foo.msgpack.xz  ->  dst_dir/foo.msgpack
    out_name = src_file.stem  # strips the last suffix (.xz)
    out_path = dst_dir / out_name
    print(f"  Decompressing {src_file.name} -> {out_path.name}")
    with lzma.open(src_file, "rb") as f_in, open(out_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)


def move_and_decompress(src_dir: Path, dst_dir: Path) -> None:
    """Move every .xz file from src_dir to dst_dir and decompress it."""
    dst_dir.mkdir(parents=True, exist_ok=True)
    xz_files = list(src_dir.glob("*.xz"))
    if not xz_files:
        print(f"  WARNING: no .xz files found in {src_dir}")
        return
    for xz_file in xz_files:
        decompress_xz(xz_file, dst_dir)
    print(f"  Done — {len(xz_files)} file(s) decompressed into {dst_dir}")


# DOWNLOAD LOGIC

def _github_api_request(url: str) -> list:
    """Fetch JSON from the GitHub contents API."""
    req = urllib.request.Request(url, headers={"User-Agent": "refgenDetector-installer"})
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            return json.loads(resp.read().decode())
    except urllib.error.HTTPError as e:
        raise RuntimeError(
            f"GitHub API error {e.code} for {url}.\n"
            "Check GITHUB_USER / GITHUB_REPO / GITHUB_BRANCH in post_install.py."
        ) from e


def download_github_msgpacks(dst_dir: Path) -> None:
    """Download every .xz file from the GitHub folder into dst_dir."""
    print(f"  Fetching file list from GitHub …")
    entries = _github_api_request(GITHUB_API_URL)

    xz_entries = [e for e in entries if e["name"].endswith(".xz")]
    if not xz_entries:
        raise RuntimeError("No .xz files found in the GitHub folder. Check the repo path.")

    dst_dir.mkdir(parents=True, exist_ok=True)
    for entry in xz_entries:
        raw_url  = f"{GITHUB_RAW_BASE}/{entry['name']}"
        out_path = dst_dir / entry["name"]
        print(f"  Downloading {entry['name']} ({entry.get('size', '?')} bytes) …")
        req = urllib.request.Request(raw_url, headers={"User-Agent": "refgenDetector-installer"})
        with urllib.request.urlopen(req, timeout=120) as resp, open(out_path, "wb") as f:
            shutil.copyfileobj(resp, f)

    print(f"  Downloaded {len(xz_entries)} file(s).")


# MAIN

def run():
    src_clone, src_pip, dst = get_paths()

    # Check if the reference files have already been downloaded and decompressed
    if is_already_setup(dst):
        print("[refgenDetector] Reference files already present — skipping setup.")
        return

    print("[refgenDetector] Setting up reference files …")

    # Intallation done by clonning the repo - github_msgpacks/ is already present
    if src_clone.exists() and any(src_clone.glob("*.xz")):
        print(f"  Detected clone install — using local {src_clone.name}/")
        move_and_decompress(src_clone, dst)

    # Installation done via pip — need to download from GitHub and decompress 
    else:
        print("  Detected pip install — downloading from GitHub …")
        try:
            download_github_msgpacks(src_pip)
            move_and_decompress(src_pip, dst)
        finally:
            # Clean up temp download dir regardless of success/failure
            if src_pip.exists():
                shutil.rmtree(src_pip, ignore_errors=True)

    print("[refgenDetector] Reference files ready.\n")


if __name__ == "__main__":
    run()