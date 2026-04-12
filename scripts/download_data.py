#!/usr/bin/env python3
"""
Riker Engine — Data Download Script

Downloads GEO series matrix files, platform annotations, and HGNC gene data
needed to reproduce the five-disease validation.

Usage:
    python scripts/download_data.py asd          # Download ASD data only
    python scripts/download_data.py all           # Download all diseases
    python scripts/download_data.py --list        # Show what would be downloaded
    python scripts/download_data.py --platforms    # Download platform annotations only
    python scripts/download_data.py --hgnc        # Download HGNC data only

Downloaded files are placed in data/geo/<disease>/, data/platforms/, and data/hgnc/.
Seed gene files are already included in data/seeds/ (committed to the repo).

Note: Some RNA-seq datasets require reconstruction from supplementary data.
These are documented but not automatically downloaded. See the notes printed
after download completes.
"""

import argparse
import gzip
import hashlib
import os
import shutil
import sys
import urllib.request
import urllib.error

# Base directory: repo root
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR)
DATA_DIR = os.path.join(REPO_ROOT, "data")

# NCBI FTP base URLs
GEO_SERIES_FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"
GEO_PLATFORM_FTP = "https://ftp.ncbi.nlm.nih.gov/geo/platforms"

# HGNC download URL
# Primary: Google Cloud mirror (more reliable). Fallback: EBI FTP.
HGNC_URLS = [
    "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt",
    "https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt",
]


# ---------------------------------------------------------------------------
# Dataset definitions per disease
# ---------------------------------------------------------------------------

DATASETS = {
    "asd": {
        "description": "Autism Spectrum Disorder — 7 brain cortex datasets",
        "seed_file": "data/seeds/asd_sfari_genes.csv",
        "symbol_column": "symbol",
        "tissue": "brain",
        "geo_datasets": [
            # Discovery
            {"id": "GSE28521", "platform": "GPL6883", "role": "discovery",
             "tissue": "brain", "note": "Voineagu et al. — brain cortex microarray"},
            {"id": "GSE28475", "platform": "GPL6883", "role": "discovery",
             "tissue": "brain", "subseries": "GSE28475-GPL6883",
             "note": "Brain cortex microarray, GPL6883 subset"},
            {"id": "GSE64018", "platform": "ensembl", "role": "discovery",
             "tissue": "brain", "reconstructed": True,
             "note": "Gupta et al. — brain cortex RNA-seq. REQUIRES RECONSTRUCTION. "
                     "See docs/DATA_RECONSTRUCTION.md"},
            # Replication
            {"id": "GSE102741", "platform": "ensembl", "role": "replication",
             "tissue": "brain", "reconstructed": True,
             "note": "Brain cortex RNA-seq. REQUIRES RECONSTRUCTION. "
                     "See docs/DATA_RECONSTRUCTION.md"},
            {"id": "GSE18123", "platform": "GPL6244", "role": "replication",
             "tissue": "blood", "subseries": "GSE18123-GPL6244",
             "note": "LCL, GPL6244 subset"},
            {"id": "GSE26415", "platform": "GPL6480", "role": "replication",
             "tissue": "blood", "note": "Peripheral blood"},
            {"id": "GSE42133", "platform": "GPL10558", "role": "replication",
             "tissue": "blood", "note": "Blood microarray"},
        ],
    },
    "t2d": {
        "description": "Type 2 Diabetes — 4 pancreatic islet datasets",
        "seed_file": "data/seeds/t2d_open_targets_genes.csv",
        "symbol_column": "symbol",
        "tissue": "islet",
        "geo_datasets": [
            # Discovery
            {"id": "GSE41762", "platform": "GPL6244", "role": "discovery",
             "tissue": "islet", "note": "Rosengren et al. — pancreatic islets"},
            {"id": "GSE25724", "platform": "GPL96", "role": "discovery",
             "tissue": "islet", "note": "Taneera et al. — pancreatic islets"},
            {"id": "GSE20966", "platform": "GPL1352", "role": "discovery",
             "tissue": "islet", "note": "Marselli et al. — beta cells (LCM)"},
            # Replication
            {"id": "GSE86468", "platform": "ensembl", "role": "replication",
             "tissue": "islet", "reconstructed": True,
             "note": "Lawlor et al. — islet RNA-seq. REQUIRES RECONSTRUCTION. "
                     "See docs/DATA_RECONSTRUCTION.md"},
        ],
    },
    "ibd": {
        "description": "Inflammatory Bowel Disease — 6 intestinal mucosa datasets",
        "seed_file": "data/seeds/ibd_open_targets_genes.csv",
        "symbol_column": "symbol",
        "tissue": "colon",
        "geo_datasets": [
            # Discovery
            {"id": "GSE75214", "platform": "GPL6244", "role": "discovery",
             "tissue": "colon", "note": "Colonic + ileal mucosa, CD + UC"},
            {"id": "GSE16879", "platform": "GPL570", "role": "discovery",
             "tissue": "colon", "note": "Colonic mucosa, CD + UC, infliximab"},
            {"id": "GSE59071", "platform": "GPL6244", "role": "discovery",
             "tissue": "colon", "note": "Colonic mucosa, CD + UC"},
            # Replication
            {"id": "GSE87466", "platform": "GPL13158", "role": "replication",
             "tissue": "colon", "note": "Colonic mucosa, UC"},
            {"id": "GSE38713", "platform": "GPL570", "role": "replication",
             "tissue": "colon", "note": "Colonic mucosa, UC"},
            {"id": "GSE36807", "platform": "GPL570", "role": "replication",
             "tissue": "colon", "note": "Colonic tissue, CD + UC"},
        ],
    },
    "ad": {
        "description": "Alzheimer's Disease — 5 brain cortex datasets",
        "seed_file": "data/seeds/ad_curated_genes.csv",
        "symbol_column": "hgnc_symbol",
        "tissue": "brain",
        "geo_datasets": [
            # Discovery
            {"id": "GSE33000", "platform": "GPL4372", "role": "discovery",
             "tissue": "brain", "filtered": True,
             "note": "Prefrontal cortex, Rosetta/Merck. REQUIRES FILTERING: "
                     "remove Huntington's disease samples. See docs/DATA_RECONSTRUCTION.md"},
            {"id": "GSE44770", "platform": "GPL4372", "role": "discovery",
             "tissue": "brain", "note": "Prefrontal cortex, Rosetta/Merck"},
            {"id": "GSE118553", "platform": "GPL10558", "role": "discovery",
             "tissue": "brain", "filtered": True,
             "note": "REQUIRES FILTERING: extract temporal cortex samples only, "
                     "exclude AsymAD. See docs/DATA_RECONSTRUCTION.md"},
            # Replication
            {"id": "GSE5281", "platform": "GPL570", "role": "replication",
             "tissue": "brain", "filtered": True,
             "note": "REQUIRES FILTERING: extract superior frontal gyrus samples only. "
                     "See docs/DATA_RECONSTRUCTION.md"},
            {"id": "GSE15222", "platform": "GPL2700", "role": "replication",
             "tissue": "brain", "note": "Cortical tissue"},
        ],
    },
    "breast_cancer": {
        "description": "Breast Cancer — 5 breast tumor datasets",
        "seed_file": "data/seeds/breast_cancer_curated_genes.csv",
        "symbol_column": "symbol",
        "tissue": "breast",
        "geo_datasets": [
            # Discovery
            {"id": "GSE10810", "platform": "GPL570", "role": "discovery",
             "tissue": "breast", "note": "Breast tumor vs. healthy"},
            {"id": "GSE42568", "platform": "GPL570", "role": "discovery",
             "tissue": "breast", "note": "Breast cancer vs. normal"},
            {"id": "GSE15852", "platform": "GPL96", "role": "discovery",
             "tissue": "breast", "note": "Breast tumor vs. normal tissue"},
            # Replication
            {"id": "GSE45827", "platform": "GPL570", "role": "replication",
             "tissue": "breast", "note": "Breast cancer diagnosis"},
            {"id": "GSE65194", "platform": "GPL570", "role": "replication",
             "tissue": "breast", "note": "TNBC + HER2 + luminal vs. healthy"},
        ],
    },
    "ipf": {
        "description": "Idiopathic Pulmonary Fibrosis — 5 lung tissue datasets + 1 held-out",
        "seed_file": "data/seeds/ipf_curated_genes.csv",
        "symbol_column": "symbol",
        "tissue": "lung",
        "geo_datasets": [
            # Discovery
            {"id": "GSE32537", "platform": "GPL6244", "role": "discovery",
             "tissue": "lung", "note": "Boon et al. — IIP cohort, IPF/UIP subset"},
            {"id": "GSE53845", "platform": "GPL6480", "role": "discovery",
             "tissue": "lung", "note": "DePianto et al. — 40 IPF vs 8 controls"},
            {"id": "GSE24206", "platform": "GPL570", "role": "discovery",
             "tissue": "lung", "note": "Meltzer & Noble — early + advanced IPF"},
            # Replication
            {"id": "GSE110147", "platform": "GPL6244", "role": "replication",
             "tissue": "lung", "note": "Cecchini et al. — 22 IPF vs 11 controls"},
            {"id": "GSE10667", "platform": "GPL4133", "role": "replication",
             "tissue": "lung", "note": "Konishi et al. — UIP + acute exacerbation"},
            # Cold replication held-out dataset (not used by the pipeline —
            # used by scripts/cold_replication_ipf.py for independent validation)
            {"id": "GSE47460", "platform": "GPL6480", "role": "held_out",
             "tissue": "lung", "subseries": "GSE47460-GPL6480",
             "note": "LGRC cohort — 122 IPF/UIP + controls, held out for cold replication"},
        ],
    },
}

# All unique platform IDs across all diseases
ALL_PLATFORMS = sorted({
    ds["platform"]
    for disease_data in DATASETS.values()
    for ds in disease_data["geo_datasets"]
    if ds["platform"] != "ensembl"
})


def geo_series_url(gse_id: str, subseries: str | None = None) -> str:
    """Build the NCBI FTP URL for a GEO series matrix file."""
    # GSE28521 -> GSE28nnn
    numeric = gse_id.replace("GSE", "")
    bucket = f"GSE{numeric[:-3]}nnn"
    filename = f"{subseries or gse_id}_series_matrix.txt.gz"
    return f"{GEO_SERIES_FTP}/{bucket}/{gse_id}/matrix/{filename}"


def geo_platform_url(gpl_id: str) -> str:
    """Build the NCBI FTP URL for a GEO platform annotation file."""
    numeric = gpl_id.replace("GPL", "")
    bucket = f"GPL{numeric[:-3]}nnn" if len(numeric) > 3 else "GPLnnn"
    return f"{GEO_PLATFORM_FTP}/{bucket}/{gpl_id}/annot/{gpl_id}.annot.gz"


def download_file(url: str, dest: str, description: str = "") -> bool:
    """Download a file with progress reporting. Returns True on success."""
    if os.path.exists(dest):
        print(f"  [skip] {os.path.basename(dest)} already exists")
        return True

    os.makedirs(os.path.dirname(dest), exist_ok=True)
    desc = description or os.path.basename(dest)
    print(f"  [downloading] {desc}")
    print(f"    URL: {url}")

    try:
        req = urllib.request.Request(url, headers={"User-Agent": "RikerEngine/0.3.0"})
        with urllib.request.urlopen(req, timeout=120) as response:
            total = response.headers.get("Content-Length")
            total = int(total) if total else None

            tmp_dest = dest + ".tmp"
            downloaded = 0
            with open(tmp_dest, "wb") as f:
                while True:
                    chunk = response.read(8192)
                    if not chunk:
                        break
                    f.write(chunk)
                    downloaded += len(chunk)
                    if total:
                        pct = downloaded * 100 // total
                        mb = downloaded / (1024 * 1024)
                        print(f"\r    {mb:.1f} MB ({pct}%)", end="", flush=True)

            print()

            # If it's a .gz file for an .annot, decompress
            if dest.endswith(".annot") and tmp_dest.endswith(".annot.tmp"):
                # Already handled below
                pass

            os.rename(tmp_dest, dest)
            mb = os.path.getsize(dest) / (1024 * 1024)
            print(f"    Saved: {dest} ({mb:.1f} MB)")
            return True

    except urllib.error.HTTPError as e:
        print(f"    ERROR: HTTP {e.code} — {e.reason}")
        print(f"    URL: {url}")
        return False
    except urllib.error.URLError as e:
        print(f"    ERROR: {e.reason}")
        return False
    except Exception as e:
        print(f"    ERROR: {e}")
        return False


def download_and_decompress(url: str, dest: str, description: str = "") -> bool:
    """Download a .gz file and decompress it."""
    if os.path.exists(dest):
        print(f"  [skip] {os.path.basename(dest)} already exists")
        return True

    gz_dest = dest + ".gz"
    if not download_file(url, gz_dest, description):
        return False

    print(f"  [decompress] {os.path.basename(gz_dest)}")
    try:
        with gzip.open(gz_dest, "rb") as f_in:
            with open(dest, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove(gz_dest)
        return True
    except Exception as e:
        print(f"    ERROR decompressing: {e}")
        return False


def download_platforms(platforms: list[str] | None = None) -> dict[str, bool]:
    """Download platform annotation files."""
    targets = platforms or ALL_PLATFORMS
    platform_dir = os.path.join(DATA_DIR, "platforms")
    os.makedirs(platform_dir, exist_ok=True)

    results = {}
    print(f"\n{'='*60}")
    print(f"Downloading {len(targets)} platform annotation files")
    print(f"{'='*60}")

    for gpl in targets:
        dest = os.path.join(platform_dir, f"{gpl}.annot")
        url = geo_platform_url(gpl)
        results[gpl] = download_and_decompress(url, dest, f"{gpl} platform annotation")

    return results


def download_hgnc() -> bool:
    """Download the HGNC complete gene set."""
    hgnc_dir = os.path.join(DATA_DIR, "hgnc")
    os.makedirs(hgnc_dir, exist_ok=True)
    dest = os.path.join(hgnc_dir, "hgnc_complete_set.txt")

    print(f"\n{'='*60}")
    print("Downloading HGNC complete gene set")
    print(f"{'='*60}")

    for url in HGNC_URLS:
        if download_file(url, dest, "HGNC complete set (~17 MB)"):
            return True
    return False


def download_disease(disease: str) -> dict:
    """Download all GEO data for a disease. Returns summary."""
    if disease not in DATASETS:
        print(f"Unknown disease: {disease}")
        print(f"Available: {', '.join(DATASETS.keys())}")
        sys.exit(1)

    config = DATASETS[disease]
    geo_dir = os.path.join(DATA_DIR, "geo", disease)
    os.makedirs(geo_dir, exist_ok=True)

    print(f"\n{'='*60}")
    print(f"{config['description']}")
    print(f"Seed genes: {config['seed_file']}")
    print(f"{'='*60}")

    downloaded = []
    skipped = []
    failed = []
    manual = []

    # Collect unique platforms needed
    platforms_needed = set()

    for ds in config["geo_datasets"]:
        gse_id = ds["id"]
        subseries = ds.get("subseries")

        # Track platform
        if ds["platform"] != "ensembl":
            platforms_needed.add(ds["platform"])

        # Check if this needs reconstruction or filtering
        if ds.get("reconstructed") or ds.get("filtered"):
            manual.append(ds)
            print(f"\n  [manual] {gse_id} — {ds['note']}")
            continue

        # Download series matrix
        url = geo_series_url(gse_id, subseries)
        filename = f"{subseries or gse_id}_series_matrix.txt.gz"
        dest = os.path.join(geo_dir, filename)

        print(f"\n  --- {gse_id} ({ds['role']}, {ds['tissue']}) ---")
        if ds.get("note"):
            print(f"  {ds['note']}")

        if download_file(url, dest, f"{gse_id} series matrix"):
            downloaded.append(gse_id)
        else:
            failed.append(gse_id)

    # Download platforms for this disease
    if platforms_needed:
        platform_results = download_platforms(sorted(platforms_needed))
        for gpl, success in platform_results.items():
            if not success:
                failed.append(f"platform:{gpl}")

    return {
        "disease": disease,
        "downloaded": downloaded,
        "skipped": skipped,
        "failed": failed,
        "manual": manual,
    }


def print_summary(results: list[dict]):
    """Print a summary of all downloads."""
    print(f"\n{'='*60}")
    print("DOWNLOAD SUMMARY")
    print(f"{'='*60}")

    all_manual = []
    any_failed = False

    for r in results:
        disease = r["disease"]
        n_ok = len(r["downloaded"])
        n_fail = len(r["failed"])
        n_manual = len(r["manual"])
        total = n_ok + n_fail + n_manual

        status = "OK" if n_fail == 0 else "INCOMPLETE"
        print(f"\n  {disease.upper()}: {n_ok}/{total} downloaded, "
              f"{n_manual} need manual prep, {n_fail} failed  [{status}]")

        if r["failed"]:
            any_failed = True
            for f in r["failed"]:
                print(f"    FAILED: {f}")

        all_manual.extend(r["manual"])

    if all_manual:
        print(f"\n{'='*60}")
        print("DATASETS REQUIRING MANUAL PREPARATION")
        print(f"{'='*60}")
        print()
        print("The following datasets use RNA-seq data or require sample filtering.")
        print("They cannot be directly downloaded as standard GEO series matrices.")
        print("See docs/DATA_RECONSTRUCTION.md for detailed instructions.")
        print()
        for ds in all_manual:
            print(f"  {ds['id']}: {ds['note']}")

    print(f"\n{'='*60}")
    print("NEXT STEPS")
    print(f"{'='*60}")
    print()
    print("1. If any datasets require manual preparation (listed above),")
    print("   follow the instructions in docs/DATA_RECONSTRUCTION.md")
    print()
    print("2. Verify your data directory structure:")
    print("   data/")
    print("   ├── geo/<disease>/GSExxxxx_series_matrix.txt.gz")
    print("   ├── platforms/GPLxxxx.annot")
    print("   ├── hgnc/hgnc_complete_set.txt")
    print("   ├── seeds/  (already included in repo)")
    print("   └── mappings/ensembl_to_symbol.txt  (already included)")
    print()
    print("3. Update your config YAML to point to the downloaded files,")
    print("   or use the example configs in configs/examples/")
    print()
    print("4. Run the pipeline:")
    print("   riker run configs/examples/asd_bulk.yaml")


def list_downloads():
    """Show what would be downloaded without downloading."""
    print("Riker Engine — Data Dependencies\n")

    total_datasets = 0
    total_platforms = set()

    for disease, config in DATASETS.items():
        print(f"\n{disease.upper()}: {config['description']}")
        print(f"  Seed genes: {config['seed_file']}")
        for ds in config["geo_datasets"]:
            auto = "" if not ds.get("reconstructed") and not ds.get("filtered") else " [MANUAL]"
            print(f"  {ds['role']:12s}  {ds['id']:12s}  {ds['platform']:8s}  "
                  f"{ds['tissue']}{auto}")
            total_datasets += 1
            if ds["platform"] != "ensembl":
                total_platforms.add(ds["platform"])

    print(f"\nTotals: {total_datasets} datasets, {len(total_platforms)} platforms, "
          f"1 HGNC file")
    print(f"\nPlatforms: {', '.join(sorted(total_platforms))}")


def main():
    parser = argparse.ArgumentParser(
        description="Download GEO data for Riker Engine validation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "diseases",
        nargs="*",
        choices=list(DATASETS.keys()) + ["all"],
        help="Disease(s) to download data for",
    )
    parser.add_argument(
        "--list", action="store_true",
        help="List what would be downloaded without downloading",
    )
    parser.add_argument(
        "--platforms", action="store_true",
        help="Download all platform annotations only",
    )
    parser.add_argument(
        "--hgnc", action="store_true",
        help="Download HGNC gene data only",
    )

    args = parser.parse_args()

    if args.list:
        list_downloads()
        return

    if args.platforms:
        download_platforms()
        return

    if args.hgnc:
        download_hgnc()
        return

    if not args.diseases:
        parser.print_help()
        print("\nExamples:")
        print("  python scripts/download_data.py asd")
        print("  python scripts/download_data.py asd t2d ibd")
        print("  python scripts/download_data.py all")
        print("  python scripts/download_data.py --list")
        return

    diseases = list(DATASETS.keys()) if "all" in args.diseases else args.diseases

    # Always download HGNC first
    download_hgnc()

    # Download each disease
    results = []
    for disease in diseases:
        results.append(download_disease(disease))

    print_summary(results)


if __name__ == "__main__":
    main()
