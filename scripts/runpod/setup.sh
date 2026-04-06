#!/bin/bash
# RunPod Setup Script — Riker Engine Benchmarks & Blind Runs
# Run this first after connecting to the pod.
#
# Prerequisites: RunPod with 32GB+ RAM, Python 3.11+
# Estimated setup time: 20-40 minutes (mostly R package compilation)
# Estimated task runtime: 4-6 hours for all tasks

set -e

echo "============================================"
echo "Riker Engine — RunPod Setup"
echo "============================================"

# 0. System dependencies for R packages
# WGCNA and MEGENA require compiled C/Fortran code
echo "[0/6] Installing system dependencies for R..."
apt-get update -qq
apt-get install -y -qq \
    r-base r-base-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libgfortran5 \
    gfortran \
    build-essential \
    2>/dev/null || echo "Some packages may already be installed"

# 1. Clone and install
echo "[1/6] Cloning repo..."
cd /root
if [ -d "Riker_Engine" ]; then
    echo "  Repo already exists, pulling latest..."
    cd Riker_Engine && git pull
else
    git clone https://github.com/RaySigmon/Riker_Engine.git
    cd Riker_Engine
fi
pip install -e ".[clustering,dev]"

# 2. Download data for all diseases
echo "[2/6] Downloading GEO data..."
python scripts/download_data.py all

# 3. Handle reconstructed datasets
echo "[3/6] Checking for reconstructed datasets..."
echo "  NOTE: AD datasets requiring filtering (GSE33000, GSE118553, GSE5281)"
echo "  and RNA-seq reconstructions (GSE64018, GSE102741, GSE86468)"
echo "  must be prepared manually. See docs/DATA_RECONSTRUCTION.md"
echo "  If these were prepared locally and committed, they should be in the repo."

# 4. Install R packages for benchmarks
echo "[4/6] Installing R packages (this takes 15-30 minutes)..."
Rscript -e '
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install BiocManager first
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install WGCNA (depends on many Bioconductor packages)
cat("Installing WGCNA...\n")
BiocManager::install("WGCNA", ask=FALSE, update=FALSE, quiet=TRUE)

# Install MEGENA
cat("Installing MEGENA...\n")
BiocManager::install("MEGENA", ask=FALSE, update=FALSE, quiet=TRUE)

# Install MetaDE (CRAN)
cat("Installing MetaDE...\n")
install.packages("MetaDE", quiet=TRUE)

# Verify all loaded
cat("\nVerifying R packages:\n")
for (pkg in c("WGCNA", "MEGENA", "MetaDE")) {
    if (requireNamespace(pkg, quietly = TRUE)) {
        cat(sprintf("  %s: OK\n", pkg))
    } else {
        cat(sprintf("  %s: FAILED\n", pkg))
    }
}
'

# 5. Verify Python installation
echo "[5/6] Verifying Python installation..."
python -m pytest tests/ -q --tb=line

# 6. Check available RAM
echo "[6/6] System check..."
echo "  RAM: $(free -h | awk '/^Mem:/ {print $2}') total, $(free -h | awk '/^Mem:/ {print $7}') available"
echo "  CPUs: $(nproc)"
echo "  Disk: $(df -h /root | awk 'NR==2 {print $4}') free"

echo ""
echo "============================================"
echo "Setup complete!"
echo ""
echo "Run all tasks:"
echo "  bash scripts/runpod/run_all.sh"
echo ""
echo "Or run individual tasks:"
echo "  riker run configs/runpod/ad_blind.yaml"
echo "  riker run configs/runpod/breast_cancer_blind.yaml"
echo "  Rscript benchmarks/wgcna_asd_benchmark.R"
echo "  Rscript benchmarks/megena_asd_benchmark.R"
echo "  Rscript benchmarks/metade_asd_benchmark.R"
echo "============================================"
