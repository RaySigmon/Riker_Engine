#!/bin/bash
# RunPod Execution Script — Run all benchmark & blind tasks
# Assumes setup.sh has already been run successfully.
#
# Tasks:
# 1. AD blind run (~30-60 min, needs 16+ GB)
# 2. Breast cancer blind run (~20-40 min, needs 16+ GB)
# 3. WGCNA benchmark re-run with timing (~3-4 hours)
# 4. MEGENA benchmark (~1-2 hours)
# 5. MetaDE benchmark (~30 min)

set -e
cd /root/Riker_Engine  # adjust if cloned elsewhere

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_DIR="runpod_logs/${TIMESTAMP}"
mkdir -p "$LOG_DIR"

echo "============================================"
echo "Riker Engine — RunPod Batch Execution"
echo "Log directory: $LOG_DIR"
echo "Started: $(date)"
echo "============================================"

# --- Task 1: AD Blind Run ---
echo ""
echo "[Task 1/5] AD Blind Run (all expressed genes)"
echo "  Expected: 14,442 study genes, needs ~16GB for consensus matrix"
time riker run configs/runpod/ad_blind.yaml 2>&1 | tee "$LOG_DIR/ad_blind.log"
echo "[Task 1] Complete: $(date)"

# --- Task 2: Breast Cancer Blind Run ---
echo ""
echo "[Task 2/5] Breast Cancer Blind Run (all expressed genes)"
time riker run configs/runpod/breast_cancer_blind.yaml 2>&1 | tee "$LOG_DIR/breast_cancer_blind.log"
echo "[Task 2] Complete: $(date)"

# --- Task 3: WGCNA Benchmark ---
echo ""
echo "[Task 3/5] WGCNA Benchmark (with timing)"
time Rscript benchmarks/wgcna_asd_benchmark.R 2>&1 | tee "$LOG_DIR/wgcna_benchmark.log"
echo "[Task 3] Complete: $(date)"

# --- Task 4: MEGENA Benchmark ---
echo ""
echo "[Task 4/5] MEGENA Benchmark"
time Rscript benchmarks/megena_asd_benchmark.R 2>&1 | tee "$LOG_DIR/megena_benchmark.log"
echo "[Task 4] Complete: $(date)"

# --- Task 5: MetaDE Benchmark ---
echo ""
echo "[Task 5/5] MetaDE Benchmark"
time Rscript benchmarks/metade_asd_benchmark.R 2>&1 | tee "$LOG_DIR/metade_benchmark.log"
echo "[Task 5] Complete: $(date)"

echo ""
echo "============================================"
echo "All tasks complete: $(date)"
echo "Logs: $LOG_DIR/"
echo ""
echo "Next steps:"
echo "  1. Review logs for errors"
echo "  2. Copy results to repo: results/ad/blind/, results/breast_cancer/blind/"
echo "  3. Copy benchmark outputs to benchmarks/results/"
echo "  4. Commit and push"
echo "============================================"
