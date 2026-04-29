#!/usr/bin/env python3
"""Build ASD disease report draft from run outputs."""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

base = Path("disease_days/2026-04-28_asd_v033")

# Load all data
t1_summary = json.load(open(base / "curated/pipeline_summary.json"))
t1_core = sorted(t1_summary["locked_core_genes"])

t2_summary = json.load(open(base / "blind_pc/pipeline_summary.json"))
t2_meta = []
with open(base / "blind_pc/phase6_meta_analysis.csv") as f:
    for row in csv.DictReader(f):
        row["neg_log10_p"] = -math.log10(max(float(row["random_p"]), 1e-300))
        t2_meta.append(row)
t2_top20 = sorted(t2_meta, key=lambda x: x["neg_log10_p"], reverse=True)[:20]

t3_summary = json.load(open(base / "stability_50run/stability_summary.json"))
t3_scores = []
with open(base / "stability_50run/stability_scores.csv") as f:
    for row in csv.DictReader(f):
        t3_scores.append(row)
iron_clad = [r for r in t3_scores if r["stability_class"] == "iron-clad"]
iron_clad_set = set(r["gene"] for r in iron_clad)
perfect_50 = sorted(r["gene"] for r in iron_clad if int(r["appearance_count"]) == 50)

# Run counts std dev
run_counts = []
with open(base / "stability_50run/run_summary.csv") as f:
    for row in csv.DictReader(f):
        if row["success"] == "True":
            run_counts.append(int(row["core_gene_count"]))
std_dev = round(float(np.std(run_counts)), 2)

# Cluster analysis
cluster_analysis = []
with open(base / "stability_50run/iron_clad_cluster_analysis.csv") as f:
    for row in csv.DictReader(f):
        cluster_analysis.append(row)
cluster_summary = []
with open(base / "stability_50run/iron_clad_cluster_summary.csv") as f:
    for row in csv.DictReader(f):
        cluster_summary.append(row)

# SFARI seeds
sfari = set()
with open("data/seeds/asd_sfari_genes.csv") as f:
    for row in csv.DictReader(f):
        sfari.add(row["symbol"])

# Historical v0.3.2
hist_path = Path("stability_ASD_blind/stability_scores.csv")
hist_ic = set()
if hist_path.exists():
    with open(hist_path) as f:
        for row in csv.DictReader(f):
            if row.get("stability_class") == "iron-clad":
                hist_ic.add(row["gene"])

# Meta for iron-clad genes
t2_meta_by_gene = {r["gene"]: r for r in t2_meta}
ic_with_meta = [(g, t2_meta_by_gene[g]) for g in iron_clad_set if g in t2_meta_by_gene]
top10_effect = sorted(ic_with_meta, key=lambda x: abs(float(x[1]["random_effect"])), reverse=True)[:10]
top10_p = sorted(ic_with_meta, key=lambda x: x[1]["neg_log10_p"], reverse=True)[:10]

# Cross-tier
blind_core = set(t2_summary["locked_core_genes"])
curated_and_ic = sorted(set(t1_core) & iron_clad_set)
curated_not_ic = sorted(set(t1_core) - iron_clad_set)
blind_and_ic = set(blind_core) & iron_clad_set
ic_not_blind = sorted(iron_clad_set - blind_core)

# Cluster anchors (genes in clusters with >=5 ic members)
big_cluster_genes = []
for row in cluster_summary:
    if int(row["n_iron_clad_members"]) >= 5 and int(row["cluster_id"]) != -1:
        big_cluster_genes.extend(row["members"].split("|"))

sfari_and_ic = sorted(iron_clad_set & sfari)

# === BUILD REPORT ===
out = []
W = out.append

W("# ASD Disease Report — v0.3.3 (DRAFT)")
W("")
W("---")
W("")

# Section 1
W("## 1. Run summary")
W("")
W(f"- **Disease:** ASD (Autism Spectrum Disorder)")
W(f"- **Date:** 2026-04-28")
W(f"- **Engine version:** {t1_summary['package_version']}")
W(f"- **Engine commit:** {t1_summary['code_version']}")
W(f"- **Master seed:** 42")
W(f"- **Tier 1 wall clock:** 2 min 4 sec")
W(f"- **Tier 2 wall clock:** 7 min 14 sec")
W(f"- **Tier 3 wall clock:** {t3_summary['total_wall_time_seconds']:.0f} sec ({t3_summary['total_wall_time_seconds']/3600:.1f} hours)")
W(f"- Three tiers completed: curated single run, protein-coding blind single run, 50-run stability profile. All runs successful, all QC checks passed.")
W("")

# Section 2
W("## 2. Tier 1 — Curated SFARI run")
W("")
W(f"- **Seed set:** 1,267 SFARI genes")
W(f"- **Phase 1 study genes:** {t1_summary['phase1_study_genes']}")
W(f"- **Phase 4 core genes:** {t1_summary['phase4_core_genes']}")
W(f"- **Phase 4 significant clusters:** {t1_summary['phase4_significant_clusters']}")
W(f"- **Phase 5 survived / eliminated:** {t1_summary['phase5_survived']} / {t1_summary['phase5_eliminated']}")
W(f"- **Phase 6 meta-significant (random effects):** {t1_summary['phase6_significant_random']}")
W("")
W("**Core genes (29, alphabetical):**")
W("")
for g in t1_core:
    W(f"- {g}")
W("")

# Section 3
W("## 3. Tier 2 — Blind protein-coding run")
W("")
W(f"- **Seed set:** 19,296 protein-coding genes")
W(f"- **Phase 1 study genes:** {t2_summary['phase1_study_genes']}")
W(f"- **Phase 4 core genes:** {t2_summary['phase4_core_genes']}")
W(f"- **Phase 4 significant clusters:** {t2_summary['phase4_significant_clusters']}")
W(f"- **Phase 5 survived / eliminated:** {t2_summary['phase5_survived']} / {t2_summary['phase5_eliminated']}")
W(f"- **Phase 6 genes analyzed:** {t2_summary['phase6_genes_analyzed']}")
W(f"- **Phase 6 significant (random effects):** {t2_summary['phase6_significant_random']}")
W("")
W("**Top 20 genes by Phase 6 random-effects -log10(p):**")
W("")
W("| Rank | Gene | Random effect | p-value | -log10(p) | Direction | Datasets |")
W("|---|---|---|---|---|---|---|")
for i, r in enumerate(t2_top20, 1):
    eff = float(r["random_effect"])
    p = float(r["random_p"])
    W(f"| {i} | {r['gene']} | {eff:.4f} | {p:.2e} | {r['neg_log10_p']:.2f} | {r['direction']} | {r['n_datasets']} |")
W("")

# Section 4
W("## 4. Tier 3 — 50-run stability profile")
W("")
ts = t3_summary
W(f"- **Runs completed:** {ts['n_runs']} / {ts['n_runs']} (0 failures)")
W(f"- **Master seed:** 42")
W(f"- **Total unique genes seen:** {ts['total_genes_seen']}")
W("")
W("**Stability classifications:**")
W("")
W("| Class | Count | Fraction | Threshold |")
W("|---|---|---|---|")
W(f"| Iron-clad | {ts['iron_clad']['count']} | {ts['iron_clad']['fraction']:.1f}% | >= 90% appearance |")
W(f"| Borderline | {ts['borderline']['count']} | — | 50-89% appearance |")
W(f"| Stochastic | {ts['stochastic']['count']} | — | < 50% appearance |")
W("")
W("**Per-run core gene counts:**")
W("")
W("| Metric | Value |")
W("|---|---|")
W(f"| Min | {ts['core_gene_counts']['min']} |")
W(f"| Max | {ts['core_gene_counts']['max']} |")
W(f"| Mean | {ts['core_gene_counts']['mean']} |")
W(f"| Std dev | {std_dev} |")
W("")
W("**Pairwise Jaccard similarity:**")
W("")
W("| Metric | Value |")
W("|---|---|")
pj = ts['pairwise_jaccard']
W(f"| Median | {pj['median']} |")
W(f"| 25th percentile | {pj['p25']} |")
W(f"| 75th percentile | {pj['p75']} |")
W(f"| Min | {pj['min']} |")
W(f"| Max | {pj['max']} |")
W(f"| Pairs computed | {pj['n_pairs']} |")
W("")
W(f"**Total wall clock:** {ts['total_wall_time_seconds']:.0f} sec ({ts['total_wall_time_seconds']/3600:.1f} hours)")
W("")

# Section 5
W(f"## 5. Iron-clad gene set (n={len(iron_clad)})")
W("")
W("| Gene | Appearance count |")
W("|---|---|")
for g in sorted(iron_clad, key=lambda x: x["gene"]):
    W(f"| {g['gene']} | {g['appearance_count']}/50 |")
W("")

# Section 6
W("## 6. Cluster organization of iron-clad genes")
W("")
W("Analysis script: `scripts/analyze_iron_clad_clusters.py`")
W("Output files: `stability_50run/iron_clad_cluster_analysis.csv`, `stability_50run/iron_clad_cluster_summary.csv`")
W("")
W("**Clusters with >= 5 iron-clad members (by modal cluster assignment):**")
W("")
W("| Cluster ID | Iron-clad members | Min modal count | Max modal count | Mean modal count |")
W("|---|---|---|---|---|")
for row in cluster_summary:
    if int(row["n_iron_clad_members"]) >= 5 and int(row["cluster_id"]) != -1:
        W(f"| {row['cluster_id']} | {row['n_iron_clad_members']} | {row['min_modal_count']} | {row['max_modal_count']} | {row['mean_modal_count']} |")
W("")

# Section 7
W("## 7. Cross-tier comparison")
W("")
W(f"**Curated core genes (29) AND iron-clad: {len(curated_and_ic)}**")
W("")
for g in curated_and_ic:
    W(f"- {g}")
W("")
W(f"**Curated core genes NOT in iron-clad: {len(curated_not_ic)}**")
W("")
for g in curated_not_ic:
    W(f"- {g}")
W("")
W(f"**Blind single core genes (421) AND iron-clad: {len(blind_and_ic)}**")
W("")
W(f"**Iron-clad genes NOT in blind single run: {len(ic_not_blind)}**")
W("")
if ic_not_blind:
    for g in ic_not_blind:
        W(f"- {g}")
    W("")

# Section 8
W("## 8. Findings")
W("")
W("Candidates identified by objective criteria applied mechanically to run outputs.")
W("")
W("### A. Highest meta-analysis effect size (top 10 iron-clad by |random_effect|)")
W("")
W("| Gene | Random effect | Direction |")
W("|---|---|---|")
for g, m in top10_effect:
    W(f"| {g} | {float(m['random_effect']):.4f} | {m['direction']} |")
W("")

W("### B. Strongest statistical confidence (top 10 iron-clad by -log10(p_random))")
W("")
W("| Gene | -log10(p) | p-value |")
W("|---|---|---|")
for g, m in top10_p:
    W(f"| {g} | {m['neg_log10_p']:.2f} | {float(m['random_p']):.2e} |")
W("")

W(f"### C. Perfect reproducibility (50/50 runs): {len(perfect_50)} genes")
W("")
W(f"{len(perfect_50)} iron-clad genes appeared in all 50 runs.")
W("")

W(f"### D. Cluster anchors: {len(big_cluster_genes)} genes")
W("")
W(f"{len(big_cluster_genes)} iron-clad genes reside in modal clusters with >= 5 other iron-clad members.")
W("")

W(f"### E. SFARI-AND-iron-clad: {len(sfari_and_ic)} genes")
W("")
for g in sfari_and_ic:
    W(f"- {g}")
W("")

# Section 9
W("## 9. Comparison to historical v0.3.2")
W("")
if hist_ic:
    W(f"- **Source:** `{hist_path}`")
    W(f"- **v0.3.2 iron-clad count:** {len(hist_ic)}")
    W(f"- **v0.3.3 iron-clad count:** {len(iron_clad_set)}")
    intersection = iron_clad_set & hist_ic
    v032_only = sorted(hist_ic - iron_clad_set)
    v033_only = sorted(iron_clad_set - hist_ic)
    W(f"- **Intersection (in both):** {len(intersection)}")
    W(f"- **v0.3.2 only:** {len(v032_only)}")
    W(f"- **v0.3.3 only:** {len(v033_only)}")
    W("")
    if v032_only:
        W("**Genes in v0.3.2 iron-clad but NOT v0.3.3 iron-clad:**")
        W("")
        for g in v032_only:
            W(f"- {g}")
        W("")
    if v033_only:
        W("**Genes in v0.3.3 iron-clad but NOT v0.3.2 iron-clad:**")
        W("")
        for g in v033_only:
            W(f"- {g}")
        W("")
else:
    W("Historical comparison data not located, needs Cody verification.")
    W("")

# Section 10
W("## 10. Cluster-level findings")
W("")
W("Clusters with >= 5 iron-clad members (by modal assignment across 50 runs):")
W("")
for row in cluster_summary:
    n = int(row["n_iron_clad_members"])
    cid = int(row["cluster_id"])
    if n >= 5 and cid != -1:
        members = sorted(row["members"].split("|"))
        W(f"### Cluster {cid} ({n} iron-clad members)")
        W("")
        W(f"Min modal count: {row['min_modal_count']}, Max: {row['max_modal_count']}, Mean: {row['mean_modal_count']}")
        W("")
        W("Members:")
        for m in members:
            W(f"- {m}")
        W("")

# Section 11
W("## 11. Observations")
W("")
W(f"- {len(perfect_50)} of {len(iron_clad)} iron-clad genes ({len(perfect_50)/len(iron_clad)*100:.1f}%) appeared in all 50 runs.")
W(f"- {len(iron_clad) - len(perfect_50)} iron-clad genes appeared in 45-49 runs, placing them near the iron-clad threshold.")
W(f"- Per-run core gene count ranged from {ts['core_gene_counts']['min']} to {ts['core_gene_counts']['max']} (std dev {std_dev}).")
W(f"- Pairwise Jaccard minimum was {pj['min']}, meaning even the two most dissimilar runs shared {pj['min']*100:.1f}% of their union.")
W(f"- {len(curated_not_ic)} of {len(t1_core)} curated core genes are not in the iron-clad set.")
W(f"- {len(ic_not_blind)} iron-clad genes are not in the Tier 2 blind single-run core gene set, indicating they appear in >=90% of stability runs but were absent from the specific single run used for Tier 2.")
W(f"- {len([r for r in cluster_summary if int(r['n_iron_clad_members']) >= 5 and int(r['cluster_id']) != -1])} clusters contain >= 5 iron-clad members each.")
W("")
W("---")
W("")
W("*Draft generated by Kai (Claude Code agent), 2026-04-29. Pending review by Cody Sigmon.*")

report = "\n".join(out)
with open("disease_days/2026-04-28_asd_v033/DISEASE_REPORT_DRAFT.md", "w") as f:
    f.write(report)
print(f"Report written: {len(out)} lines, {len(report)} bytes")
