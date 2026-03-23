#!/usr/bin/env python3
"""Compare WGCNA vs Riker Engine results on ASD brain cortex data."""

import csv
import os
from collections import defaultdict
from datetime import datetime

BENCHMARKS_DIR = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(BENCHMARKS_DIR, "results")
RIKER_CORE = "/home/kai001/asd_full_run/output/bulk/phase4_core_genes.csv"
RIKER_VERDICTS = "/home/kai001/asd_full_run/output/bulk/phase5_verdicts.csv"
DATASETS = ["GSE28521", "GSE28475", "GSE64018"]


def load_csv(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def load_riker_core_genes():
    rows = load_csv(RIKER_CORE)
    genes = {}
    for r in rows:
        genes[r["gene"]] = {
            "cluster": int(r["cluster_id"]),
            "direction": r["direction"],
            "mean_log2fc": float(r["mean_log2fc"]),
        }
    return genes


def load_wgcna_modules(dataset_id):
    path = os.path.join(RESULTS_DIR, f"{dataset_id}_gene_modules.csv")
    rows = load_csv(path)
    gene_to_module = {}
    module_to_genes = defaultdict(set)
    for r in rows:
        gene = r["gene"]
        mod = int(r["module_number"])
        gene_to_module[gene] = mod
        module_to_genes[mod].add(gene)
    return gene_to_module, dict(module_to_genes)


def load_wgcna_trait(dataset_id):
    path = os.path.join(RESULTS_DIR, f"{dataset_id}_module_trait.csv")
    rows = load_csv(path)
    results = {}
    for r in rows:
        mod_name = r["module"]
        mod_num = int(mod_name.replace("ME", ""))
        results[mod_num] = {
            "correlation": float(r["correlation"]),
            "pvalue": float(r["pvalue"]),
            "significant": r["significant"] == "TRUE",
            "n_genes": int(r["n_genes"]) if r["n_genes"] else 0,
        }
    return results


def main():
    core_genes = load_riker_core_genes()
    core_gene_set = set(core_genes.keys())
    n_core = len(core_gene_set)

    # Load Riker cluster assignments
    riker_clusters = defaultdict(set)
    for gene, info in core_genes.items():
        riker_clusters[info["cluster"]].add(gene)

    lines = []
    lines.append("# WGCNA vs Riker Engine — ASD Brain Cortex Benchmark")
    lines.append(f"\nGenerated: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    lines.append(f"\n## Riker Engine Results (reference)")
    lines.append(f"\n- **Core genes**: {n_core}")
    lines.append(f"- **Clusters**: {len(riker_clusters)}")
    for cid in sorted(riker_clusters):
        genes = sorted(riker_clusters[cid])
        lines.append(f"  - Cluster {cid} ({len(genes)} genes): {', '.join(genes)}")

    # --- Per-dataset WGCNA analysis ---
    all_wgcna = {}
    all_trait = {}
    core_in_any_sig = set()

    for ds in DATASETS:
        gene_to_mod, mod_to_genes = load_wgcna_modules(ds)
        trait = load_wgcna_trait(ds)
        all_wgcna[ds] = (gene_to_mod, mod_to_genes)
        all_trait[ds] = trait

        # Which core genes are found in this dataset?
        found = core_gene_set & set(gene_to_mod.keys())
        missing = core_gene_set - found
        grey = {g for g in found if gene_to_mod[g] == 0}
        assigned = found - grey

        # Module distribution of core genes
        core_module_dist = defaultdict(set)
        for g in assigned:
            core_module_dist[gene_to_mod[g]].add(g)

        # Significant modules
        sig_modules = {m: t for m, t in trait.items() if t["significant"]}
        sig_module_nums = set(sig_modules.keys())

        # Core genes in significant modules
        core_in_sig = set()
        for g in assigned:
            if gene_to_mod[g] in sig_module_nums:
                core_in_sig.add(g)
                core_in_any_sig.add(g)

        lines.append(f"\n## {ds} — WGCNA Results")
        lines.append(f"\n- **Total modules found**: {len(mod_to_genes) - (1 if 0 in mod_to_genes else 0)}")
        lines.append(f"- **Significant modules (p<0.05)**: {len(sig_modules)}")
        if sig_modules:
            for m in sorted(sig_modules):
                t = sig_modules[m]
                lines.append(f"  - Module {m}: r={t['correlation']:.3f}, p={t['pvalue']:.4f}, {t['n_genes']} genes")
        lines.append(f"- **Core genes found in dataset**: {len(found)}/{n_core}")
        lines.append(f"- **Core genes in grey (unassigned)**: {len(grey)}")
        lines.append(f"- **Core genes assigned to modules**: {len(assigned)}")
        lines.append(f"- **Core genes in significant modules**: {len(core_in_sig)}/{n_core}")

        if missing:
            lines.append(f"- **Core genes absent from dataset**: {', '.join(sorted(missing))}")

        lines.append(f"\n### Core gene module distribution in {ds}")
        lines.append(f"\n| WGCNA Module | # Core Genes | Core Genes | Significant? |")
        lines.append(f"|---|---|---|---|")
        for mod_num in sorted(core_module_dist):
            genes = sorted(core_module_dist[mod_num])
            is_sig = "Yes" if mod_num in sig_module_nums else "No"
            total_in_mod = len(mod_to_genes.get(mod_num, []))
            lines.append(f"| {mod_num} ({total_in_mod} total) | {len(genes)} | {', '.join(genes)} | {is_sig} |")
        if grey:
            lines.append(f"| grey (unassigned) | {len(grey)} | {', '.join(sorted(grey))} | N/A |")

    # --- Cross-dataset consistency ---
    lines.append(f"\n## Cross-Dataset Consistency")

    # For each core gene, what module is it in across datasets?
    lines.append(f"\n### Core gene placement across datasets")
    lines.append(f"\n| Gene | Riker Cluster | {' | '.join(DATASETS)} | Consistent? |")
    lines.append(f"|---|---|{'---|' * len(DATASETS)}---|")

    # Check if WGCNA places each gene consistently
    # Since module numbers aren't comparable across datasets, we check whether
    # core genes that Riker groups together also land in the same WGCNA module
    for gene in sorted(core_gene_set):
        riker_cl = core_genes[gene]["cluster"]
        placements = []
        for ds in DATASETS:
            g2m = all_wgcna[ds][0]
            if gene in g2m:
                mod = g2m[gene]
                is_sig = mod in {m for m, t in all_trait[ds].items() if t["significant"]}
                suffix = "*" if is_sig else ""
                placements.append(f"{mod}{suffix}" if mod != 0 else "grey")
            else:
                placements.append("—")
        # "Consistent" = assigned to a non-grey module in all datasets where present
        non_missing = [p for p in placements if p != "—"]
        consistent = all(p != "grey" for p in non_missing) and len(non_missing) == len(DATASETS)
        lines.append(f"| {gene} | {riker_cl} | {' | '.join(placements)} | {'Yes' if consistent else 'No'} |")

    lines.append(f"\n*Module numbers marked with * are significant (p<0.05)*")

    # --- Riker cluster cohesion in WGCNA ---
    lines.append(f"\n### Riker cluster cohesion in WGCNA")
    lines.append(f"\nDo genes that Riker groups into the same cluster also land in the same WGCNA module?")

    for ds in DATASETS:
        g2m = all_wgcna[ds][0]
        lines.append(f"\n**{ds}**:")
        for cid in sorted(riker_clusters):
            cluster_genes = riker_clusters[cid]
            mods_seen = defaultdict(list)
            for g in sorted(cluster_genes):
                if g in g2m:
                    mods_seen[g2m[g]].append(g)
                else:
                    mods_seen["absent"].append(g)
            if len(mods_seen) == 1 and "absent" not in mods_seen:
                mod_id = list(mods_seen.keys())[0]
                status = "COHESIVE" if mod_id != 0 else "all grey"
            else:
                status = "SCATTERED"
            mod_summary = "; ".join(
                f"mod {m}: {', '.join(gs)}" if m != "absent" else f"absent: {', '.join(gs)}"
                for m, gs in sorted(mods_seen.items(), key=lambda x: str(x[0]))
            )
            lines.append(f"  - Riker cluster {cid} ({len(cluster_genes)} genes): **{status}** — {mod_summary}")

    # --- Summary statistics ---
    lines.append(f"\n## Summary")

    # Core genes in ANY significant module across ALL datasets
    core_never_sig = core_gene_set - core_in_any_sig
    lines.append(f"\n### Recovery of Riker's 35 core genes by WGCNA")
    lines.append(f"\n- Core genes in a significant WGCNA module in **at least one** dataset: "
                 f"**{len(core_in_any_sig)}/{n_core}** ({100*len(core_in_any_sig)/n_core:.1f}%)")
    lines.append(f"- Core genes **never** in any significant WGCNA module: "
                 f"**{len(core_never_sig)}/{n_core}** ({100*len(core_never_sig)/n_core:.1f}%)")
    if core_never_sig:
        lines.append(f"  - Genes WGCNA misses: {', '.join(sorted(core_never_sig))}")

    # Cross-dataset consistency: genes in non-grey module in all 3 datasets
    consistent_count = 0
    for gene in core_gene_set:
        all_assigned = True
        for ds in DATASETS:
            g2m = all_wgcna[ds][0]
            if gene not in g2m or g2m[gene] == 0:
                all_assigned = False
                break
        if all_assigned:
            consistent_count += 1
    lines.append(f"- Core genes assigned to a non-grey module in **all 3 datasets**: "
                 f"**{consistent_count}/{n_core}** ({100*consistent_count/n_core:.1f}%)")

    # WGCNA's significant module gene counts vs Riker
    lines.append(f"\n### Scale comparison")
    for ds in DATASETS:
        trait = all_trait[ds]
        sig_genes = set()
        mod_to_genes = all_wgcna[ds][1]
        for m, t in trait.items():
            if t["significant"]:
                sig_genes |= mod_to_genes.get(m, set())
        lines.append(f"- {ds}: WGCNA significant module genes = **{len(sig_genes)}** "
                     f"(vs Riker's {n_core} core genes)")

    lines.append(f"\n### Key takeaways")
    lines.append(f"\n1. **Specificity**: Riker Engine identifies {n_core} core genes through "
                 f"progressive filtering. WGCNA's significant modules contain orders of magnitude "
                 f"more genes, making biological interpretation harder.")
    lines.append(f"2. **Cross-dataset consistency**: Riker requires genes to survive across "
                 f"multiple datasets by design. WGCNA runs independently per dataset with no "
                 f"built-in cross-dataset validation.")
    lines.append(f"3. **Cluster cohesion**: See per-dataset analysis above for whether Riker's "
                 f"gene groupings are preserved or scattered by WGCNA.")

    lines.append(f"\n---\n*Benchmark run on Raspberry Pi 5 (Ghost)*")

    report = "\n".join(lines)
    report_path = os.path.join(BENCHMARKS_DIR, "benchmark_report.md")
    with open(report_path, "w") as f:
        f.write(report)
    print(f"Report written to: {report_path}")
    print(f"\n{'='*60}")
    print(report)


if __name__ == "__main__":
    main()
