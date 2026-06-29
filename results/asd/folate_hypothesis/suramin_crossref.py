#!/usr/bin/env python3
"""
Suramin Trial × Riker Engine Cross-Reference
=============================================
Tests whether ASD iron-clad genes overlap with CDR metabolic pathways
identified by Naviaux's suramin trial metabolomics.
"""

import csv
import os
from collections import defaultdict
from math import comb, factorial

OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))
TOTAL_GENES = 20000  # approximate human protein-coding genes

def load_iron_clad():
    """Load iron-clad genes with stability data."""
    genes = {}
    with open('/home/kai001/riker-engine/stability_ASD_blind/stability_scores.csv') as f:
        for row in csv.DictReader(f):
            if row['stability_class'] == 'iron-clad':
                genes[row['gene']] = row
    return genes

def load_core_genes():
    """Load core genes for direction info."""
    genes = {}
    with open('/home/kai001/riker-engine/output/asd_blind/phase4_core_genes.csv') as f:
        for row in csv.DictReader(f):
            genes[row['gene']] = row
    return genes

def load_meta():
    """Load meta-analysis results."""
    genes = {}
    with open('/home/kai001/riker-engine/output/asd_blind/phase6_meta_analysis.csv') as f:
        for row in csv.DictReader(f):
            genes[row['gene']] = row
    return genes

def load_suramin_genes():
    """Load suramin-response gene set organized by CDR domain."""
    genes = {}
    domains = defaultdict(list)
    with open(os.path.join(OUTPUT_DIR, 'suramin_response_genes.csv')) as f:
        for row in csv.DictReader(f):
            genes[row['gene']] = row
            domains[row['CDR_domain']].append(row['gene'])
    return genes, domains

def hypergeometric_p(k, n, K, N):
    """
    P(X >= k) for hypergeometric distribution.
    k = observed overlap
    n = iron-clad set size
    K = suramin gene set size (or domain size)
    N = total genes
    """
    p = 0.0
    for i in range(k, min(n, K) + 1):
        try:
            numerator = comb(K, i) * comb(N - K, n - i)
            denominator = comb(N, n)
            p += numerator / denominator
        except (OverflowError, ValueError):
            break
    return p

def run_crossref():
    iron_clad = load_iron_clad()
    core = load_core_genes()
    meta = load_meta()
    suramin_genes, domains = load_suramin_genes()

    n_iron = len(iron_clad)
    n_suramin = len(suramin_genes)

    print(f"Iron-clad genes: {n_iron}")
    print(f"Suramin-response genes: {n_suramin}")
    print(f"Background genome: {TOTAL_GENES}")
    print()

    # Overall overlap
    overlap = set(iron_clad.keys()) & set(suramin_genes.keys())
    print(f"=== OVERALL OVERLAP: {len(overlap)} genes ===")

    # Expected by chance
    expected = n_iron * n_suramin / TOTAL_GENES
    print(f"Expected by chance: {expected:.1f}")
    print(f"Fold enrichment: {len(overlap)/expected:.2f}x" if expected > 0 else "N/A")

    # Hypergeometric test
    p_overall = hypergeometric_p(len(overlap), n_iron, n_suramin, TOTAL_GENES)
    print(f"Hypergeometric p-value: {p_overall:.4e}")
    print()

    # List overlapping genes with details
    results = []
    for gene in sorted(overlap):
        s = suramin_genes[gene]
        direction = core.get(gene, {}).get('direction', 'N/A')
        log2fc = core.get(gene, {}).get('mean_log2fc', 'N/A')
        predicted = s['direction_predicted_in_ASD']
        cluster = core.get(gene, {}).get('cluster_id', 'N/A')

        # Check direction match
        if predicted == 'variable' or direction == 'N/A':
            match = 'undetermined'
        elif predicted == direction:
            match = 'MATCH'
        else:
            match = 'MISMATCH'

        results.append({
            'gene': gene,
            'CDR_domain': s['CDR_domain'],
            'metabolite': s['metabolite_or_pathway'],
            'role': s['role'],
            'predicted_direction': predicted,
            'observed_direction': direction,
            'log2FC': log2fc,
            'cluster': cluster,
            'match': match
        })

        if log2fc != 'N/A':
            log2fc = f"{float(log2fc):+.3f}"
        print(f"  {gene:<12} {s['CDR_domain']:<22} {direction:<6} (pred: {predicted:<10}) {match:<12} log2FC={log2fc}")

    print()

    # Domain-by-domain scoring
    print("=== DOMAIN SCORING ===")
    domain_results = {}
    for domain, genes_in_domain in sorted(domains.items()):
        n_domain = len(genes_in_domain)
        domain_overlap = set(genes_in_domain) & set(iron_clad.keys())
        n_hits = len(domain_overlap)
        expected_domain = n_iron * n_domain / TOTAL_GENES
        fold = n_hits / expected_domain if expected_domain > 0 else 0

        p_domain = hypergeometric_p(n_hits, n_iron, n_domain, TOTAL_GENES)

        # Direction consistency for this domain
        matches = 0
        mismatches = 0
        undetermined = 0
        for g in domain_overlap:
            r = [x for x in results if x['gene'] == g][0]
            if r['match'] == 'MATCH':
                matches += 1
            elif r['match'] == 'MISMATCH':
                mismatches += 1
            else:
                undetermined += 1

        determined = matches + mismatches
        dir_pct = matches / determined * 100 if determined > 0 else 0

        enriched = "YES" if p_domain < 0.05 and n_hits >= 2 else "no"

        domain_results[domain] = {
            'n_genes': n_domain,
            'n_hits': n_hits,
            'expected': expected_domain,
            'fold': fold,
            'p': p_domain,
            'enriched': enriched,
            'dir_match': matches,
            'dir_mismatch': mismatches,
            'dir_undetermined': undetermined,
            'dir_pct': dir_pct
        }

        hit_genes = ', '.join(sorted(domain_overlap)) if domain_overlap else '(none)'
        print(f"  {domain:<22} {n_hits:>2}/{n_domain:<3} hits  expected={expected_domain:.1f}  "
              f"fold={fold:.1f}x  p={p_domain:.3e}  enriched={enriched}  "
              f"dir={matches}/{determined} ({dir_pct:.0f}%)  genes: {hit_genes}")

    print()

    # Summary statistics
    n_enriched = sum(1 for d in domain_results.values() if d['enriched'] == 'YES')
    total_matches = sum(1 for r in results if r['match'] == 'MATCH')
    total_mismatches = sum(1 for r in results if r['match'] == 'MISMATCH')
    total_determined = total_matches + total_mismatches
    overall_dir_pct = total_matches / total_determined * 100 if total_determined > 0 else 0

    print("=== SCORING ===")
    print(f"Domains enriched (p<0.05, ≥2 hits): {n_enriched} / {len(domains)}")
    print(f"Direction consistency: {total_matches}/{total_determined} = {overall_dir_pct:.0f}%")
    print()

    # Apply rubric
    print("=== DECISION RUBRIC ===")
    if n_enriched >= 3 and overall_dir_pct > 65:
        decision = "PROCEED"
        print(f"≥3 domains enriched ({n_enriched}) + >65% direction match ({overall_dir_pct:.0f}%) → PROCEED")
    elif n_enriched >= 2 and overall_dir_pct > 50:
        decision = "MAYBE"
        print(f"2 domains enriched ({n_enriched}) + >50% direction match ({overall_dir_pct:.0f}%) → MAYBE")
    elif n_enriched <= 1 or overall_dir_pct < 50:
        decision = "DROP"
        print(f"≤1 domain enriched ({n_enriched}) or <50% direction match ({overall_dir_pct:.0f}%) → DROP")
    else:
        decision = "DROP"
        print(f"Does not meet any threshold → DROP")

    print(f"\n*** DECISION: {decision} ***\n")

    # Save results
    with open(os.path.join(OUTPUT_DIR, 'suramin_crossref_results.csv'), 'w') as f:
        writer = csv.DictWriter(f, fieldnames=['gene', 'CDR_domain', 'metabolite', 'role',
                                                'predicted_direction', 'observed_direction',
                                                'log2FC', 'cluster', 'match'])
        writer.writeheader()
        writer.writerows(results)

    # Save domain summary
    with open(os.path.join(OUTPUT_DIR, 'suramin_domain_summary.csv'), 'w') as f:
        writer = csv.DictWriter(f, fieldnames=['domain', 'n_genes', 'n_hits', 'expected',
                                                'fold', 'p', 'enriched', 'dir_match',
                                                'dir_mismatch', 'dir_undetermined', 'dir_pct'])
        writer.writeheader()
        for domain, data in sorted(domain_results.items()):
            writer.writerow({'domain': domain, **data})

    return decision, n_enriched, overall_dir_pct, results, domain_results


if __name__ == '__main__':
    decision, n_enriched, dir_pct, results, domains = run_crossref()
