# Independent Validation Results

These results were produced by independent AI agents with **no author involvement** in disease selection, seed gene curation, dataset choice, or configuration. Each agent was given only the GitHub URL and a skeptical PI-style prompt.

## Psoriasis (Gemini CLI Agent)

- **Agent:** Google Gemini CLI ("Jim"), zero prior context
- **Seed genes:** 96 (agent-curated from literature)
- **Datasets:** 5 skin tissue datasets, all GPL570 (3 discovery + 2 replication)
- **Core genes:** 50
- **Replication:** 50/50 survived (100%)
- **Meta-significant:** 28 (random effects)
- **Biology:** Textbook psoriasis — S100A7/A9/A12, DEFB4B, IL36G, CCL20, KRT16/17, BCL2
- **Blind recovery:** 47/50 (94.0%) — 3 missing (KYNU, PI3, TGM1)

Results: `psoriasis/independent_validation/`

## Colorectal Cancer (Claude Code CLI Agent)

- **Agent:** Claude Code CLI, zero prior context
- **Seed genes:** 515 (agent-curated from COSMIC, Open Targets, KEGG, literature)
- **Datasets:** 6 colon tissue datasets, all GPL570 (3 discovery + 3 replication)
- **Core genes:** 263
- **Replication:** 244/263 survived (92.8%)
- **Meta-significant:** 218 (random effects)
- **Biology:** APC, TP53, AXIN1, AURKA, TOP2A, ERBB2, MYC, CDH1, MMP7, VEGFA
- **Blind recovery:** 258/264 (97.7%) — 6 missing (BAK1, CCL20, CXCL2, DKK3, EPHB2, UGT1A1)

Results: `crc/independent_validation/`

## Significance

These are the first Riker Engine validations performed entirely by third parties. The agents independently selected diseases, curated seed genes from their own knowledge, chose GEO datasets, built YAML configs, and evaluated biological plausibility of results — all without coaching or access to expected outcomes.

Full validation report: see the project's independent validation documentation.
