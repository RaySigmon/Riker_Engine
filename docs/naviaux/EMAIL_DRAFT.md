# Naviaux Re-engagement Email — DRAFT

**STATUS: DRAFT — DO NOT SEND until experimental proposal is complete**
**Last updated:** 2026-05-10

---

Subject: ASD mitochondrial cluster — updated validation + collaboration proposal

---

Dr. Naviaux,

I hope this finds you well. I wanted to follow up on our March correspondence about the mitochondrial/OxPhos cluster in ASD brain transcriptomics, with an update on the validation work and a specific collaboration proposal.

**Engine update:** Since we last spoke, the Riker Engine underwent a version update (v0.3.2 → v0.3.3.2) that corrected a documented fold-change calculation issue. I re-validated the ASD analysis under the corrected engine: 40 of 41 genes in the energy metabolism cluster survive (HSPD1 was the single loss); CYC1 improved from borderline to iron-clad stability; all 8 FDR-significant genes remain iron-clad; the overall iron-clad count increased from 376 to 394. The direction pattern (23 up, 18 down) is unchanged. The CDR-consistent interpretation we discussed is fully intact.

**Extended validation:** I've now completed the full 8-disease validation across six tissue types, with 50-run stability profiling on each. The engine produces consistent, reproducible results across neurodevelopmental (ASD), neurodegenerative (AD), fibrotic (IPF), cancer (breast, colorectal), autoimmune (psoriasis, IBD), and metabolic (T2D) diseases. An empirical model characterizing when the engine's discovery mode produces candidate gene sets vs transcriptomic-scale measurements has been validated across 7 diseases.

**Collaboration proposal:** In our earlier exchange, I noted that a critical next question is whether the mitochondrial cluster genes normalize when the CDR is pharmacologically resolved. I've prepared a staged experimental proposal built around this question. The first stage requires only data sharing — if your lab has existing transcriptomic data from suramin-treated mouse brain (MIA model or similar), I can run the Riker Engine analysis myself and test whether the 15 mito-core genes show direction-reversal with treatment. This costs your lab no experimental effort and produces a publishable result within weeks. Later stages, contingent on the initial results, would involve targeted experimental validation.

I've attached the full proposal document with pre-specified directional predictions, power analysis, and publication strategy for each outcome scenario. The analysis is designed so that every result direction — positive, negative, or mixed — produces a publishable finding for Molecular Autism.

I would welcome the opportunity to discuss this further at your convenience.

Respectfully,
Ray Sigmon
Alpha Research Labs
alphalabsincorp@gmail.com
GitHub: github.com/RaySigmon/Riker_Engine

---

## Notes for revision before sending

- [ ] Confirm whether to mention the 7 additional diseases by name or keep it brief
- [ ] Decide whether to attach the full proposal or offer to send on request
- [ ] Verify Naviaux's current email is still rnaviaux@health.ucsd.edu
- [ ] Check if there's a better time to send (avoid holidays, conference season)
- [ ] Have proposal document complete and reviewed before sending
- [ ] Sit on final draft 48 hours before sending
