# Significance Index (Si) — PCOS Gene Prioritization Framework

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![Preprint](https://img.shields.io/badge/Preprint-Research%20Square-blue)](https://www.researchsquare.com/article/rs-8610143/v1)
[![Status: Under Review](https://img.shields.io/badge/Status-Under%20Review-yellow)](https://www.researchsquare.com/article/rs-8610143/v1)

This repository contains the analysis code associated with the preprint:

> **Rajavelu S & De D.** *Evidence-based genetic variants to gene mapping and prioritization uncovers distinct molecular pathophysiology and therapeutic landscape in polycystic ovary syndrome patients of different ethnicities.* Research Square (2025). https://doi.org/10.21203/rs.3.rs-8610143/v1

---

## Overview

Polycystic ovary syndrome (PCOS) is a complex endocrine disorder with strong genetic underpinnings, yet the majority of GWAS-identified variants reside in non-coding regions, making effector gene assignment and therapeutic translation challenging. This study implements an integrative, **population-stratified gene prioritization framework** applied to East Asian (EAS) and European (EUR) PCOS GWAS cohorts across seven disease-relevant tissues.

The framework maps non-coding GWAS variants to candidate effector genes through three complementary regulatory layers, integrates protein–protein interaction (PPI) network diffusion, and scores all genes using a unified **Significance Index (Si)** — a composite score ranging from 0 to 5. The approach uncovers population-specific molecular phenotypes: **metabolic dysregulation** predominates in EAS PCOS, while **inflammatory and immune-related signatures** are more prominent in EUR PCOS.

---

## Methods Summary

### 1. GWAS Data & SNP Processing

- PCOS GWAS summary statistics were retrieved from the [NHGRI-EBI GWAS Catalog](https://www.ebi.ac.uk/gwas) (EFO:0000660)
- 6 EAS studies (Han Chinese & Korean) and 5 EUR studies were included
- Lead SNPs at genome-wide significance (P < 5×10⁻⁵) were expanded via linkage disequilibrium (LD) using [LDProxy](https://ldlink.nci.nih.gov/) (r² ≥ 0.5, ±500 kb window; 1000 Genomes Phase 3)
- A composite SNP score integrating association P-value and LD r² was computed for each variant

### 2. Core Gene Identification

Three complementary approaches identify **core genes** — those with direct genomic evidence linking them to PCOS-associated SNPs:

| Category | Description |
|----------|-------------|
| **rGenes** | Genes in physical proximity to SNPs, weighted by chromatin state (ChromHMM / Roadmap Epigenomics) and constrained by topologically associating domain (TAD) boundaries |
| **qGenes** | Genes associated with SNPs via quantitative trait loci (eQTL, pQTL, etc.) from GTEx and QTLbase |
| **cGenes** | Genes connected to SNP-overlapping enhancers through chromatin conformation (Hi-C / ABC model) |

Core genes receive additional weight via functional annotation (fGenes from PCOSKB), phenotype evidence (pGenes from HPO), and disease evidence (dGenes from DisGeNet).

### 3. Peripheral Gene Identification

**Random Walk with Restart (RWR)** on the integrated PPI network (Cheng et al. interactome + STRING v12.0; ~17,500 genes) propagates influence from core genes to retrieve **peripheral genes** — functionally connected candidates lacking direct genomic evidence.

### 4. Significance Index (Si) Scoring

Scores from all six predictor categories (rGene, qGene, cGene, fGene, pGene, dGene) are combined using **Fisher's meta-analysis method** (via empirical CDF transformation) into a single composite P-value, then rescaled to a **0–5 Si score** for each gene across all tissues and populations.

### 5. Benchmarking & Validation

- **Benchmarking:** AUC-based comparison against Open Targets (text mining & genetic association) using PCOS Gold Standard Positives/Negatives (GSPs/GSNs) from ChEMBL v34
- **Validation 1 (TSEA):** Target Set Enrichment Analysis (GSEA) of known PCOS drug targets (Phase 2+) in the Si-ranked gene list
- **Validation 2 (CREEDS):** Gene Set Enrichment Analysis of PCOS patient DEGs from GEO (GSE5090, GSE6798, GSE8157, GSE98595, GSE155489, GSE138518, GSE106724)

### 6. Downstream Analyses

- **Pathway Enrichment:** Reactome pathways (MSigDB C2, v2023.2) via one-sided Fisher's exact test with Benjamini-Hochberg FDR correction
- **Pathway Crosstalk:** Minimum-scoring subgraph extraction (prize-collecting Steiner tree) on MSigDB-annotated PPI subnet using the `dnet` R package
- **Node Removal Analysis:** Single and combinatorial node removal to assess network robustness and drug repurposing potential
- **Modular Analysis:** Community detection using spin-glass model (`clusterspinglass`, igraph v2.0.3) across cross-tissue merged crosstalk networks
- **Cluster Analysis:** Self-organizing map (supraHex) on Si scores of crosstalk genes across tissues to define therapeutic space; druggability overlay using fpocket-annotated PDB structures
---

## Contact

**Dr. Debojyoti De** (Corresponding Author)  
Assistant Professor, Department of Biotechnology  
National Institute of Technology Durgapur, West Bengal – 713209, India  
📧 dde.bt@nitdgp.ac.in

---
## Citation

If you use this code, please cite:

```bibtex
@article{rajavelu2025pcos,
  title   = {Evidence-based genetic variants to gene mapping and prioritization uncovers 
             distinct molecular pathophysiology and therapeutic landscape in polycystic 
             ovary syndrome patients of different ethnicities},
  author  = {Rajavelu, Sindhuja and De, Debojyoti},
  journal = {Research Square (Preprint)},
  year    = {2025},
  doi     = {10.21203/rs.3.rs-8610143/v1},
  url     = {https://www.researchsquare.com/article/rs-8610143/v1}
}
---

## License

This work is licensed under a [Creative Commons Attribution 4.0 International License (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/).
