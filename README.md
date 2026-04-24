# Significance-Index: Evidence-based Genetic Variant to Gene Mapping & Prioritization

## About the Project
This repository contains the source code for the Significance-Index framework, as detailed in the study: **"Evidence-based genetic variants to gene mapping and prioritization uncovers distinct molecular pathophysiology and therapeutic landscape in polycystic ovary syndrome patients of different ethnicities."** 📄 **Preprint:** [10.21203/rs.3.rs-8610143/v1](https://www.researchsquare.com/article/rs-8610143/v1)

The Significance-Index is an integrative, population-aware framework designed to bridge the gap between non-coding GWAS variants and clinically actionable targets. By systematically combining tissue-specific regulatory functional genomics, long-range chromatin interactions, genome-wide quantitative trait loci (QTLs), and protein–protein interaction (PPI) networks, the pipeline accurately ranks disease-associated genes. 

At the core of the methodology is a percentile rank-based scoring method—utilizing an empirical cumulative distribution function (eCDF)—to calculate the **Significance Index (Si rating)** for gene prioritization.

## Key Features
* **Multi-Predictor Variant Mapping:** Integrates functional annotations and tissue-specific regulatory data (e.g., ovarian tissue) to map variants to target genes.
* **Significance Index (Si) Computation:** Uses an eCDF framework to calculate quantitative traits (qGenes) and functional evidence scores.
* **Population-Specific Analytics:** Built to uncover ethnic heterogeneity (e.g., distinguishing between East Asian and European molecular phenotypes).
* **Advanced Network Topologies:** Features community detection, node removal analysis, and inter-tissue cross-talk modularity.
* **Therapeutic Repurposing:** Maps genetics to current therapeutics (G2CT) to identify novel druggable targets.

## Prerequisites & Environment
The pipeline relies heavily on **R** and **Python** for statistical analysis and network generation. 
* **R:** Utilized for core statistical tests (one-sided Fisher's exact tests, Benjamini-Hochberg FDR) and network topology analysis (e.g., `igraph` for pathway crosstalk and community detection).
* **Python:** Utilized for data processing and generating visualizations.

*Note: The visualization scripts in this repository are pre-configured to export high-resolution, publication-quality graphics (600 DPI JPEGs), including binary heatmaps, ridge plots of Si ratings, and hierarchical circular subgraphs.*

## Repository Structure & Usage
The codebase is organized into modular steps matching the study's analytical workflow:

1.  **GWAS Processing (`/01_gwas_processing`):** Scripts for handling summary-level data and identifying core genes under genomic influence.
2.  **Gene Annotation & Scoring (`/02_si_scoring`):** Implementation of the eCDF formulas to calculate the final Si rating for qGenes and cGenes.
3.  **Validation & Benchmarking (`/03_benchmarking`):** Code to evaluate prioritization performance via AUC and compare against external datasets like Open Targets.
4.  **Network & Pathway Analysis (`/04_network_analysis`):** Scripts to build hierarchical Reactome subgraphs (Signal Transduction, Immune System, Disease) and evaluate pathway enrichment odds ratios.
5.  **Therapeutics (`/05_drug_repurposing`):** Algorithms for evaluating druggability and therapeutic landscape mapping.

*(Note: Detailed execution instructions and command-line arguments are provided within each subdirectory's local README).*

## Citation
If you utilize this framework or code in your research, please cite our preprint:
> Rajavelu, S., & De, D. (2026). *Evidence-based genetic variants to gene mapping and prioritization uncovers distinct molecular pathophysiology and therapeutic landscape in polycystic ovary syndrome patients of different ethnicities*. Research Square. [https://doi.org/10.21203/rs.3.rs-8610143/v1](https://doi.org/10.21203/rs.3.rs-8610143/v1)

## License
[Insert License Here - e.g., MIT License]
