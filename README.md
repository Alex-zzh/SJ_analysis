# SJ_analysis

This repository contains code and resources designed to assist readers in reproducing specific results presented in our study and applying our methods to their own personalized datasets.

**Note:** Standard analyses (such as basic differential expression analysis, pathway enrichment, and specific filtering thresholds) are detailed in the Methods section of our manuscript and are not repeated here. We aim to provide unique scripts and custom workflows that facilitate the replication of our specific findings.

## Repository Contents

### 1. Ortholog Analysis
*   **Script:** `c1_ortholog_analysis.R`
*   **Description:** This script demonstrates how to generate homology heatmaps for *Schistosoma japonicum* based on homology results obtained from [**BioMart**](https://parasite.wormbase.org/biomart/martview/42c9211101263c9f47df0a8c88e13db3).
*   **Data:** The ortholog identity results for *Schistosoma japonicum* and other species are available on Zenodo.

### 2. m<sup>6</sup>A Peak Identification and Visualization
*   **Scripts:**
    *   `c2.1_m6a_peak.R`: Code for identifying m6A-enriched peaks based on BAM data.
    *   `c2.2_m6a_peak_visualization.R`: Visualization pipelines for the identified peaks.
*   **Features:**
    *   Metagene distribution plots
    *   Pie charts
    *   Motif analysis

## Dependencies & Data Availability

### Custom BSgenome Package
The standard *Schistosoma japonicum* BSgenome package is not currently available on Bioconductor. We have created a custom package based on the **WormBase reference (PRJNA520774)**.

### Downloads
Both the **BSgenome package** required for the m6A analysis and the **orthology identity results** for the ortholog analysis can be downloaded from Zenodo:

[**Download Data & Packages on Zenodo**](https://zenodo.org/records/18427709)

## Citation

If you use the code or data from this repository, please cite our paper (currently unpublished):

```bibtex
@article{SJ_analysis,
    title={m6A RNA modification controls reproductive development and metabolic adaptation in the human parasite Schistosoma japonicum},
    author={Bikash Ranjan Giri#, Zihao Zhang#, Lu Liu#, Chuantao Fang, Xiaoxu Wang, Shiyu Zhu, Mengting Qiu, Xiaoli Yan, Yan Ge, Zhuo Lan, Omri Wurtzel, Cizhong Jiang*, Guofeng Cheng*},
    journal={XX},
    year={2026},
    doi={xx}
}
```

## Contacts

If you encounter any issues with the code or have questions, please feel free to contact:

*   **Zihao Zhang:** 2110819@tongji.edu.cn
*   **Guofeng Cheng:** chengguofeng@tongji.edu.cn
