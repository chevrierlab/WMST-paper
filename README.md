# WMST-paper

Code repository associated with the manuscript  
**“Molecular and cellular cartography of the laboratory mouse using whole-body sections.”**  
Manuscript currently under review.

This repository serves as a snapshot of the analysis code used in the study.

## Repository structure

- **`LABEL/`**  
  Code for the **LABEL pipeline** (Spatially Aware Histology-Based Classification), used for histology-based tissue, tissue region, and cell type annotation.

- **`analysis/spatial transcriptomics/`**  
  Analysis scripts and notebooks for whole-body spatial transcriptomics, including preprocessing, subregion annotation, differential expression, and downstream analyses used to generate figures and tables in the manuscript.


## Computing environment

All analyses were performed using:
- **R 4.3.3**
- **Python** (Scanpy/AnnData-based workflows)
- **Command-line tools** for data handling and preprocessing where appropriate

## Data availability

Processed data objects and large image files are not stored directly in this repository due to size constraints. Data availability and access details are described in the manuscript. LABEL can be accessed for reproducible runs on Code Ocean at https://doi.org/10.24433/CO.3645054.v1

## Notes

- For questions related to the code or analyses, please contact the corresponding authors.