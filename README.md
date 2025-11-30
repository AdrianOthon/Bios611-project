Hi, this is my 611 Data Science  Project. More to come.

# Human Carotid Single-Cell Analysis  
## BIOS611 Final Project — Adrian Othon

This repository contains a fully reproducible single-cell RNA-seq analysis of 
human carotid endarterectomy tissue. The workflow demonstrates how to organize, 
containerize, and automate a complete data-science analysis using Docker, 
Make, and R.

All computational steps—including dimensionality reduction, clustering, marker 
identification, and figure generation—are performed inside a controlled Docker 
environment and orchestrated with a Makefile. The final report is built using 
R Markdown and can be reproduced from scratch with a single command.

---

## Content

bios611-project/
├── code/
│   └── carotid_analysis.R
├── data/
│   ├── adrian_smc_small.rds        # input Seurat object (ignored in git)
│   └── README.md                    # description & data source
├── figures/                         # automatically generated plots
├── results/                         # processed intermediate objects
├── report.Rmd                       # main report source
├── Dockerfile                       # reproducible R environment
├── makefile                         # project automation
└── README.md                        # this file


---

## Dataset

The input data derive from a publicly available single-cell RNA-seq dataset of 
human carotid plaque tissue (endarterectomy samples). The dataset includes:

- smooth muscle cells  
- endothelial cells  
- macrophages and monocytes  
- T cells and NK cells  
- additional stromal and immune populations  

A cleaned and preprocessed Seurat object is stored under `data/`. Because the 
file is large, it is **not included** in the GitHub repository; instead, `data/README.md`
provides instructions for obtaining it.

---

## Analysis Summary

The analysis script (`code/carotid_analysis.R`) performs:

- loading of the input Seurat object  
- normalization and feature selection  
- PCA computation  
- UMAP embedding  
- nearest-neighbor graph construction  
- Leiden clustering  
- QC visualization  
- differential expression to identify cluster markers  
- canonical marker validation  
- generation of all figures used in the final report

All output figures are saved to `figures/`. Intermediate objects are stored in `results/`.

The final R Markdown report (`report.Rmd`) integrates these figures with interpretive text.

---

## How to Build and Run the Project

### 1. Clone the repository

```bash
git clone <your-repo-url>
cd bios611-project

### 2. Build the Docker Image 

docker build -t bios611-project .

### 3. Run the container with RStudio Server 

docker run --rm -p 8787:8787 \
  -e PASSWORD=rstudio \
  -v "$(pwd)":/home/rstudio/bios611-project \
  bios611-project

Then open http://localhost:8787 in your browser

Login with username: rstudio, password: rstudio

### 4. Reproduce analysis and final report 

Inside the RStudio terminal run: 

cd ~/bios611-project
make

This will:

1. Run carotid_analysis.R
2. Generate all figures (figues/)
3. Knit the final report (report.html)

# Makefile Summary 

The makefile automates the entire workflow. Key targets:
1. make — runs everything
2. make figures — re-generates all analysis plots 
3. make report.html — rebuilds the R Markdown report
4. make clean — removes temporary files
5. make veryclean — removes all generated files (figures, results, report)

This ensures a transparent and dependency-tracked pipeline.

# Developer Instructions 

Developers can extend the workflow by:

1. Modifying or expanding code/carotid_analysis.R
2. Adding new analysis steps and figures
3. Updating the Makefile to track new dependencies
4. Modifying Dockerfile if additional R packages are needed

All computational steps should remain containerized so results are reproducible
regardless of host operating system or package versions.

# Contact 

adrian_othon@med.unc.edu
