# Default target
all: report.html

## -----------------------
## Paths
## -----------------------
SEURAT_OBJ = data/adrian_smc_small.rds
ANALYSIS   = code/carotid_analysis.R
REPORT_RMD = report.Rmd

FIGS = \
  figures/01_pca.png \
  figures/02_umap.png \
  figures/03_qc_violins.png \
  figures/04_marker_heatmap.png \
  figures/05_canonical_markers_dotplot.png \
  figures/06_cluster_summary.png

## -----------------------
## Build figures
## -----------------------
$(FIGS): $(SEURAT_OBJ) $(ANALYSIS)
	Rscript $(ANALYSIS)

figures: $(FIGS)

## -----------------------
## Build report
## -----------------------
report.html: $(REPORT_RMD) $(FIGS)
	Rscript -e "rmarkdown::render('$(REPORT_RMD)', output_format = 'html_document')"

## -----------------------
## Cleaning
## -----------------------
.PHONY: clean veryclean figures all

clean:
	rm -f report.html

veryclean: clean
	rm -f figures/*.png results/*.rds