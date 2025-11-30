#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

## paths and output dirs
seurat_in  <- "data/adrian_smc_small.rds"
seurat_out <- "data/adrian_smc_small_with_dr.rds"
stopifnot(file.exists(seurat_in))

dir.create("figures", showWarnings = FALSE, recursive = TRUE)
dir.create("results", showWarnings = FALSE, recursive = TRUE)

## load object
obj <- readRDS(seurat_in)
DefaultAssay(obj) <- "RNA"

## ensure PCA / UMAP / clusters
need_pca   <- !"pca"  %in% Reductions(obj)
need_umap  <- !"umap" %in% Reductions(obj)
need_clust <- is.null(obj$seurat_clusters)

if (need_pca || need_umap || need_clust) {
  obj <- NormalizeData(obj, verbose = FALSE)
  if (length(VariableFeatures(obj)) == 0) {
    obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
  }
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:30, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:30, verbose = FALSE)
  saveRDS(obj, seurat_out)
}

group_var <- if ("seurat_clusters" %in% colnames(obj@meta.data)) "seurat_clusters" else NULL

## FIGURE 1: PCA
stopifnot("pca" %in% Reductions(obj))
p_pca <- DimPlot(obj, reduction = "pca",
                 group.by = group_var, pt.size = 0.1) +
  ggtitle(if (is.null(group_var)) "PCA" else "PCA by cluster")
ggsave("figures/01_pca.png", p_pca, width = 6, height = 5, dpi = 300)

## FIGURE 2: UMAP
stopifnot("umap" %in% Reductions(obj))
p_umap <- DimPlot(obj, reduction = "umap",
                  group.by = group_var, label = TRUE, pt.size = 0.1) +
  ggtitle(if (is.null(group_var)) "UMAP" else "UMAP by cluster")
ggsave("figures/02_umap.png", p_umap, width = 6, height = 5, dpi = 300)

## FIGURE 3: QC violins (faceted)
qc_feats <- intersect(c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                      colnames(obj@meta.data))
if (length(qc_feats) > 0 && !is.null(group_var)) {
  qc_df <- FetchData(obj, vars = c(qc_feats, group_var))
  qc_df$cluster <- qc_df[[group_var]]
  qc_long <- qc_df |>
    pivot_longer(cols = all_of(qc_feats),
                 names_to = "metric", values_to = "value")
  p_qc <- ggplot(qc_long, aes(x = factor(cluster), y = value)) +
    geom_violin(fill = "grey85", color = "black", scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.12, outlier.size = 0.2) +
    facet_wrap(~ metric, scales = "free_y", nrow = 1) +
    labs(title = "QC metrics by cluster", x = "Cluster", y = "Value") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.minor = element_blank())
  ggsave("figures/03_qc_violins.png", p_qc, width = 11, height = 4, dpi = 300)
}

## FIGURE 4: marker heatmap (top 5 per cluster, balanced cells)
if (!is.null(group_var)) {
  obj <- JoinLayers(obj)
  markers <- FindAllMarkers(
    obj,
    only.pos = TRUE,
    logfc.threshold = 0.25,
    min.pct = 0.25,
    max.cells.per.ident = 500
  )
  if (nrow(markers) > 0) {
    saveRDS(markers, "results/all_markers.rds")
    top5 <- markers |>
      group_by(cluster) |>
      slice_max(order_by = avg_log2FC, n = 5)

    set.seed(123)
    cells_by_cluster <- split(colnames(obj), obj[[group_var, drop = TRUE]])
    n_per_cluster <- pmin(60, sapply(cells_by_cluster, length))
    cells_heat <- unlist(mapply(
      function(cells, n) if (n > 0) sample(cells, n) else character(0),
      cells_by_cluster, n_per_cluster,
      SIMPLIFY = FALSE
    ))

    p_heat <- DoHeatmap(
      obj,
      features = top5$gene,
      group.by = group_var,
      cells = cells_heat
    ) +
      NoLegend() +
      ggtitle("Top 5 markers per cluster")
    ggsave("figures/04_marker_heatmap.png", p_heat,
           width = 8, height = 10, dpi = 300)
  }
}

## FIGURE 5: canonical marker DotPlot
canonical <- intersect(
  c("ACTA2","TAGLN","PECAM1","VWF","CD68","LYZ","CD3E","NKG7"),
  rownames(obj)
)
if (length(canonical) > 0 && !is.null(group_var)) {
  p_dot <- DotPlot(obj, features = canonical, group.by = group_var) +
    RotatedAxis() +
    scale_color_gradient2(low = "navy", mid = "white",
                          high = "firebrick", midpoint = 0) +
    ggtitle("Canonical cell-type markers")
  ggsave("figures/05_canonical_markers_dotplot.png", p_dot,
         width = 8, height = 5, dpi = 300)
}

## FIGURE 6: cluster size barplot
if (!is.null(group_var)) {
  cl_counts <- table(obj[[group_var, drop = TRUE]])
  cl_df <- as.data.frame(cl_counts)
  colnames(cl_df) <- c("cluster", "count")
  p_sizes <- ggplot(cl_df, aes(x = cluster, y = count, fill = cluster)) +
    geom_col() +
    labs(title = "Cluster sizes", x = "Cluster", y = "Number of cells") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave("figures/06_cluster_summary.png", p_sizes,
         width = 8, height = 5, dpi = 300)
}
