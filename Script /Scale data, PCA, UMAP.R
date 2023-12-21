m <-
  # Principal Component Analysis (PCA) reduces the dimensionality of the data  
  RunPCA(m, verbose = F)

ElbowPlot(m, ndims = 50)


# Find the nearest neighbor for each cell based on the first 40 dimensions 
m <- FindNeighbors(m, dims = 1:40)

# Clustering
# Resolution controls the granularity of clustering. Higher values = more clusters, while lower values group cells together more broadly
m <- FindClusters(m, 
                    resolution = 0.5, 
                    random.seed = 20220613)

# UMAP (Uniform Manifold Approximation and Projection) To reduce dimensionality 
m <-
  RunUMAP(
    m,
    dims = 1:40,
    n.neighbors = 50,
    min.dist = 0.3,
    seed.use = 20220613
  )



annCols <-
  c(
    "seurat_clusters",
    "Treatment", 
    "Sample_name",
    "Automated_cell_type_annotation_complete",
    "sc_type_classification",
    "Phase"
  )

for (a in annCols) {
  png(
  paste0(figPath, dataName, "_", cellName,
         "_DimPlot_UMAPs_Annot_", a, ".png"),
  height = 7.5,
  width = 5, 
  unit = "in", 
  res = 300
)
  if (a == "seurat_clusters") {
    print(
      DimPlot(
        m,
        group.by = a,
        reduction = "umap",
        label = T
      ) +
        scale_color_manual(values = cols) +
        theme_bw() +
        theme(legend.position = "top",
              legend.direction = "vertical") +
        guides(color = guide_legend(
          nrow = 8,
          override.aes = list(size = 3)
        ))
    )
  } else{
    print(
      DimPlot(m, group.by = a, reduction = "umap") +
        scale_color_manual(values = cols) +
        theme_bw() +
        theme(legend.position = "top",
              legend.direction = "vertical") +
        guides(color = guide_legend(
          nrow = 8,
          override.aes = list(size = 3)
        ))
    )
  }
  dev.off()
}
# elbow" point in the plot is where variance starts to level off. This indicates the optimal number of PCs to use for dimensionality reduction.
