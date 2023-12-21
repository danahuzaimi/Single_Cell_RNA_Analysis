# Normalise data to account for differences in sequencing depth between cells 
m <- NormalizeData(m)

# Detection of variable genes across the single cells
# vst = Variance Stabilizing Transformation
# Selects top 2000 variable features
m <- FindVariableFeatures(m, selection.method = "vst", nfeatures = 2000)

# Show the top 20 genes 
top10 <- head(VariableFeatures(m), 20)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(m)
plot2 <- LabelPoints(plot = plot1,
                     points = top10,
                     repel = TRUE)

# HVGs = Highly Variable Genes
png(
  paste0(figPath, dataName, "_", cellName, "_HVGs_vst_logNorm.png"),
  height = 5,
  width = 7,
  res = 300,
  unit = "in"
)
plot2
dev.off()
