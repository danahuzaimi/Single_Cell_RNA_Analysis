```r
library(sctransform)

m <- SCTransform(
  object = m,
  assay = "RNA",
  new.assay.name = "SCT",
  do.correct.umi = TRUE,
  ncells = NULL,
  variable.features.n = 3000,
  # variable.features.rv.th = 1.3,
  vars.to.regress = c("percent.mt", "G2M.Score", "S.Score"),
  do.scale = FALSE,
  do.center = TRUE,
  conserve.memory = FALSE,
  return.only.var.genes = TRUE,
  seed.use = 20211101,
  verbose = TRUE
)

top10 <- head(VariableFeatures(m), 20)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(m)
plot2 <- LabelPoints(plot = plot1,
                     points = top10,
                     repel = TRUE)

png(
  paste0(figPath, dataName, "_", cellName, "_HVGs_vst_SCT.png"),
  height = 5,
  width = 7,
  res = 300,
  unit = "in"
)
print(plot2)
dev.off()
```
