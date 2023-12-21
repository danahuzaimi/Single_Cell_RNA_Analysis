#-----------RLE--------------
# RLE is used to check for batch effects or outliers.

RLE1 <- t(apply(logCPM1, 1, function(x) x - median(x)))

# Create a boxplot
par(mar=c(12, 4, 4, 2), cex.axis=0.8)
boxplot(RLE1, outline = FALSE, las = 2, cex.axis = 0.4)

#-----------MDS--------------

# MDS is used to visualize the distances or dissimilarities between samples. 
# Use the top 500 genes with the highest variability for the plot
plotMDS(logCPM1, top=500, labels=colnames(logCPM1), col=c("red","blue"))
                
#-----------PCA--------------
# PCA is used to emphasize variation and bring out strong patterns in a dataset (dimensionality reduction). 
# # Perform Singular Value Decomposition (SVD) 
Annotation1 <- dgeFiltNorm1$samples
s<- svd(apply(logCPM1 , 1, function(x) scale(x, scale=FALSE , center = TRUE)))

# Add new coll'SampleID'.
Annotation1 <- dgeFiltNorm1$samples
Annotation1$SampleID<- rownames(Annotation1)
source("path_to_function")

PCA_plot (expr = logCPM1,
  clin = Annotation1,
  group_to_test = 'gene1_Exp_Group',
  svdResult = s,
  plot_cex = 1,
  legend_cex = 1,
  labeled = FALSE,
  group_to_shape = NULL,
  cols =
    c(brewer.pal(8, "Dark2")[-5],
      brewer.pal(10, "Paired"),
      brewer.pal(12, "Set3"))
  ) 
#-----------# Library size assesment -------------

dge_sample1 <- dgeFiltNorm1$samples
lib_size1 <- dge_sample1[c("lib.size", "gene1_Exp_Group")]

# Find median values after applying RLE normalization 
median1 <- (apply(RLE1, 2, median))
median1 <- as.data.frame(median1)
median1$gene1_Exp_Group = rownames(median1) 

# Merge data 
med_lib1 <- merge(median1, lib_size1, by = "gene1_Exp_Group")
# Line plot 
Median_values1 <- med_lib1$median
library_size1 <- med_lib1$lib.size
plot(Median_values1, library_size1, type = "p", main = "median vs library size in gene1 cells")
