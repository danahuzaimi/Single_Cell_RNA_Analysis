library(statmod)

# Estimate the dispersion to understand gene expression variability
# Dispersion is a measure of variability of gene expression counts around the mean. 
eDisp1 <- estimateDisp(dgeFiltNorm1, design1, robust = T)

# Taking the square root of common.dispersion
# common.dispersion is the average dispersion (variability) across all genes. 
sqrt(eDisp1$common.dispersion)
#[1] NA

# Voom is used to transform data to make it suitable for the linear modelling. 
v1 <- voom(eDisp1, design1, plot = TRUE)

# Adjust columns names 
colnames(design1) <-  gsub("Group_IL12RB", "", colnames(design1))

# Create a contrast matrix for comparisons between groups in the design matrix.
# Switch control and treatment 
contr.matrix1 <- makeContrasts(
  gene1_vs_none__dc_C = gene1_Dendritic_Tumor_Control - none_gene1_Dendritic_Tumor_Control,
  gene1_vs_none__dc_T = gene1_Dendritic_Tumor_NHS_rmIL12 - none_gene1_Dendritic_Tumor_NHS_rmIL12,
  gene1_vs_none__mono_C = gene1_Monocyte_Tumor_Control - none_gene1_Monocyte_Tumor_Control, 
  gene1_vs_none__mono_T = gene1_Monocyte_Tumor_NHS_rmIL12 - none_gene1_Monocyte_Tumor_NHS_rmIL12, 
  T_vs_C_dc_gene1 = gene1_Dendritic_Tumor_NHS_rmIL12 - gene1_Dendritic_Tumor_Control,
  T_vs_C_mono_gene1 = gene1_Monocyte_Tumor_NHS_rmIL12 - gene1_Monocyte_Tumor_Control,
  levels=design1)

#-------------------------------------------------------------------------------
# Git linear models for each gene.
vfit1 <- lmFit(v1, design1)

# Modify the linear models based on specific comparisons of interest.
vfit1 <- contrasts.fit(vfit1, contrasts = contr.matrix1)

# Use empirical Bayes methods to improve the reliability and power of the statistical tests.
efit1 <- eBayes(vfit1)

# Plot residual standard deviation versus average log expression for a fitted model. This plot assess the quality of the normalization process in RNA-seq data
# x_axis is avarage expression of each gene.
# y-axis is variability of gene expression from sample to sample. 
# Each dot is a gene. The blue line is the relationship between gene expression and variance.
plotSA(efit1)


# Find top-ranked genes from a linear model fit.
dt1 <- decideTests(efit1)
summary(dt1)
table(dt1)

# To find the upregulated genes 
# If dt1 is 1, it will return the corresponding row name from efit. 

ncol(contr.matrix1)

# There are 6 contrasts. we can look at up and down genes in each contrast

# Contrast 1
up_genes <- rownames(dt1)[dt1[,1] == 1]
down_genes <- rownames(dt1)[dt1[,1] == -1]

# complete the rest of the 6 contrasts following the above method
