unique_groups1 <- unique(md$gene1_Exp_Group)

# Aggregate data 
pb1 <- sapply(unique_groups1, function(x) {
    cell_groups1 <- md$Barcode[md$gene1_Exp_Group == x]
    
    # Subset countData to include only the columns (i.e cells) in that group.
    groupCountData1 <- countData[, cell_groups1]
    
    # Sum counts for each gene across all cells in the group.
    # This step creats pseudobulk data 
    summedCounts1 <- rowSums(groupCountData1)
    
    # The result is a matrix where each column is a unique group, and each row is a gene
    return(summedCounts1)
})

dim(pb1)
# [1] 19157    16

#To know the number of cells in each pb sample 
table(md$gene1_Exp_Group)

md_pb1 <- md %>%
  dplyr::select(
    Donor,
    Baseline,
    Treatment,
    Automated_cell_type_annotation_complete,
    gene1_Exp_Group,
    exp_status_RB,
    cell_group,
    Group_gene1
  ) %>%
  # Remove any duplicate row based on IL12RB_Exp_Group
  filter(!duplicated(gene1_Exp_Group))

# Set the row names from barcodes to the values in the Group column
row.names(md_pb1) <- md_pb1$gene1_Exp_Group

## Filter genes
library(edgeR)

all(row.names(md_pb1) == colnames(pb1))

# Create a DGEList object
# DGEList stands for "Differential Gene Expression List"
dge1 <- DGEList(counts = pb1, samples = md_pb1)

saveRDS(dge1, paste0(outPath, "DGEList_Pseudobulk_Hong_IL12_Mouse_IL12RB.RDS"))

# Design matrix modeling the difference between cell types from different subgroups, also models the variation from Treatment # "~ 0" means would not include an intercept column. 
design1 <- model.matrix(~ 0 + Group_IL12RB , data = md_pb1)

# Filter lowly expressed genes
keep.exprs1 <- filterByExpr(dge1, design = design1)

# Refine the dataset by including only  genes that are likely to provide meaningful information.
dgeFilt1 <- dge1[keep.exprs1, ] 

# Apply TMM normalisation
dgeFiltNorm1 <- calcNormFactors(dgeFilt1)

# Calculate log CPM(counts per million) 
# logCPM is a matrix with genes in rows and samples in columns
logCPM1 <- cpm(dgeFiltNorm1, log =T)
dim(logCPM1) 
# [1] 9877    8

# Visualize the distributions of counts and logCPM values. 
# The histogram shows gene expression levels across samples.
# The x-axis represent the levels of genes expression . 
# The y-axis: represent how many genes falls at each expression level. 
# Most genes have moderate activity levels. Only few genes with v.low or v.high activity level (0 vs 15)

hist(logCPM1, main="Histogram of logCPM1 Values", xlab="logCPM1")

#Boxplot to visualize the distribution of the logCPM values to see the effects of normalization.
# To adjust the margin of the plot 
png(
  paste0(figPath, dataName, "_", cellName, "_boxplot_logCPM1.png"),
  height = 5,
  width = 7,
  res = 300,
  unit = "in"
)

par(mar=c(12, 4, 4, 2), cex.axis=0.8)  # Decrease axis font size to 80% of the default
boxplot(logCPM1, las=2, col="lightblue", main="logCPM1 values across cells")

dev.off()


#-----------------------------------------------------------------------------

all(row.names(md_pb2) == colnames(pb2))

# Create a DGEList object
# DGEList stands for "Differential Gene Expression List"
dge2 <- DGEList(counts = pb2, samples = md_pb2)

saveRDS(dge2, paste0(outPath, "DGEList_Pseudobulk_Hong_IL12_Mouse_IL12A_B.RDS"))

# Design matrix modeling the difference between cell types from different subgroups, also models the variation from Treatment # "~ 0" means would not include an intercept column. 
design2 <- model.matrix(~ 0 + Group_IL12 , data = md_pb2)

# Filter lowly expressed genes
keep.exprs2 <- filterByExpr(dge2, design = design2)

# Refine the dataset by including only  genes that are likely to provide meaningful information.
dgeFilt2 <- dge2[keep.exprs2, ] 

# Apply TMM normalisation
dgeFiltNorm2 <- calcNormFactors(dgeFilt2)

# Calculate log CPM(counts per million) 
# logCPM is a matrix with genes in rows and samples in columns
logCPM2 <- cpm(dgeFiltNorm2, log =T)
dim(logCPM2) 
# [1] 9877    8

# Visualize the distributions of counts and logCPM values. 
# The histogram shows gene expression levels across samples.
# The x-axis represent the levels of genes expression . 
# The y-axis: represent how many genes falls at each expression level. 
# Most genes have moderate activity levels. Only few genes with v.low or v.high activity level (0 vs 15)

hist(logCPM2, main="Histogram of logCPM2 Values", xlab="logCPM2")

#Boxplot to visualize the distribution of the logCPM values to see the effects of normalization.
# To adjust the margin of the plot 
png(
  paste0(figPath, dataName, "_", cellName, "_boxplot_logCPM2.png"),
  height = 5,
  width = 7,
  res = 300,
  unit = "in"
)

par(mar=c(12, 4, 4, 2), cex.axis=0.8)  # Decrease axis font size to 80% of the default
boxplot(logCPM2, las=2, col="lightblue", main="logCPM2 values across cells")

dev.off()
