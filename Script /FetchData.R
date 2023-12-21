Idents(m) = "seurat_clusters"
table(m@meta.data$Automated_cell_type_annotation_complete, m@meta.data$seurat_clusters)

#gene1
gene1 <- FetchData(m, vars = "gene1")
cells_expressing_gene1<- sum(gene1 > 0)
print(paste("Number of cells expressing gene1:", cells_expressing_gene1))

#gene2
gene2 <- FetchData(m, vars = "gene2")
cells_expressing_gene2 <- sum(gene2 > 0)
print(paste("Number of cells expressing gene2:", cells_expressing_gene2))

# Co-expression
gene12 <- FetchData(m, vars = c("gene1", "gene2"))
coexpression_gene12 <- gene12 > 0 
coexpression_gene12 <- as.data.frame(coexpression_IL12)
coexpression_gene12 <- mutate(coexpression_gene12, Sum = rowSums(coexpression_gene12))
#"Number 2 mean coexpression"
table(coexpression_gene12$Sum)
