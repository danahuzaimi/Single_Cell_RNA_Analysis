# Shows the slots in seurat obj
str(m)

# Access one of the slot in seurate obj
dim(m@meta.data)

# Access a variable in the seurat obj 
#m$Treatment

# Return the number of genes (rows) and number of cells (cols)
dim(m@assays$RNA@counts)

# To see cells types in each sample 
table(m@meta.data$Automated_cell_type_annotation_complete, m@meta.data$Sample_name)
table(m@meta.data$Automated_cell_type_annotation_complete, m@meta.data$Treatment)
