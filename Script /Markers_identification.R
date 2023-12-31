# Retrive the cell identity 
head(Idents(object = m))
# The table shows there are 10 clusters (0 - 9) Numbers represent the number of cells in each cluster. 
table(Idents(m))


Idents(m) = 'Treatment'
# The table shows number of cells for each treatment group. 
table(Idents(m))

# To find the number of clusters 
cluster_ids <- Idents(m)
number_of_clusters <- length(unique(cluster_ids))

cat("Number of clusters:", number_of_clusters)

head(m$seurat_clusters)
#FindAllMarkers(m)

table(m@meta.data$seurat_clusters, m@meta.data$Treatment)

# Only positive markers are considered 
# 0.25 is the min percentage threshold for a marker to be considered
# Select top 2 markers with the highest avg_Log2FC from each cluster 
m.markers <- FindAllMarkers(m, only.pos = TRUE, min.pct = 0.25,logfc.threshold = 0.25)
m.markers %>% group_by(cluster) %>% slice_max( n = 2, order_by = avg_log2FC) 

# Find all markers that distinguish cluster 0 from cluster 8
#cluster0markers = FindMarkers(m, ident.1 = 0, ident.2 = 8 , min.pct = 0.25)
#head(cluster0markers, n = 3)
# pct.1 and pct.2 represent the percentage of cells expressing a particular gene within two clusters being compared.
```
```r 
# Find genes that are differential expressed between the two treatment groups
Idents(m) = 'Treatment'
# Finding differential expression markers
DE_Treatment <- FindMarkers(m, ident.1 = "NHS-rmIL12 treated", ident.2 = "PBS treated", min.pct = 0.25)

# Subset the top 5 significant genes
top5_Selection <- top_n(DE_Treatment, -5, p_val_adj)
# Printing the top 5 rows
print(top5_Selection)

# Create a volcano plot
# Assign color of each point based on avg_log2FC and significance
colors_points <- function(log2FC, p_value) {
  # Changed the significant value to 1x10^(-100) to only visualize genes that are highly significant on the plot 
  significant <- -log10(1e-100) 
  if(p_value < significant) 
    {return("gray")} 
  else 
    if(log2FC < 0) {return("blue")} # Down
  else {return("red")} # Up
}

# Determine colors for each point in the dataset
DE_Treatment$colors <- mapply(colors_points, DE_Treatment$avg_log2FC, -log10(DE_Treatment$p_val_adj))

# Identify genes to label: genes with colors either blue or red
genes_label <- DE_Treatment[DE_Treatment$colors %in% c("blue", "red"), ]
genes_label

# Plot
volcano_plot <- ggplot(DE_Treatment, aes(x = avg_log2FC, y = -log10(p_val), color = colors)) +  # Use the 'colors' column for the color aesthetic
  geom_point(alpha = 0.5) + 
  scale_color_identity() +  #  tells ggplot to use the colors we've specified directly
  scale_color_manual(name = "Expression",
                     values = c("gray" = "gray", "blue" = "blue", "red" = "red"),
                     breaks = c("blue", "red", "gray"),
                     labels = c("Down", "Up", "Not Significant")) +
  labs(
    title = "Volcano Plot of Differentially Expressed Genes Between Treatment Groups",
    x = "avg_log2FC",
    y = "-log10(p-value)"
  ) +
  theme(legend.position = "right") + 
  theme_bw(9) +  
  geom_hline(yintercept = -log10(1e-100), linetype = "dashed", color = "black") +  
  geom_vline(xintercept = 0 , linetype = "dashed", color = "black") + 
  geom_text(data = genes_label, aes(label = rownames(genes_label)), vjust = -1, hjust = 1, size = 2)

#Show the plot 
volcano_plot

# Visualizes the expression of the selected genes using VlnPlot

png(
  paste0(
    figPath,
    dataName,
    "_",
    cellName,
    "_Violin_Plot_DE_by_Treatment_Top5.png"
  ),
  height = 15,
  width = 28,
  res = 300,
  unit = "in"
)
# Access the row indices of top5_Selection
row_indices <- row.names(top5_Selection)

# Print the row indices
print(row_indices)

VlnPlot(m, features = row_indices, group.by = "Treatment", assay = "RNA", pt.size = 0.1, ncol = 5)

dev.off()

# Visualizes the expression of the selected genes using Heatmap

png(
  paste0(
    figPath,
    dataName,
    "_",
    cellName,
    "_Heatmap_DE_by_Treatment_Top5.png"
  ),
  height = 15,
  width = 28,
  res = 300,
  unit = "in"
)
# Access the row indices of top5_Selection
row_indices <- row.names(top5_Selection)

# Heatmap using row_indices as features
DoHeatmap(m, features = row_indices, size = 3)

dev.off()

```
```r
# Task: Group cells based on cell type x donor to know the percentage of each cells from each donor 
# Make a table 
cells_dt <- table(m@meta.data$Donor, m@meta.data$Automated_cell_type_annotation_complete)
# Find the proportion of each cell type per row 
cells_d <- prop.table(cells_dt, 1)*100
# add sum to the table and show percentage 
addmargins(cells_d)

# Another way
# find sum of all cells 
total_cell <- sum(cells_dt)

# Divide each element in the table by total * 100 
per_cell <- (cells_dt/total_cell) * 100

# Convert table to df 
per_cell <- data.frame(rbind(per_cell))

# add sum
# Percentage of cells from each donor. total is 100
row1 <- print(rowSums(per_cell))
