
# Heatmap of Percentage of cells per cell type based on expression of Gene1 and Gene2 
# Calculates the percent of cells that express a given set of features by various grouping factors
#install.packages("scCustomize")
#library(scCustomize)

# Percentage of cells that express 'Gene1'and 'Gene2' from each donor across the entire dataset 
percent_by_donor <- Percent_Expressing(
  m,
  features = c('Gene1', 'Gene2'),
  threshold = 0,
  split_by = "Donor",
  entire_object = TRUE,
  slot = "data")


# Percentage of cells expressing genes of interest by ident= Automated_cell_type_annotation_complete 

percent1 <- Percent_Expressing(
  m,
  features = c('Gene1', 'Gene2'),
  threshold = 0,
  entire_object = FALSE,
  slot = "data")


# Convert df to a matrix
percent1_matrix <- as.matrix(percent1)

#Heatmap 
pheatmap(percent1_matrix, angle_col = 45, legend = TRUE, main = "% Percentage of cells expressing Gene1 and 'Gene2, display_numbers = T)
