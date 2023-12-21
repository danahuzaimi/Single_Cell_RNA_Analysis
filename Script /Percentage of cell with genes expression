#install.packages("scCustomize")
library(scCustomize)

percent1 <- Percent_Expressing(
  m,
  features = c('gene1', 'gene2'),
  group_by = 'Automated_cell_type_annotation_complete', 
  threshold = 0,
  entire_object = FALSE)


# Convert df to a matrix
percent1_matrix <- as.matrix(percent1)

#Heatmap 
# Define a color palette from light to dark
color_palette <- colorRampPalette(c("lightyellow", "red"))(n = 100)  # n defines the number of colors

# Create the heatmap with the custom color palette
#install.packages("pheatmap")
library(pheatmap)

pheatmap(percent1_matrix, 
         angle_col = 45, 
         legend = TRUE, 
         main = "% Percentage of cells expressing IL12RB1 and IL12RB2", 
         display_numbers = T, 
         color = color_palette)
# -----------------------------------------
percent2 <- Percent_Expressing(
  m,
  features = c('gene1', 'gene2'),
  group_by = 'Treatment', 
  threshold = 0,
  entire_object = FALSE)


# Convert df to a matrix
percent2_matrix <- as.matrix(percent2)

#Heatmap 
# Define a color palette from light to dark
color_palette <- colorRampPalette(c("lightyellow", "red"))(n = 100)  # n defines the number of colors

# Create the heatmap with the custom color palette
#install.packages("pheatmap")
library(pheatmap)

pheatmap(percent2_matrix, 
         angle_col = 45, 
         legend = TRUE, 
         main = "% Percentage of cells expressing IL12RB1 and IL12RB2", 
         display_numbers = T, 
         color = color_palette)
