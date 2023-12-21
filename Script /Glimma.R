library(Glimma)

cell_surface <- read.delim("Path_to_Surface_annotations.tsv") %>% 
  dplyr::rename(SYMBOL = hgnc)

#--------------------------------------------------------

# Save each top table and run GSEA analysis on each contrast
for(i in 1:length(colnames(contr.matrix1))){
  
  # Set the grouping factor for glimma based on the contrast in question
  contrast1 <- colnames(contr.matrix1)[i]
 
  MA_fac1 <- factor(md_pb1$exp_status_RB, levels = unique(md_pb1$exp_status_RB)[order(unique(md_pb1$exp_status_RB))])
  
  # #Save the glimma plots
  htmlwidgets::saveWidget(glimmaMA(efit1, coef = contrast1, main = gsub("_", " ",contrast1), dge1 = dgeFiltNorm1, groups = MA_fac1),
                          paste0(MA, contrast1, "_RNA_MA.html"))
  
  # Rounding is a hack to make expression work
  htmlwidgets::saveWidget(glimmaVolcano(efit1, coef = contrast1, main = gsub("_"," ",contrast1), dge1 = dgeFiltNorm1, groups = MA_fac1),
                          paste0(VOL, contrast1, "_RNA_Volcano.html"))
}

cont1 <- colnames(contr.matrix1)

##----- Extract DE stats and DEGs:
DEStatList_limma<- list()

for (cc in cont){
  DEStatList_limma[[cc]] <- topTable(efit1, n = Inf, coef = cc)

saveRDS(DEStatList_limma_RB, paste0(outPath, "DEStatList_limma_RB.RDS"))
