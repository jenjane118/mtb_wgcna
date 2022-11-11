# make new heatmap for trait/module correlation

# either grey out non-significant correlations or select for certain traits/modules

library(WGCNA)

# load trait_df
trait_df <- readRDS(here("R_data/trait_df.RData"))
# load MEs, module colours, etc
load(here("R_data/TB_modules_network_construction.RData"))


moduleTraitBicor.data <- bicorAndPvalue(MEs,
                                        trait_df,
                                        maxPOutliers=0.05,
                                        robustY = F)
# create dataframe of module bicorrelation with trait
moduleTraitBicor <- moduleTraitBicor.data$bicor
# create dataframe of p-values for bicor 
moduleTraitBicorPvalue <- as.data.frame(moduleTraitBicor.data$p)

# fdr for multiple testing applied for each trait
traitNames<-colnames(trait_df)
modNames <- substring(colnames(MEs),3)
# create dataframe for p_adj
p_adj_bicor<-data.frame(matrix(0, nrow=length(modNames),
                               ncol=length(traitNames)), 
                        row.names = modNames)
colnames(p_adj_bicor) <- traitNames
# for each trait find adjusted p-value
for (i in 1:length(traitNames)){
  p_adj_bicor[,i]<-p.adjust(moduleTraitBicorPvalue[,i],method="fdr")
}
#convert df to matrix
p_adj_moduleTraitBicorPvalue <- as.matrix(p_adj_bicor)

# make heatmap of module/trait correlations
# create text matrix for what text appears in each box
# this includes bicor correlation calc / adjusted p-value
textMatrix <- paste(signif(moduleTraitBicor, 2),
                    "\n(",
                    signif(p_adj_moduleTraitBicorPvalue, 1),
                    ")",
                    sep = "")
dim(textMatrix) <- dim(moduleTraitBicor)

#create color matrix for what color appears in each box
# mutate case_when signif > 0.05, color == grey
#colMatrix <- 
#make new moduleTraitBicor where if corresponding pvalue in p_adj_bicor is >0.05, label NA

new_moduleTraitBicor <- moduleTraitBicor
for (i in 1:length(colnames(moduleTraitBicor))){
  index <- p_adj_bicor[,i] > 0.05
  new_moduleTraitBicor[index,i] <- NA
}


# #heatmap plot with greyed out non-significant boxes
#png(here("Images/labeled_heatmap_grey.png"), width = 2400, height = 1920)
par(mar = c(8, 12, 2, 2))
labeledHeatmap(Matrix <- new_moduleTraitBicor,
               xLabels = colnames(trait_df),
               yLabels = colnames(MEs),
               ySymbols = colnames(MEs),
               colorLabels = TRUE,
               colors = blueWhiteRed(50),  #don't use if using colorMatrix
               #showRows = c(),  # numeric vector with indices of rows to be shown
               #showCols = c(),  # numeric vector with indices of cols to be shown
               naColor="grey",
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1,
               cex.lab.x = 1.5,
               cex.lab.y = 0.75,
               zlim = c(-1,1),
               main = paste("Module-trait relationships FDR adjusted p-values: net"))
#dev.off()


# greyed out non-sig boxes AND eliminating conditions with no significant correlations

my_cols <- c(4,5,6,7,12,14,15)
#png(here("Images/heatmap_grey_selected.png"), width = 1200, height = 1920)
par(mar = c(8, 12, 2, 2))
labeledHeatmap(Matrix <- new_moduleTraitBicor,
               xLabels = colnames(trait_df),
               yLabels = colnames(MEs),
               ySymbols = colnames(MEs),
               colorLabels = TRUE,
               colors = blueWhiteRed(50),  #don't use if using colorMatrix
               #showRows = c(),  # numeric vector with indices of rows to be shown
               showCols = my_cols,  # numeric vector with indices of cols to be shown
               naColor="grey",
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1,
               cex.lab.x = 1.5,
               cex.lab.y = 0.75,
               zlim = c(-1,1),
               main = paste("Module-trait relationships and FDR adjusted p-values for selected conditions"))
#dev.off()


# use plotDendroAndColors() to plot heatmap of traits with dendrogram on side?

sum_conn<-signif(cor(MEs, use="p"), 2)
# use this to cluster the eigengenes
dissimME=(1-t(cor(MEs, method="p")))/2 
hclustdatME=hclust(as.dist(dissimME), method="average" )
dend <- as.dendrogram(hclustdatME)

#heatmap with dendrogram
#png(here("Images/network_heatmap_dendro.png"), width = 2400, height = 2200)
par(mar = c(1, 1, 2, 1), mai = c(1.02, 0.82, 0.82, 0.42), cex.main=4)
heatmap(Matrix <- moduleTraitBicor,
        Rowv = dend,
        Colv = NA,
        labCol = colnames(trait_df),
        labRow = colnames(MEs),
        col = blueWhiteRed(50),
        RowSideColors = sub("ME", "", colnames(MEs)),
        cexCol = 3,
        cexRow = 3,
        margins  = c(20,30),
        main = "Module-trait relationships")
#dev.off()


# with selected cols and white for non-sig values
#png(here("Images/network_heatmap_dendro_select.png"), width = 2400, height = 2200)
par(bg="lightgrey")
par(mar = c(1, 1, 2, 1), mai = c(1.02, 0.82, 0.82, 0.42), cex.main=4)
heatmap(Matrix <- new_moduleTraitBicor[,my_cols],
        Rowv = dend,
        Colv = NA,
        labCol = colnames(trait_df[,my_cols]),
        labRow = colnames(MEs),
        col = blueWhiteRed(50),
        RowSideColors = sub("ME", "", colnames(MEs)),
        cexCol = 3,
        cexRow = 3,
        margins  = c(20,30),
        na.rm = T,
        main = "Module-trait relationships")
#dev.off()


