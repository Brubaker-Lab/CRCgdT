library(slingshot)
library(SingleCellExperiment)
library(grDevices)
library(dplyr)
library(ggrepel)
library(ggplot2)
library(viridisLite)
setwd()
gdT = readRDS()
sce <- SingleCellExperiment(
  assays = list(counts = gdT[["SCT"]]@counts),  # Add the count matrix to the assay slot
  colData = DataFrame(cluster = gdT@meta.data$cell.type)  # Add the clustering labels to colData
)
reducedDims(sce) <- SimpleList(PCA =  gdT@reductions$pca@cell.embeddings)
start <- slingshot(sce, clusterLabels = 'cluster', reducedDim = 'PCA', start.clus ='Tpex')
saveRDS(start, 'pseudotime_start_Tpex.rds')
my_colors = c( 'Poised Teff 1' = "#A7E1F1",'Poised Teff 2' = "#4c9bd4",'IL7R+ TRM' = "#3B749C",'IL7R- TRM' = "#1452a3",'Tpex' = '#6E51E5','Tex1' = "#132876",'Tex2' = '#4539A3')
curves <- slingCurves(start)
lineages <- slingLineages(start)
clustering = gdT@meta.data$cell.type

pca_data <- as.data.frame(reducedDims(start)$PCA)
pca_data$CellType <- factor(clustering)

centroids <- pca_data %>%
  group_by(CellType) %>%
  summarize(
    Centroid_PC1 = mean(PC_1),  # Centroid for PC1
    Centroid_PC2 = mean(PC_2)   # Centroid for PC2
  )

lincol = c('red', 'orange', 'yellow')
linwd = c(6,3,2)
cexs = c (2,1.5,1)
p <- ggplot(pca_data, aes(x = PC_1, y = PC_2)) +
  geom_point(aes(color = CellType), size = 2) +
  scale_color_manual(values = my_colors)

for (i in seq_along(lineages)) {
  lineage <- lineages[[i]]
  
  centroids_lineage <- centroids %>%
    filter(CellType %in% lineage) %>%
    mutate(CellType = factor(CellType, levels = lineage)) %>%  # Ensure CellType follows lineage order
    arrange(CellType)
  
  # Add centroids and connections to the plot
  p <- p + 
    geom_point(data = centroids_lineage, aes(x = Centroid_PC1, y = Centroid_PC2), 
               color = lincol[i], size = cexs[i]*2) +  # Manually assign color and size outside of aes()
    geom_path(data = centroids_lineage, aes(x = Centroid_PC1, y = Centroid_PC2, group = 1), 
              color = lincol[i], size = cexs[i])  # Connect centroids with a line
}

p = p + geom_text_repel(data = centroids, color ="white", aes(Centroid_PC1 , Centroid_PC2, label = CellType),
                        size = 4,
                        fontface = 'bold',
                        bg.color = 'black')
# Print the plot
print(p)
ggsave('Pseudotime.png', plot = p, width = 8, height = 6.5)

start = readRDS('pseudotime_gdT.rds')
library(tradeSeq)
cellweights = slingCurveWeights(start)
pt = slingPseudotime(start)
for  (i in seq_along(lineages)) {
  pseudotime = pt[,i]
  cellweight = cellweights[,i]
  pseudotime = pseudotime[cellweight>0]
  counts = gdT[["SCT"]]@counts[row.names(gdT@assays$integrated),cellweight>0]
  metadata = gdT@meta.data[colnames(counts),]['cell.type']
  cellweight = cellweight[cellweight>0]
  gam.pval <- apply(counts,1,function(z){
    d <- data.frame(z=z, t=pseudotime)
    tmp <- gam(z ~ lo(t), data=d)
    pse <- summary(tmp)[4][[1]][1,5]
    pse
  })
  saveRDS(gam.pval,paste0('lin',as.character(i),'.rds'))
  #Select the top 50 most significant genes that change over pseudotime
  
  gam.pval_sig = gam.pval[gam.pval<0.05]
  topgenes <- c(names(sort(gam.pval_sig, decreasing = FALSE)))
  topgenes = topgenes[topgenes != ""]
  topgenes = topgenes[1:50]
  
  pst.ord <- order(pseudotime)
  
  heatdata <- counts[topgenes,pst.ord]
  heatdata <- t(scale(t(heatdata)))
  
  heatdata[heatdata > 3] = 3
  heatdata[heatdata < -2] = -2
  
  Cell.Types_df <- data.frame(Cell.Types = metadata[pst.ord,])
  
  Cell.Types_df$Cell.Types <- factor(Cell.Types_df$Cell.Types, levels = names(my_colors))
  Cell.Types_df$pseudotime <- sort(pseudotime)
  row.names(Cell.Types_df) = row.names(metadata)[pst.ord]
  
  pseudotime_colors <- colorRampPalette(c("blue", "white", "red"))
  breaks <- seq(from = min(Cell.Types_df$pseudotime, na.rm = TRUE), 
                to = max(Cell.Types_df$pseudotime, na.rm = TRUE), 
                length.out = 100)
  annotation_colors = list(
    Cell.Types = my_colors,
    pseudotime = pseudotime_colors(length(breaks) - 1)
  )
  
  color_gradient <- viridis(100)

  ptime_heatmap = pheatmap(heatdata,
                           annotation_col = Cell.Types_df, 
                           annotation_colors = annotation_colors, 
                           cluster_rows = TRUE, 
                           cluster_cols = FALSE, 
                           show_rownames = TRUE, 
                           show_colnames = FALSE,
                           border_color = NA,
                           color = color_gradient, 
                           treeheight_row = 10,
                           width = 7, height = 8, fontsize = 8,fontsize_col = 12) 
  
  ggsave(paste0('ptime_map',as.character(i),'.png'), plot = ptime_heatmap, width = 4, height = 6,dpi = 300)
}
