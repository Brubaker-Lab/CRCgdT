main_dir <- ""

folders <- list.dirs(path = main_dir, full.names = TRUE, recursive = FALSE)
gsm_folders <- folders[grep("^GSM", basename(folders))] #only the 36-39 in base folder
seurat_objects = list()
for (folder in gsm_folders) {
    data <- Read10X(data.dir = folder)
    seurat_object <- CreateSeuratObject(counts = data, project = basename(folder))
    seurat_objects[[basename(folder)]] <- seurat_object
}


for (i in seq_along(seurat_objects)) {
  mito.genes <- grep("^MT-", rownames(seurat_objects[[i]]), value = TRUE)
  seurat_objects[[i]][["percent.mt"]] <- PercentageFeatureSet(seurat_objects[[i]], pattern = "^MT-")
  seurat_objects[[i]] <- subset(seurat_objects[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA > 500 & percent.mt < 20)
}
seurat_objects <- lapply(seurat_objects, function(x) {
  x <- NormalizeData(x, verbose = FALSE)
})
seurat_objects <- lapply(seurat_objects, function(x) {
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})
seurat_objects <- lapply(seurat_objects, function(x) {
  x <- ScaleData(x)
})
seurat_objects <- lapply(seurat_objects, function(x) {
  x <- RunPCA(x, verbose = FALSE)
})

features_to_integrate <- SelectIntegrationFeatures(object.list = seurat_objects, nfeatures = 3000)

anchors <- FindIntegrationAnchors(object.list = seurat_objects, normalization.method = "LogNormalize", anchor.features = features_to_integrate)
integrated_data <- IntegrateData(anchorset = anchors, normalization.method = "LogNormalize")

integrated_data <-ScaleData(integrated_data)
integrated_data <- RunPCA(integrated_data)
integrated_data <- FindNeighbors(integrated_data, dims = 1:50)
integrated_data <- RunUMAP(integrated_data, dims = 1:50)

integrated_data <- FindClusters(integrated_data, resolution = 0.5)
DimPlot(integrated_data, reduction = "umap", group.by = "seurat_clusters")
FeaturePlot(integrated_data, features = c('CD3G','CD3E','CD3D','TRDC','ZBTB16','IKZF2'), reduction = 'umap', min.cutoff = 0, max.cutoff = 5)




saveRDS(seurat_objects, 'seurat_objects.rds')
seurat_objects = readRDS('seurat_objects.rds')
saveRDS(anchors, 'anchors.rds')
anchors = readRDS('anchors.rds')
saveRDS(integrated_data, 'integrated_data.rds')
integrated_data = readRDS('integrated_data.rds')


main_dir <- ""
folders <- list.dirs(path = main_dir, full.names = TRUE, recursive = FALSE)
gsm_folders <- folders[grep("^GSM", basename(folders))]
cond_names = list()
for (folder in gsm_folders) {
  cond_names = c(cond_names,basename(folder)) 
}

rna_assay <- integrated_data@assays$RNA
counts_slots <- names(rna_assay@layers)[grepl("^counts", names(rna_assay@layers))]
cell_counts <- sapply(counts_slots, function(slot) ncol(rna_assay@layers[[slot]]))
names(cell_counts) = cond_names
cell_origins <- rep(names(cell_counts), times = cell_counts)
integrated_data@meta.data$orig.ident <- cell_origins

integrated_data = readRDS('integrated_data.rds')
barcode_data <- read.csv("GSE221575_gd.csv")
merged = JoinLayers(integrated_data@assays$RNA)
wanted_mat <- merged@layers$counts[, colnames(merged) %in% barcode_data$X]
keep_cell_barcodes = colnames(merged)[colnames(merged) %in% barcode_data$X]
colnames(wanted_mat) = keep_cell_barcodes
rownames(wanted_mat) = rownames(merged)
new_seurat_object <- CreateSeuratObject(wanted_mat, project = 'GSE221575' )
new_seurat_object = AddMetaData(new_seurat_object, metadata = barcode_data )
saveRDS(new_seurat_object, 'GSE221575_gd.rds')

