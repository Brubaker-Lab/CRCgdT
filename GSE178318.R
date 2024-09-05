library(Seurat)
data <- Read10X(data.dir = 'F:/CRC/AA_DONE/GSE178318/')
seurat_object <- CreateSeuratObject(counts = data)
saveRDS(seurat_object,'F:/CRC/GSE178318/srat.rds')
# Data integrated, calling it integrated data from now on
integrated_data = readRDS('F:/CRC/AA_Done/GSE178318/srat.rds')
#original author encode the condition in the cell name, so we need to retrieve it
obs_names = colnames(integrated_data)
barcodes = list()
patient =list()
tissue = list()
# Initialize vectors
barcodes <- vector("character", length(obs_names))
patient <- vector("character", length(obs_names))
tissue <- vector("character", length(obs_names))

for (i in seq_along(obs_names)) {
  parts <- strsplit(obs_names[i], split='_', fixed=TRUE)[[1]]
  if (length(parts) >= 3) {
    barcodes[i] <- parts[1]
    patient[i] <- parts[2]
    tissue[i] <- parts[3]
  } else {
    warning(sprintf("obs_names[%d] does not contain three parts: %s", i, obs_names[i]))
  }
}

# Assign to integrated_data's metadata
integrated_data@meta.data[['barcodes']] <- barcodes
integrated_data@meta.data[['patient']] <- patient
integrated_data@meta.data[['tissue']] <- tissue

integrated_data <- subset(integrated_data, subset = patient %in% list("COL07", "COL12","COL16"))
integrated_data <- subset(integrated_data, subset = tissue %in% list("CRC"))

#Start preprocessing
mito.genes <- grep("^MT-", rownames(integrated_data), value = TRUE)
# Calculate the percentage of mitochondrial gene expression per cell
integrated_data[["percent.mt"]] <- PercentageFeatureSet(integrated_data, pattern = "^MT-")
integrated_data <- subset(integrated_data, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA > 500 & percent.mt < 20)
saveRDS(integrated_data,'F:/CRC/AA_Done/GSE178318/processed_GSE178318.rds')
