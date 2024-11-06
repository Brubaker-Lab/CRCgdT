## First, back in python, save the SCT corrected count of ONLY the cell type you want to run the DEG with the 'tissue' specification


srat = readRDS("your_isoated_seurat_object_with_2_tissue_types")
cds <- as.cell_data_set(srat)
cds <- estimate_size_factors(cds)
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(srat[["RNA"]])

gene_fits <- fit_models(cds, model_formula_str = "~tissue", cores = 6)
fit_coefs <- coefficient_table(gene_fits)

# Adjust condition in the filter
terms <- fit_coefs %>% filter(term == "tissueCarcinoma")
terms = terms %>% select(gene_short_name, term, q_value, normalized_effect)