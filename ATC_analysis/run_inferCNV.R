library(SingleCellExperiment)
sce <- readRDS("ATC_RDR.rds")
sce

cell_annot <- matrix("tumor", ncol(sce), 1)
rownames(cell_annot) <- colnames(sce)
cell_annot[colData(sce)$is_norm, ] <- "normal"

infercnv_obj <- infercnv::CreateInfercnvObject(raw_counts_matrix = assays(sce)[[counts]],
                                               gene_order_file = rowData(sce)$gene_order,
                                               annotations_file = cell_annot,
                                               ref_group_names = c("normal"))

system.time({
  infercnv_res <- infercnv::run(infercnv_obj,
                                cutoff = .1, 
                                cluster_by_groups = TRUE, 
                                denoise = TRUE,
                                HMM = TRUE)
})

saveRDS(infercnv_res, "inferCNV.rds")
