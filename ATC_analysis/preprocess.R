expr = readRDS(url('http://pklab.med.harvard.edu/teng/data/count_mat_ATC2.rds'))
colnames(expr) <- gsub("-1$", "", gsub("^ATC2_", "", colnames(expr)))

################################# gene order ################################# 
mart.obj <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = 'hsapiens_gene_ensembl', host = "jul2015.archive.ensembl.org")
gene_order <- biomaRt::getBM(values = rownames(expr),
                             attributes = c("hgnc_symbol", "chromosome_name","start_position","end_position"),
                             filters = c("hgnc_symbol"), mart = mart.obj)
gene_order <- gene_order[gene_order$chromosome_name %in% c(1:22, "X"), ]
gene_order <- gene_order[order(gene_order$chromosome_name, gene_order$start_position, gene_order$end_position), ]
print(length(intersect(gene_order$hgnc_symbol, rownames(expr))))
expr <- expr[gene_order$hgnc_symbol, ]

library(SingleCellExperiment)
sce <- SingleCellExperiment(assays=list(counts = expr))
rowData(sce)$gene_order <- gene_order

#################################  Gene QC ################################# 
keep_genes <- rowSums(counts(sce) > 0) >= 20 ## expr >= 1 in >= 5 cells
sce <- sce[keep_genes, ]

#################################  Get RDR ################################# 
res_BAF <- readRDS("Chloris_BAFref.rds")
imbal_cnt <- apply(res_BAF$state_est, 1, function(x) sum(x != 2))
norm_cluster <- which.min(imbal_cnt)
colData(sce)$is_norm <- res_BAF$cluster_est == norm_cluster

## Normalization
tmp_data <- counts(sce)
cs = colSums(tmp_data)
tmp_data <- sweep(tmp_data, STATS = cs, MARGIN = 2, FUN="/")  
normalize_factor = median(cs)  # normalize_factor
normalized <- tmp_data * normalize_factor

## Get RDR from gene expr
normalized <- log2(normalized + 1)
ref <- rowMeans(normalized[, norm_BAF])
RDR <- normalized - ref

## Window smoothing
## smooth_helper is an internal function in infercnv: https://github.com/broadinstitute/infercnv/
RDR <- apply(RDR, 2, smooth_helper, window_length = 101)

## Recentering
col_median <- apply(RDR, 2, function(x) { median(x, na.rm=TRUE) })
RDR <- t(apply(RDR, 1, "-", col_median))

assays(sce)$RDR <- RDR
saveRDS(sce, "ATC_RDR.rds")
