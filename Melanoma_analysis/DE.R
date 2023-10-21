load("M71.Rdata")
res <- readRDS("Chloris_fullmode.rds")

K <- length(unique(res$cluster_est))
## Adjust cluster order:
cluster_order <- order(apply(res$state_est, 1, function(x) sum(x != 2)))
cluster_est <- NULL
for(k in 1:K) cluster_est[res$cluster_est == cluster_order[k]] <- k

#### filtering 
sce$Chloris_clustering <- paste0("clone", cluster_est)
gene_use_DE <- rownames(sce)[rowMeans(counts(sce)) > 0.5] ## remove low expressed genes
sce_filt <- sce[gene_use_DE,]
rm(sce)
############################################### QL F-test #################################################
dge_list <- edgeR::DGEList(round(counts(sce_filt)))
dge_list <- edgeR::calcNormFactors(dge_list, method = "TMM")
sce_filt$cdr <- colSums(counts(sce_filt) > 0) / nrow(sce_filt)
design_list <- model.matrix(~cdr + Chloris_clustering, data = colData(sce_filt))
dge_list <- edgeR::estimateDisp(dge_list, design_list)
fit_list <- edgeR::glmQLFit(dge_list, design_list)
qlf_list <- edgeR::glmQLFTest(fit_list,
                              coef = (ncol(design_list) - K + 2):ncol(design_list))
sum(p.adjust(qlf_list$table$PValue, method = "BH") <= 0.05, na.rm = TRUE)
print(summary(edgeR::decideTestsDGE(qlf_list)))

### pairwise clones ========================================
cat("....calculating DE \n")
num_covars <- ncol(design_list) - (K - 1)
num_contrs <- choose(K, 2)
out_list <- list()
out_list$base_clone <- levels(factor(sce_filt$Chloris_clustering))[1]
for (j in 2:K) {
    coef_idx <- (ncol(design_list) - K + j)
    tmp <- edgeR::glmQLFTest(fit_list, coef = coef_idx)
    comp <- paste0(
        gsub("assigned", "", colnames(design_list)[coef_idx]),
        "_",
        out_list$base_clone)
    out_list[[comp]] <- tmp
}

if (K > 2) {
    contrasts <- list()
    n <- 1
    for (k in seq_len(K - 1)) {
        for (m in seq_len(K - 1)) {
            if (m > k) {
                contr <- rep(0, ncol(design_list))
                contr[num_covars + m] <- -1
                contr[num_covars + k] <- 1
                contrasts[[n]] <- contr
                n <- n + 1
            }
        }
    }
    for (p in seq_len(length(contrasts))) {
        tmp <- edgeR::glmQLFTest(fit_list, contrast = contrasts[[p]])
        tmp2 <- strsplit(tmp$comparison, split = " ")[[1]]
        comp <- paste0(
            gsub("1\\*assigned", "", tmp2[2]),
            "_",
            gsub("-1\\*assigned", "", tmp2[1]))
        out_list[[comp]] <- tmp
    }
}
qlf_pairwise_list <- out_list


############################################### Gene set testing #################################################
library(limma)
library(org.Hs.eg.db) 
xx <- as.list(org.Hs.egENSEMBL2EG)

get_camera <- function(qlf_pairwise_list, Hs){
    ret <- list()
    for (j in 2:length(qlf_pairwise_list)) {
        idx <- ids2indices(Hs, identifiers = qlf_pairwise_list[[j]]$table$entrezid)
        ret[[names(qlf_pairwise_list)[j]]] <- list()
        ret[[names(qlf_pairwise_list)[j]]][["logFC"]] <-
            cameraPR(statistic = qlf_pairwise_list[[j]]$table$logFC, idx,
                     inter.gene.cor = 0.01)
        ret[[names(qlf_pairwise_list)[j]]][["signF"]] <-
            cameraPR(statistic = (sign(qlf_pairwise_list[[j]]$table$logFC) *
                                      qlf_pairwise_list[[j]]$table$F), idx,
                     inter.gene.cor = 0.01)
    }
    ret
}

for (j in 2:length(qlf_pairwise_list)) {
    qlf_pairwise_list[[j]]$table$ensembl_gene_id <-
        strsplit2(rownames(qlf_pairwise_list[[j]]$table), split = "_")[,1]
    mm <- match(qlf_pairwise_list[[j]]$table$ensembl_gene_id,
                rowData(sce_filt)$gene_id)
    qlf_pairwise_list[[j]]$table$hgnc_symbol <-
        rowData(sce_filt)$symbol[mm]
    qlf_pairwise_list[[j]]$table$entrezid <- NA
    for (k in seq_len(nrow(qlf_pairwise_list[[j]]$table))) {
        if (qlf_pairwise_list[[j]]$table$ensembl_gene_id[k] %in% names(xx))
            qlf_pairwise_list[[j]]$table$entrezid[k] <- xx[[qlf_pairwise_list[[j]]$table$ensembl_gene_id[k]]][1]
    }
}

## Hallmark gene sets (H) 
load("human_H_v5p2.rdata")  # Hs.H
camera_msigdb_H_list <- get_camera(qlf_pairwise_list, Hs.H)

############################################### Save results #################################################
de_results_list <- list()
de_results_list[["camera"]] <- list()
de_results_list[["camera"]][["H"]] <- camera_msigdb_H_list
de_results_list[["design_list"]] <- design_list
de_results_list[["dge_list"]] <- dge_list
de_results_list[["fit_list"]] <- fit_list
de_results_list[["qlf_list"]] <- qlf_list
de_results_list[["qlf_pairwise"]] <- qlf_pairwise_list
de_results_list[["sce_filt"]] <- sce_filt

saveRDS(de_results_list, file = "de_results.rds")
