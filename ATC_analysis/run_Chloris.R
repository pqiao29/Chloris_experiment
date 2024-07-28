library(SingleCellExperiment)
sce <- readRDS("ATC.rds")
sce

#### (optional) exclude chr X
chrs <- rowData(sce)$gene_order$chromosome_name
chrs[chrs == "X"] <- 23
chrs <- as.numeric(chrs)
break_idx <- which(chrs[-1]  - chrs[-length(chrs)] > 0) + 1

#### run model
library(Chloris)
start <- Sys.time()
res <- Chloris(RDR = assays(sce)$RDR, A = NULL, D = NULL, min_cluster_size = 2, 
               S = 4, CN_to_neutral = 100, prior_Q_diag = 500, break_idx = break_idx, 
               burnin_tol = 500, Gibbs_tol = 500, verbose = F)
end <- Sys.time()
print(end - start)

saveRDS(res, file = "Chloris_RDRmode.rds")

#### plot result
gg <- plot_inout(assays(sce)[["RDR"]], cluster_labels = list(res$cluster_est), 
                 CN_states = list(res$state_est))
ggsave(gg, width = 20, height = 10, file = "Chloris_RDRmode.png")

### clonal tree
DD <- dist(res$state_est)
tre <- ape::nj(DD)
tre <- ape::ladderize(tre)
plot(tre, cex = 0.6)
