
load("M71.Rdata")

##### some model parameters
chrs <- rowData(sce)$gene_order$chr
break_idx <- which(chrs[-1]  - chrs[-length(chrs)] > 0) + 1

#### run model
res <- Chloris(RDR = RDR, A = A, D = D, break_idx = break_idx,
               S = 4, prior_mu = prior_mu, prior_Q_diag = 500, 
               burnin_tol = 500, Gibbs_tol = 500, min_cluster_size = 2)

saveRDS(res,  file = "Chloris_fullmode.rds")

#### plot result
U <- nrow(RDR)
cluster_order <- order(apply(res$state_est, 1, function(x) sum(x != 2)))
cluster_est <- NULL
for(k in 1:K) cluster_est[res$cluster_est == cluster_order[k]] <- k
res$cluster_est <- cluster_est
res$state_est <- res$state_est[cluster_order, ]

gg <- plot_inout(RDR[U:1, ], "RDR", list(res$cluster_est), res$state_est[, U:1], lim = c(-1, 1))
ggsave(gg, file = "Chloris_fullmode.png", width = 200, height = 600, unit = "mm")
