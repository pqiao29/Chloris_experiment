import::from(dplyr, `%>%`)
library(ggplot2)

de_res <- readRDS(file = "de_results.rds")
camera_msigdb_H_list <- de_res[["camera"]][["H"]]
sce_filt <- de_res[["sce_filt"]]

## Extract significant Hallmark gene sets
idx_base_clone <- grep("*_clone1", names(camera_msigdb_H_list))
sig_gene_sets_ls <- list()
for(i in idx_base_clone){
    logFC <- camera_msigdb_H_list[[i]]$logFC
    logFC_filt <- logFC[logFC$FDR <= 0.05, ]
    sig_gene_sets_ls[[i]] <- rownames(logFC_filt)
}
sig_gene_sets <- unique(unlist(sig_gene_sets_ls))
y <- gsub('HALLMARK_*', '', sig_gene_sets)

### Extract direction for each gene set
logP <- matrix(NA, length(y), length(idx_base_clone))
direction <- logP
is_Sig <- logP
for(i in idx_base_clone){
    logP[, i] <- camera_msigdb_H_list[[i]]$logFC[sig_gene_sets, ]$PValue
    direction[, i] <- camera_msigdb_H_list[[i]]$logFC[sig_gene_sets, ]$Direction
    is_Sig[, i] <- sig_gene_sets %in% sig_gene_sets_ls[[i]]
}

### From gene sets to gene ID
load("human_H_v5p2.rdata") 
library(org.Hs.eg.db) 
xx <- as.list(org.Hs.egENSEMBL2EG)

genes_in_data <- rownames(sce_filt)
idx <- Hs.H[sig_gene_sets]
idx_ensembl_gene_id <- list()
for(gene_set in names(idx)){
    tmp_gene_id <- NULL
    for(gene in 1:length(idx[[gene_set]])){
        tmp_gene_id <- c(tmp_gene_id, xx[xx == idx[[gene_set]][gene]] %>% names)
    }
    idx_ensembl_gene_id[[gene_set]] <- intersect(tmp_gene_id, genes_in_data)
    cat(length(idx_ensembl_gene_id[[gene_set]]), " out of ", length(tmp_gene_id), " genes in gene set ", gene_set, "\n")
}

### DE_pval for plot
K <- length(unique(sce_filt$Chloris_clustering))
DE_dir <- matrix(NA, K - 1, nrow(sce_filt))
rownames(direction) <- sig_gene_sets
for(gene_set in sig_gene_sets){
    id_genes <- idx_ensembl_gene_id[gene_set][[1]]
    idx_genes <- which(rownames(sce_filt) %in% id_genes)
    DE_dir[, idx_genes] <- direction[gene_set, ]
}
print(table(DE_dir))

#### plot
sig_num <- nrow(direction)
x <- paste0("clone", 2:K, "--ref")
plot_df <- expand.grid(X = x, Y = rev(y))
plot_df$Dir <- c(t(direction[sig_num:1, ]))
plot_df$is_Sig <- c(t(is_Sig[sig_num:1, ]))
plot_df$logP <- c(t(log10(logP)[sig_num:1, ]))*(-1)
plot_df$blured_logP = cut(plot_df$logP, breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 25))

gg_base_clone <- ggplot(plot_df, aes(X, Y, fill= factor(Dir))) +
    geom_tile(aes(alpha = plot_df$blured_logP)) +
    geom_point(data = plot_df[plot_df$is_Sig, ]) +
    scale_fill_manual(values = c("#2686A0", "#C7522B")) +
    #scale_y_discrete(position = "right") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none",
          panel.background = element_blank(),
          axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          axis.title.y = element_blank())
ggsave(gg_base_clone, file = "Hallmark_toclone1.png", height = 6, width = 4)

save(direction, idx, idx_ensembl_gene_id, file = "de_index.rds")
