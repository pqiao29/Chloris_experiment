##### get gene_record
gene_record_raw <- readxl::read_excel("Melanoma_genes.xlsx", col_names = FALSE)
colnames(gene_record_raw)[1:3] <- c("gene", "CNA", "ref")
gene_record_raw <- gene_record_raw[gene_record_raw$gene %in% rowData(sce)$gene_name, ]

################################################ Curate references ################################################
## Obtain the number of references
refs <- unique(unlist(lapply(gene_record_raw[, 3][[1]], function(x) strsplit(x, ","))))
refs <- sort(unique(as.numeric(refs)))
## Expand gene_record s.t. each reference takes a column
gene_record <- matrix(0, nrow(gene_record_raw), (2 + max(refs)))
gene_record[, 1] <- gene_record_raw[, 1][[1]]
gene_record[, 2] <- gene_record_raw[, 2][[1]]
colnames(gene_record) <- c("gene", "CNA", paste0("ref", 1:max(refs)))
## Add indicator of gene--reference 
for(g in 1:nrow(gene_record)){
    g_ref <- gene_record_raw[g, 3][[1]]
    g_ref <- strsplit(g_ref, ",")[[1]]
    g_ref <- as.numeric(g_ref)
    gene_record[g, g_ref + 2] <- 1
}
## sort references from more->less genes mentioned
tmp_refs <- gene_record[, -(1:2)]
tmp_refs <- matrix(as.numeric(tmp_refs), nrow(gene_record))

ref_order <- sort(colSums(tmp_refs), index.return = T, decreasing = T)$ix
gene_record <- gene_record[, c(1, 2, ref_order + 2)]

keep_refs <- colSums(tmp_refs)[ref_order] > 0 # delete empty references
gene_record <- gene_record[, c(T, T, keep_refs)]
cat(sum(keep_refs), "reference papers.\n")

selected_genes <- gene_record[rowSums(tmp_refs) >= 2, "gene"]
cat(length(selected_genes), "genes appear in more than one paper.\n")
######################################################################################################################


######################################################### plot #########################################################
#### gene names --> index ------------------------------------------------------------------------------
gene_idx <- NULL
for(g in selected_genes) gene_idx <- c(gene_idx, which(rowData(sce)$gene_name == g))
stopifnot(length(unique(gene_idx)) == length(gene_idx))

gene_order <- sort(gene_idx, decreasing = T, index.return = T)$ix ## reverse gene order so chr1 is plotted on top
gene_idx <- gene_idx[gene_order]
selected_genes <- selected_genes[gene_order]
stopifnot(identical(selected_genes, rowData(sce[gene_idx, ])$gene_name))
rownames(gene_record) <- gene_record[, 1]
gene_record <- gene_record[selected_genes, ]

#### barplot: ref count ------------------------------------------------------------------------------
d <- data.frame(x = nrow(gene_record):1, y = rowSums(matrix(as.numeric(gene_record[, -(1:2)]), nrow(gene_record))))
# interpolate values from zero to y and create corresponding number of x values
vals <- lapply(d$y, function(y) seq(0, y, by = 0.01))
y <- unlist(vals)
mid <- rep(d$x, lengths(vals))
d2 <- data.frame(x = mid - 0.4,
                 xend = mid + 0.4,
                 y = y,
                 yend = y)

gg_refcount <- ggplot(data = d2, aes(x = x, xend = xend, y = y, yend = yend, color = y)) +
    geom_segment(linewidth = 2) +
    scale_color_gradient2(low = "#FDD262", high = "#5B1A18", mid = "#F98400", midpoint = 7) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", 
          panel.background = element_blank(),
          axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.title.y = element_blank(), 
          axis.text.y = element_blank()) 
ggsave(gg_refcount, file = "Melanoma_gene_count.pdf", width = 20, height = 5)

