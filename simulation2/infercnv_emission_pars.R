import::here("infercnv_utils.R", sim_meanvar, smooth_helper)

infercnv_emission_pars <- function(obs, ref, window_length){
  
  if(nrow(obs) != nrow(ref)) stop("obs and ref group need to have the same number of genes!")
  U <- nrow(obs)
  
  num_genes_per_chr = floor(U/6)
  num_remaining = U - 5*num_genes_per_chr
  cnv_info = list(list(name='chrA', cnv = 1, ngenes = num_genes_per_chr),
                  list(name='chrB', cnv = 0.01, ngenes = num_genes_per_chr),
                  list(name='chrC', cnv = 0.5, ngenes = num_genes_per_chr),
                  list(name='chrD', cnv = 1.5, ngenes = num_genes_per_chr),
                  list(name='chrE', cnv = 2, ngenes = num_genes_per_chr),
                  list(name='chrF', cnv = 3, ngenes = num_remaining))
  cnv_idx <- c(rep(c(1, 0.01, 0.5, 1.5, 2), each = num_genes_per_chr), rep(3, num_remaining))
  
  hspike <- infercnv_add_spike(obs, ref, cnv_info, cnv_idx)
  hspike <- infercnv_substract_ref(hspike$obs, hspike$ref)
  hspike$obs <- infercnv_smoothing_recenter(window_length, hspike$obs)
  hspike$ref <- infercnv_smoothing_recenter(window_length, hspike$ref)
  
  return(infercnv_est_emission_pars(hspike$obs, hspike$ref, cnv_info, cnv_idx))
}

infercnv_add_spike <- function(obs, ref, cnv_info, cnv_idx){
    
    num_cells = 100
    num_total_genes = nrow(obs)
    gene_names <- rownames(obs)
    
    gene_means_orig = rowMeans(ref)
    genes_means_use_idx = sample(x = seq_len(num_total_genes), size = num_total_genes, replace = TRUE) 
    gene_means = gene_means_orig[genes_means_use_idx]
    gene_means[gene_means==0] = 1e-3 # just make small nonzero values
    names(gene_means) = gene_names       
    
    ## simulate normal cells via mean-var trend
    if(ncol(ref) > 1){
        sim_normal_matrix <- sim_meanvar(ref_gexp = ref, gene_means, num_cells)
    }else{
        sim_normal_matrix <- sim_meanvar(ref_gexp = obs, gene_means, num_cells)
    }
    
    colnames(sim_normal_matrix) = paste0("simnorm_cell_", seq_len(num_cells))
    rownames(sim_normal_matrix) = gene_names
    ## simulate cells with CNV via mean-var trend
    hspike_gene_means = gene_means
    for (info in cnv_info) {
        cnv = info$cnv
        tmp_cnv_idx = which(cnv_idx == cnv)
        hspike_gene_means[tmp_cnv_idx] =  hspike_gene_means[tmp_cnv_idx] * cnv
    }
    
    if(ncol(ref) > 1){
        sim_spiked_cnv_matrix <- sim_meanvar(ref_gexp = ref, hspike_gene_means, num_cells) 
    }else{
        sim_spiked_cnv_matrix <- sim_meanvar(ref_gexp = obs, hspike_gene_means, num_cells) 
    }
    
    colnames(sim_spiked_cnv_matrix) = paste0("spike_tumor_cell_", seq_len(num_cells))
    rownames(sim_spiked_cnv_matrix) = gene_names
    
    list("obs" = sim_spiked_cnv_matrix, "ref" = sim_normal_matrix)
}

infercnv_substract_ref <- function(obs, ref){
    
    log_obs <- log2(obs + 1)
    log_ref <- log2(ref + 1)
    Ave_ref <- rowMeans(log_ref)
    log_obs <- log_obs - Ave_ref
    log_ref <- log_ref - Ave_ref
    
    list("obs" = 2^log_obs, "ref" = 2^log_ref)
}

infercnv_smoothing_recenter <- function(window_length, expr){
    smoothed_data <- apply(expr, 2, smooth_helper, window_length = window_length)
    row_median <- apply(smoothed_data, 2, function(x) { median(x, na.rm=TRUE) } )
    expr_data <- t(apply(smoothed_data, 1, "-", row_median))
    return(expr_data)
}

infercnv_est_emission_pars <- function(obs, ref, cnv_info, cnv_idx){
    
    num_cells = 100
    gene_expr_by_cnv = list()
    for (info in cnv_info) {
        cnv <- info$cnv
        cnv_name = sprintf("cnv:%g", cnv)
        gene_expr_by_cnv[[cnv_name]] = c(obs[cnv_idx == cnv, ])
    }
    
    cnv_mean_sd = list()
    for (cnv_level in names(gene_expr_by_cnv) ) {
        gene_expr = gene_expr_by_cnv[[ cnv_level ]]
        gene_expr_mean = mean(gene_expr)
        gene_expr_sd = sd(gene_expr)
        cnv_mean_sd[[ cnv_level ]] = list(mean = gene_expr_mean, sd = gene_expr_sd)
    } 
    
    cnv_level_to_mean_sd = list()
    for (cnv_level in names(gene_expr_by_cnv) ) {
        expr_vals = gene_expr_by_cnv[[ cnv_level ]]
        nrounds = 100
        sds <- vapply(seq_len(100), function(ncells) {
            vals <- replicate(nrounds, sample(expr_vals, size=ncells, replace=TRUE)) ## What's the difference? sample(expr_vals, size=ncells*nrounds, replace=TRUE)
            if ("matrix" %in% is(vals)) { means <- Matrix::rowMeans(vals) } else { means <- mean(vals) }
            sd(means) }, numeric(1))
        cnv_level_to_mean_sd[[ cnv_level ]] <- sds
    }
    
    tmp_names <- names(cnv_level_to_mean_sd)
    cnv_level_to_mean_sd_fit <- lapply(tmp_names, function(cnv_level) {
        sd_vals = cnv_level_to_mean_sd[[ cnv_level ]]
        num_cells = seq_along(sd_vals)
        lm(log(sd_vals) ~ log(num_cells))
    })
    names(cnv_level_to_mean_sd_fit) <- tmp_names    
    
    for (cnv_level in names(cnv_mean_sd)) {
        fit <- cnv_level_to_mean_sd_fit[[cnv_level]]
        sd <- exp(predict(fit, newdata = data.frame(num_cells=num_cells))[[1]])
        cnv_mean_sd[[cnv_level]]$sd <- sd
    }
    
    list(mean=c(cnv_mean_sd[["cnv:0.01"]]$mean, cnv_mean_sd[["cnv:0.5"]]$mean,
                cnv_mean_sd[["cnv:1"]]$mean, cnv_mean_sd[["cnv:1.5"]]$mean,
                cnv_mean_sd[["cnv:2"]]$mean, cnv_mean_sd[["cnv:3"]]$mean),
         sd=c(cnv_mean_sd[["cnv:0.01"]]$sd, cnv_mean_sd[["cnv:0.5"]]$sd,
              cnv_mean_sd[["cnv:1"]]$sd, cnv_mean_sd[["cnv:1.5"]]$sd,
              cnv_mean_sd[["cnv:2"]]$sd, cnv_mean_sd[["cnv:3"]]$sd) )
    
}
