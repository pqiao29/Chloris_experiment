# import::here("infercnv_HMM.R", infercnv_HMM)
# import::here("infercnv_clustering_randomtree.R", single_tumor_subclustering_recursive_random_smoothed_trees)


run_infercnv <- function(obs, mu, sd, neutral_idx){
    # Clustering =================================================================================================================
    p_val = 0.1
    window_size = 101
    max_recursion_depth = 3
    min_cluster_size_recurse = 10
    colnames(obs) <- paste0("cell", 1:ncol(obs))
    
    ## smooth and median-center
    sm_obs = apply(obs, 2, caTools::runmean, k = window_size)
    row_median <- apply(sm_obs, 2, function(x) { median(x, na.rm=TRUE) } )
    sm_obs <- t(apply(sm_obs, 1, "-", row_median))
    
    tumor_subcluster_info = list()    
    ## Hierachical clustering 
    hc <- hclust(dist(t(sm_obs)), method = 'ward.D2')
    tumor_subcluster_info$hc = hc
    heights = hc$height
    
    grps <- rep(sprintf("cluster.%d", 1), ncol(obs))
    names(grps) <- colnames(obs)
    
    grps <- single_tumor_subclustering_recursive_random_smoothed_trees(obs, hclust_method = 'ward.D2', p_val, grps, window_size,
                                                                       max_recursion_depth, min_cluster_size_recurse)
    
    
    infercnv_clust <- as.integer(factor(grps, levels = unique(grps), labels = 1:length(unique(grps))))
    names(infercnv_clust) <- NULL
    
    # HMM =======================================================================================================================
    K_est <- length(unique(infercnv_clust))
    
    infercnv_clust <- factor(infercnv_clust, levels = unique(infercnv_clust), labels = 1:length(unique(infercnv_clust)))
    infercnv_clust <- as.integer(infercnv_clust)
    
    states <- infercnv_HMM(infercnv_clust, mu, sd, neutral_idx, obs = obs, K = K_est)
    
    return(list("clusters" = infercnv_clust, "states" = states))
    
}