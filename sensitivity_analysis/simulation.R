run <- function(min_cluster_size, cluster_shrink_tol, seed){
  
  N = 300
  U = 200
  K = 5
  
  ret <- matrix(NA, 2, 2)
  rownames(ret) <- c("RDR", "both")
  colnames(ret) <- c("K_hat", "SysTime")
  
  sims <- get_sim_data(K = K, N = N, U = U, expr = T, BAF_missing_percent = 0.9, Q_to_neutral = 1)
  ref_idx <- sims$cluster_true == 1
  RDR <- preprocess(sims$expr, ref_idx)
  
  K_hat = max(min(round(N/min_cluster_size) + 2, 12), K)
  
  #### Chloris: RDR 
  start <- Sys.time()
  res_RDR <- Chloris(RDR, burnin_tol = 200, Gibbs_tol = 200, cluster_shrink_tol = cluster_shrink_tol, min_cluster_size = min_cluster_size, K = K_hat, verbose = F)
  end <- Sys.time()
  
  ret[1, 1] <- length(unique(res_RDR$cluster_est))
  ret[1, 2] <- end - start
  
  #### Chloris: both 
  start <- Sys.time()
  res_both <- Chloris(RDR, A = sims$A, D = sims$D, burnin_tol = 200, Gibbs_tol = 200, 
                      cluster_shrink_tol = cluster_shrink_tol, min_cluster_size = min_cluster_size, K = K_hat, verbose = F)
  end <- Sys.time()

  ret[2, 1] <- length(unique(res_both$cluster_est))
  ret[2, 2] <- end - start
  
  print(ret)
  return(ret)
}