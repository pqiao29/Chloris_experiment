compare <- function(seed){
  
  set.seed(seed)
  ret <- matrix(NA, 7, 4)
  colnames(ret) <- c("Clusters", "States", "K", "time")
  rownames(ret) <- c("Chloris RDR", "Chloris both(10%)",
                     "inferCNV", "HB RDR", "HB both(10%)", "HB both(50%)", "HB both(100%)")
  ret[7, 3] <- seed
  
  K = 5
  K_hat = 12
  N = 300
  U = 200
  
  #### Chloris: RDR + BAF
  for(m in 3:1){
    
    set.seed(seed)
    sims <- get_sim_data(K = K, N = N, U = U, expr = T, BAF_missing_percent = c(0.9, 0.5, 0)[m], Q_to_neutral = 1)
    ref_idx <- sims$cluster_true == 1
    RDR <- preprocess(sims$expr, ref_idx)
  
    #### HB: both 
    start <- Sys.time()
    hb_res <- run_HoneyBADGER(sims$expr, sims$cluster_true == 1, sims$A, sims$D, K)
    end <- Sys.time()
    ret[m + 4, 1] <- aricode::ARI(hb_res$clusters_both, sims$cluster_true)*100
    ret[m + 4, 4] <- end - start
    states_true_for_HB <- sims$states_true
    states_true_for_HB[sims$states_true == 4] <- 3
    ret[m + 4, 2] <- state_align(hb_res$states_both, states_true_for_HB)
  }
  
  #### HB: RDR 
  start <- Sys.time()
  hb_res <- run_HoneyBADGER(sims$expr, sims$cluster_true == 1, NULL, NULL, K)
  end <- Sys.time()
  ret[4, 1] <- aricode::ARI(hb_res$clusters_RDR, sims$cluster_true)*100
  ret[4, 2] <- state_align(hb_res$states_RDR, states_true_for_HB)
  ret[4, 4] <- end - start 
  
  #### Chloris: RDR 
  set.seed(seed)
  start <- Sys.time()
  res_RDR <- Chloris(RDR, burnin_tol = 200, Gibbs_tol = 200, cluster_shrink_tol = 20, min_cluster_size = 1, K = K_hat, verbose = F)
  end <- Sys.time()
  ret[1, 1] <- aricode::ARI(res_RDR$cluster_est, sims$cluster_true)*100
  ret[1, 2] <- state_align(res_RDR$state_est, sims$states_true)
  ret[1, 3] <- length(unique(res_RDR$cluster_est))
  ret[1, 4] <- end - start 
  
  set.seed(seed)
  #### Chloris: both 
  start <- Sys.time()
  res_both <- Chloris(RDR, A = sims$A, D = sims$D, burnin_tol = 200, Gibbs_tol = 200, 
                      cluster_shrink_tol = 20, min_cluster_size = 1, K = K_hat, verbose = F)
  end <- Sys.time()
  ret[2, 1] <- aricode::ARI(res_both$cluster_est, sims$cluster_true)*100
  ret[2, 2] <- state_align(res_both$state_est, sims$states_true)
  ret[2, 3] <- length(unique(res_both$cluster_est))
  ret[2, 4] <- end - start
  
  #### inferCNV
  start <- Sys.time()
  state_emission_params <- infercnv_emission_pars(sims$expr, sims$expr[, ref_idx], 51)
  idx_mu_neutral <- which.min(abs(state_emission_params$mean))
  res_infercnv <- run_infercnv(RDR, state_emission_params$mean[idx_mu_neutral + -1:2], state_emission_params$sd[idx_mu_neutral + -1:2], idx_mu_neutral)
  end <- Sys.time()
  ret[3, 1] <- aricode::ARI(res_infercnv$clusters, sims$cluster_true)*100
  ret[3, 2] <- state_align(t(res_infercnv$states), sims$states_true)
  ret[3, 3] <- length(unique(res_infercnv$clusters))
  ret[3, 4] <- end - start
  
  print(ret)
  return(ret)
}