run <- function(RDR_var, K, seed, burnin_tol, Gibbs_tol){
  
  set.seed(seed)
  ret <- matrix(NA, 4, 2)
  colnames(ret) <- c("Clusters", "States")
  rownames(ret) <- c("RDR", "both", "hclust_Ktrue", "hclust_Khat")
  ret[4, 2] <- seed
  
  K_hat = round(K*2.5)       ## Input cluster number, same for Chloris and hclust
  N = 300
  U = 200
  
  sims <- get_sim_data(K = K, N = N, U = U, expr = F, RDR_var = RDR_var, Q_to_neutral = 1, BAF_missing_percent = 0.95)
  
  #### Chloris
  set.seed(seed)
  res_RDR <- Chloris(sims$RDR, burnin_tol = burnin_tol, Gibbs_tol = Gibbs_tol, cluster_shrink_tol = 20, min_cluster_size = 1, K = K_hat, verbose = F)
  
  set.seed(seed)
  res_both <- Chloris(sims$RDR, A = sims$A, D = sims$D, burnin_tol = burnin_tol, Gibbs_tol = Gibbs_tol, 
                      cluster_shrink_tol = 20, min_cluster_size = 1, K = K_hat, verbose = F)
  
  #### hclust
  hc <- hclust(dist(t(sims$RDR)), method = 'ward.D2')
  hc_Khat <- cutree(hc, k = K_hat)
  hc_Ktrue <- cutree(hc, k = K)
  
  
  ret[1, 2] <- state_align(res_RDR$state_est, sims$states_true, res_RDR$par_record)
  ret[2, 2] <- state_align(res_both$state_est, sims$states_true, res_both$par_record)

  ret[1, 1] <- aricode::ARI(res_RDR$cluster_est, sims$cluster_true)*100
  ret[2, 1] <- aricode::ARI(res_both$cluster_est, sims$cluster_true)*100
  ret[3, 1] <- aricode::ARI(hc_Ktrue, sims$cluster_true)*100
  ret[4, 1] <- aricode::ARI(hc_Khat, sims$cluster_true)*100
  
  
  cat("\n")
  print(ret)
  ret
  
}