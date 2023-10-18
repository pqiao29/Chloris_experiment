run_HoneyBADGER <- function(expr, ref_idx, A = NULL, D = NULL, K){
  
  ret <- list()
  U <- nrow(expr)
  ### =============================================================== RDR ===============================================================
  #hb_RDR <- my_setGexpMats(expr, ref)$gexp.norm   ## Does not perform well
  hb_RDR <- preprocess(expr, ref_idx)
  mvFit <- HoneyBADGER::setMvFit(hb_RDR, num.genes = seq(5, 100, by = 10), verbose = FALSE)
  dev <- HoneyBADGER::setGexpDev(hb_RDR, alpha = 0.05, verbose = FALSE)
  
  hb_potentialCnvs <- my_calcGexpCnvBoundaries(hb_RDR, m = dev)
  potentialCnvs_RDR <- list()
  potentialCnvs_RDR <- c(potentialCnvs_RDR, lapply(hb_potentialCnvs$amp, as.numeric))
  potentialCnvs_RDR <- c(potentialCnvs_RDR, lapply(hb_potentialCnvs$del, as.numeric))
  
  probs_RDR <- do.call(rbind, lapply(potentialCnvs_RDR, function(x){
    my_calcGexpCnvProb(gexp.norm = hb_RDR, mvFit = mvFit, m = dev, region = as.numeric(x))
  }))
  probs_RDR <- do.call(rbind, probs_RDR)
  probs_for_clustering_RDR <- probs_RDR[apply(probs_RDR, 1, max) > 0.5,  , drop = FALSE]
  
  ## cluster cells on posterior probability
  hc_RDR <- hclust(dist(t(probs_for_clustering_RDR)), method='ward.D')
  hb_cluster_RDR <- cutree(hc_RDR, K)
  
  ret$"clusters_RDR" <- hb_cluster_RDR
  
  ## States on posterior probability
  hb_state_RDR <- matrix(2, K, U)
  for(j in 1:length(potentialCnvs_RDR)){
    amp_k <- unlist(lapply(1:K, function(k)  median(probs_RDR[2*(j - 1) + 1, hb_cluster_RDR == k])))
    del_k <- unlist(lapply(1:K, function(k)  median(probs_RDR[2*j, hb_cluster_RDR == k])))
    hb_state_RDR[amp_k >= 0.5, potentialCnvs_RDR[[j]]] <- 3
    hb_state_RDR[del_k >= 0.5, potentialCnvs_RDR[[j]]] <- 1
  }
  
  ret$"states_RDR" <- hb_state_RDR
  
  if(!is.null(A)){
    ### =============================================================== Combined ===============================================================
    hb <- my_setAlleleMats(r.init = A, n.sc.init = D, het.deviance.threshold = 0.1, n.cores = 1, filter = FALSE)
    potential_cnvs <- my_calcAlleleCnvBoundaries(r.maf = A, n.sc = D, verbose = FALSE)
    
    bound.snps.cont <- potential_cnvs$bound.snps.cont
    tbv <- potential_cnvs$tbv
    
    potentialCnvs_BAF <- lapply(names(tbv), function(ti){
      bound.snps.new <- names(bound.snps.cont)[bound.snps.cont == ti]
      bound.snps.new <- as.numeric(bound.snps.new[1:round(length(bound.snps.new)-length(bound.snps.new)*0.1)]) # trim
      bound.snps.new
    })
    
    potentialCnvs_comb <- c(potentialCnvs_RDR, potentialCnvs_BAF)
    potentialCnvs_comb <- potentialCnvs_comb[!duplicated(potentialCnvs_comb)]
    
    probs_BAF <- do.call(rbind, lapply(potentialCnvs_comb, function(x){
      calcCombCnvProb(r.sub = A, 
                      n.sc.sub = D, 
                      l.sub = hb$l.maf, 
                      n.bulk.sub = hb$n.bulk, 
                      gexp.norm.sub = hb_RDR, 
                      mvFit = mvFit, 
                      m = dev, 
                      region = x, 
                      n.iter = 100, 
                      quiet=TRUE)
    }))
    
    if(!is.null(probs_BAF)) probs_BAF <- do.call(rbind, probs_BAF)
    probs_comb <- rbind(probs_RDR, probs_BAF)
    probs_for_clustering_comb <- probs_comb[apply(probs_comb, 1, max) >= 0.45,  , drop = FALSE]
    
    ## cluster cells on posterior probability
    hc_comb <- hclust(dist(t(probs_for_clustering_comb)), method='ward.D')
    hb_cluster_both <- cutree(hc_comb, K)
    ret$"clusters_both" <- hb_cluster_both
    
    ## States on posterior probability
    honey_state_both <- matrix(2, K, U)
    for(j in 1:length(potentialCnvs_comb)){
      amp_k <- unlist(lapply(1:K, function(k)  median(probs_comb[2*(j - 1) + 1, hb_cluster_both = k])))
      del_k <- unlist(lapply(1:K, function(k)  median(probs_comb[2*j, hb_cluster_both = k])))
      honey_state_both[amp_k >= 0.5, potentialCnvs_comb[[j]]] <- 3
      honey_state_both[del_k >= 0.5, potentialCnvs_comb[[j]]] <- 1
    }
    ret$"states_both" <- honey_state_both
  }
  
  ret
}