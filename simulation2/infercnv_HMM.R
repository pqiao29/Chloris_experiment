infercnv_HMM <- function(infercnv_clust, mu, sd, neutral_idx, obs, K){
    
    state_emission_params = list(mean = mu, sd = sd)
    
    S <- length(mu)
    
    Q_off_d = 1e-6
    Q <- matrix(Q_off_d, S, S)
    diag(Q) <- 1- (S - 1)*Q_off_d
    
    delta = rep(Q_off_d, S)
    delta[neutral_idx] <- 1 - (S - 1)*Q_off_d
    
    HMM_info = list(state_transitions = Q, delta = delta)
    
    states <- matrix(NA, nrow(obs), K)
    hmm.data = obs
    hmm.data[,] = -1 #init to invalid state
    
    for(k in 1:K){
        
        cluster_idx <- which(infercnv_clust == k)
        cluster_expr = rowMeans(obs[ , cluster_idx, drop=FALSE])
        
        
        hmm <- HiddenMarkov::dthmm(cluster_expr,
                                   HMM_info[['state_transitions']],
                                   HMM_info[['delta']],
                                   "norm",
                                   state_emission_params)
        
        # hmm_trace <- Viterbi.dthmm.adj(hmm)
        # ---------------------------------------------------------------------------------------
        object <- hmm
        x <- object$x
        if (length(x) < 2) stop("Markov Chain is too short")
        dfunc <- HiddenMarkov:::makedensity(object$distn)
        n <- length(x)
        m <- nrow(object$Pi) # transition matrix
        nu <- matrix(NA, nrow = n, ncol = m)  # scoring matrix
        hmm_trace <- rep(NA, n) # final trace
        pseudocount = 1e-20
        emissions <- matrix(NA, nrow = n, ncol = m) 
        object$pm$sd = median(object$pm$sd) ## restrict to constant variance to avoid nonsensical results
        
        ## init first row
        emission <- pnorm(abs(x[1] - object$pm$mean)/object$pm$sd, log.p=TRUE, lower.tail=FALSE)
        emission <- 1 / (-1 * emission)
        emission <- emission / sum(emission)
        emissions[1,] <- log(emission)
        nu[1, ] <- log(object$delta) + emissions[1, ] # start probabilities
        logPi <- log(object$Pi) # convert transition matrix to log(p)
        
        for (i in 2:n) {
            matrixnu <- matrix(nu[i - 1, ], nrow = m, ncol = m)
            emission <- pnorm(abs(x[i]-object$pm$mean)/object$pm$sd, log.p=TRUE, lower.tail=FALSE)
            emission <- 1 / (-1 * emission)
            emission <- emission / sum(emission)
            emissions[i, ] <- log(emission)
            nu[i, ] <- apply(matrixnu + logPi, 2, max) + emissions[i, ] 
        }
        if (any(nu[n, ] == -Inf)) stop("Problems With Underflow")
        
        ## traceback
        hmm_trace[n] <- which.max(nu[n, ])
        for (i in seq(n - 1, 1, -1)) hmm_trace[i] <- which.max(logPi[, hmm_trace[i + 1]] + nu[i, ])
        # --------------------------------------------------------------------------------------- (End of Viterbi.dthmm.adj)
        
        hmm.data[, cluster_idx] <- hmm_trace
        states[, k] <- hmm_trace
    }
    
    return(states)
}
