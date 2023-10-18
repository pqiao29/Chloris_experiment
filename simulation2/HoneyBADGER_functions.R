my_calcGexpCnvProb = function(gexp.norm, mvFit, m=0.15, region=NULL, verbose=FALSE) {
    
    gexp <- gexp.norm
    fits <- mvFit
    quiet <- !verbose
    
    if(!is.null(region)) gexp <- gexp[region, ]
    
    ## smooth
    mu0 <- apply(gexp, 2, mean)
    ng <- nrow(gexp)
    sigma0 <- unlist(lapply(fits, function(fit) sqrt(10^predict(fit, newdata=data.frame(x=ng), interval="predict")[, 'fit'])))
    
    ## Model
    if(verbose) {
        cat('Aggregating data to list ... \n')
    }
    data <- list(
        'K' = length(mu0),
        'JJ' = nrow(gexp),
        'gexp' = gexp,
        'sigma0' = sigma0,
        'mag0' = m
    )
    modelFile <-  system.file("bug", "expressionModel.bug", package = "HoneyBADGER")
    
    if(verbose) {
        cat('Initializing model ... \n')
    }
    ##model <- jags.model(modelFile, data=data, n.chains=4, n.adapt=300, quiet=quiet)
    ##update(model, 1000, progress.bar=ifelse(quiet,"none","text"))
    inits <- list(
        list(S = rep(0, ncol(gexp)), dd = 0),
        list(S = rep(1, ncol(gexp)), dd = 0),
        list(S = rep(0, ncol(gexp)), dd = 1),
        list(S = rep(1, ncol(gexp)), dd = 1)
    )
    model <- rjags::jags.model(modelFile, data=data, inits=inits, n.chains=4, n.adapt=100, quiet=quiet)
    update(model, 100, progress.bar=ifelse(quiet,"none","text"))
    
    parameters <- c('S', 'dd')
    samples <- rjags::coda.samples(model, parameters, n.iter=1000, progress.bar=ifelse(quiet,"none","text"))
    samples <- do.call(rbind, samples) # combine chains
    
    if(verbose) {
        cat('...Done!')
    }
    
    snpLike <- samples
    v <- colnames(snpLike)
    S <- snpLike[,grepl('S', v)]
    dd <- snpLike[,grepl('dd', v)]
    delcall <- apply(S*(1-dd), 2, mean)
    ampcall <- apply(S*dd, 2, mean)
    names(ampcall) <- names(delcall) <- colnames(gexp)
    
    return(list('posterior probability of amplification'=ampcall,
                'posterior probability of deletion'=delcall))
}

my_setGexpMats <- function(gexp.sc, gexp.ref){
    
    if (!("Matrix" %in% class(gexp.ref))) {
        gexp.ref <- as.matrix(gexp.ref)
    }
    
    # cat(paste0("Normalizing gene expression for ", nrow(gexp.sc), 
    #            " genes and ", ncol(gexp.sc), " cells ... \n"))
    refmean <- rowMeans(gexp.ref)
    gexp.norm <- gexp.sc - refmean
    
    return(list(gexp.sc = gexp.sc, gexp.red = gexp.ref, gexp.norm = gexp.norm))
}


my_calcGexpCnvBoundaries <- function (gexp.norm, m = 0.15,  
                                      min.traverse = 3, t = 1e-06, min.num.genes = 3, trim = 0.1, 
                                      verbose = FALSE) {
    # genes <- genes[rownames(gexp.norm)]
    # gos <- as.data.frame(genes)
    # rownames(gos) <- names(genes)
    # gos <- gos[rownames(gexp.norm), ]
    # tl <- tapply(1:nrow(gos), as.factor(gos$seqnames), function(ii) {
    #     na.omit(gexp.norm[rownames(gos)[ii[order((gos[ii, ]$start + 
    #                                                   gos[ii, ]$end)/2, decreasing = F)]], ])
    # })
    # tl <- tl[chrs]
    # gexp.norm <- do.call(rbind, lapply(tl, function(x) x))
    if(is.null(rownames(gexp.norm))) rownames(gexp.norm) <- 1:nrow(gexp.norm)
    k = 101
    mat.smooth <- apply(gexp.norm, 2, caTools::runmean, k)
    d <- dist(t(mat.smooth))
    d[is.na(d)] <- 0
    d[is.nan(d)] <- 0
    d[is.infinite(d)] <- 0
    hc <- hclust(d, method = "ward.D2")
    heights <- seq_len(min(min.traverse, ncol(gexp.norm)))
    boundgenes.pred <- lapply(heights, function(h) {
        ct <- cutree(hc, k = h)
        cuts <- unique(ct)
        boundgenes.pred <- lapply(cuts, function(group) {
            if (sum(ct == group) > 1) {
                mat.smooth <- apply(gexp.norm[, ct == group], 1, mean)
                delta <- c(0, 1, 0)
                t <- t
                pd <- -m
                pn <- 0
                pa <- m
                sd <- sd(mat.smooth)
                z <- HiddenMarkov::dthmm(mat.smooth, matrix(c(1 - 2 * t, t, t, t, 1 - 2 * t, t, t, t, 1 - 2 * 
                                                                  t), byrow = TRUE, nrow = 3), delta, "norm", 
                                         list(mean = c(pd, pn, pa), sd = c(sd, sd, sd)))
                results <- HiddenMarkov::Viterbi(z)
                ampgenes <- which(results ==  3)
                delgenes <- which(results ==  1)
                boundgenes <- list(amp = ampgenes, del = delgenes)
                return(boundgenes)
            }
        })
    })
    boundgenes.pred <- unlist(boundgenes.pred, recursive = FALSE)
    getTbv <- function(boundgenes.pred) {
        foo <- rep(0, nrow(gexp.norm))
        names(foo) <- rownames(gexp.norm)
        foo[unique(unlist(boundgenes.pred))] <- 1
        vote <- rep(0, nrow(gexp.norm))
        names(vote) <- rownames(gexp.norm)
        lapply(boundgenes.pred, function(b) {
            vote[b] <<- vote[b] + 1
        })
        if (verbose) {
            cat(paste0("max vote:", max(vote), "\n"))
        }
        if (max(vote) == 0) {
            if (verbose) {
                cat("Exiting; no genes found.\n")
            }
            return()
        }
        vote[vote > 0] <- 1
        mv <- 1
        cs <- 1
        bound.genes.cont <- rep(0, length(vote))
        names(bound.genes.cont) <- names(vote)
        for (i in 2:length(vote)) {
            if (vote[i] >= mv & vote[i] == vote[i - 1]) {
                bound.genes.cont[i] <- cs
            }
            else {
                cs <- cs + 1
            }
        }
        tb <- table(bound.genes.cont)
        tbv <- as.vector(tb)
        names(tbv) <- names(tb)
        tbv <- tbv[-1]
        tbv[tbv < min.num.genes] <- NA
        tbv <- na.omit(tbv)
        if (length(tbv) == 0) {
            if (verbose) {
                cat(paste0("Exiting; fewer than ", min.num.genes, 
                           " new bound genes found.\n"))
            }
            return()
        }
        boundgenes.info <- lapply(names(tbv), function(ti) {
            bound.genes.new <- names(bound.genes.cont)[bound.genes.cont == ti]
            bound.genes.new <- bound.genes.new[1:round(length(bound.genes.new) - 
                                                           length(bound.genes.new) * trim)]
        })
        
        return(boundgenes.info)
    }
    amp.info <- getTbv(lapply(boundgenes.pred, function(x) x[["amp"]]))
    del.info <- getTbv(lapply(boundgenes.pred, function(x) x[["del"]]))
    
    return(list(amp = amp.info, del = del.info))
}

my_setAlleleMats <- function(r.init, n.sc.init, l.init = NULL, n.bulk.init = NULL, filter = TRUE, het.deviance.threshold = 0.05, min.cell = 3, n.cores = 1, verbose = FALSE) {
    
    if(verbose) cat("Initializing allele matrices ... \n")
    
    if(is.null(l.init) | is.null(n.bulk.init)) {
        if(verbose) {
            cat("Creating in-silico bulk ... \n")
            cat(paste0("using ", ncol(r.init), " cells ... \n"))
        }
        l <- rowSums(r.init > 0)
        n.bulk <- rowSums(n.sc.init > 0)
    } else {
        l <- l.init
        n.bulk <- n.bulk.init
    }
    
    if(filter) {
        if(verbose) {
            cat("Filtering for putative heterozygous snps ... \n")
            cat(paste0("allowing for a ", het.deviance.threshold, " deviation from the expected 0.5 heterozygous allele fraction ... \n"))
        }
        E <- l/n.bulk
        vi <- names(which(E > het.deviance.threshold & E < 1-het.deviance.threshold))
        if(length(vi) < 0.01*length(l)) {
            cat("WARNING! CLONAL DELETION OR LOH POSSIBLE! \n")
        }
        r <- r.init[vi,]
        n.sc <- n.sc.init[vi,]
        l <- l[vi]
        n.bulk <- n.bulk[vi]
        
        ## must have coverage in at least 5 cells
        if(verbose) cat(paste0("must have coverage in at least ", min.cell, " cells ... \n"))
        vi <- rowSums(n.sc > 0) >= min.cell
        cat(paste0(length(vi), " heterozygous SNPs identified \n"))
        r <- r[vi,]
        n.sc <- n.sc[vi,]
        l <- l[vi]
        n.bulk <- n.bulk[vi]
    } else {
        r <- r.init
        n.sc <- n.sc.init
    }
    if(!is.null(r.init) | !is.null(l.init)) {
        if(verbose) cat("Setting composite lesser allele count ... \n")
        E <- l/n.bulk
        n <- nrow(r)
        m <- ncol(r)
        mat <- do.call(rbind, parallel::mclapply(1:n, function(i) {
            do.call(cbind, lapply(1:m, function(j) {
                ri <- r[i,j]
                n.sci <- n.sc[i,j]
                Ei <- E[i]
                if(is.na(Ei)) {
                    mut.frac <- 0
                }
                else if(Ei <= 0.5) {
                    mut.frac <- ri
                }
                else if(Ei > 0.5) {
                    mut.frac <- n.sci-ri
                }
                else {
                    mut.frac <- 0
                }
                
                ## f will be high if inconsistent
                ## f will be low if consistent
                ## f will be NaN if no coverage
                ## use colorRamp from green to red
                f <- mut.frac
                return(f)
            }))
        }, mc.cores=n.cores))
        r.maf <- mat
        
        mat <- sapply(1:n, function(i) {
            li <- l[i]
            n.bulki <- n.bulk[i]
            Ei <- E[i]
            if(is.na(Ei)) {
                mut.frac <- 0
            }
            else if(Ei <= 0.5) {
                mut.frac <- li
            }
            else if(Ei > 0.5) {
                mut.frac <- n.bulki-li
            }
            else {
                mut.frac <- 0
            }
            
            ## f will be high if inconsistent
            ## f will be low if consistent
            ## f will be NaN if no coverage
            ## use colorRamp from green to red
            f <- mut.frac
            return(f)
        })
        l.maf <- mat
    }
    
    if(verbose) cat("Done setting initial allele matrices! \n")
    
    return(list(r = r, r.maf = r.maf, n.sc = n.sc, 
                l = l, l.maf = l.maf, n.bulk = n.bulk))
}

my_calcAlleleCnvBoundaries <- function(r.maf, n.sc, min.traverse = 5, t=1e-6, pd=0.1, pn=0.45, min.num.snps=3, verbose=FALSE) {
    
    if(is.null(colnames(r.maf))) colnames(r.maf) <- 1:ncol(r.maf)
    if(is.null(rownames(r.maf))) rownames(r.maf) <- 1:nrow(r.maf)
    if(is.null(colnames(n.sc))) colnames(n.sc) <- 1:ncol(n.sc)
    if(is.null(rownames(n.sc))) rownames(n.sc) <- 1:nrow(n.sc)
    bound.snps.final <- NULL
    
    pred.snps.r <- matrix(0, nrow(r.maf), ncol(r.maf))
    rownames(pred.snps.r) <- rownames(r.maf)
    colnames(pred.snps.r) <- colnames(r.maf)
    bound.snps.final <- list()
    
    if(verbose) cat('ignore previously identified CNVs ... ')
    
    ## lesser allele fraction
    mat.tot <- r.maf/n.sc
    
    mat.smooth <- apply(mat.tot, 2, caTools::runmean, k=31)
    d <- dist(t(mat.smooth))
    d[is.na(d)] <- 0
    d[is.nan(d)] <- 0
    d[is.infinite(d)] <- 0
    hc <- hclust(d, method="ward.D2")
    
    if(verbose) cat('iterative HMM ... ')
    
    ## iterative HMM
    heights <- 1:min(min.traverse, ncol(r.maf))
    ## cut tree at various heights to establish groups
    boundsnps.pred <- lapply(heights, function(h) {
        
        ct <- cutree(hc, k = h)
        
        cuts <- unique(ct)
        
        ## look at each group, if deletion present
        boundsnps.pred <- lapply(cuts, function(group) {
            if(sum(ct==group)>1) {
                mafl <- rowSums(r.maf[, ct==group]) 
                sizel <- rowSums(n.sc[, ct==group])
                
                ## change point
                delta <- c(0, 1)
                z <- HiddenMarkov::dthmm(mafl, matrix(c(1-t, t, t, 1-t), byrow=TRUE, nrow=2), delta, "binom", list(prob=c(pd, pn)), list(size=sizel), discrete=TRUE)
                results <- HiddenMarkov::Viterbi(z)
                
                ## Get boundaries from states
                boundsnps <- rownames(r.maf)[results == 1]
                return(boundsnps)
            }
        })
    })
    
    foo <- rep(0, nrow(r.maf)); names(foo) <- rownames(r.maf)
    foo[unique(unlist(boundsnps.pred))] <- 1
    ## vote
    vote <- rep(0, nrow(r.maf))
    names(vote) <- rownames(r.maf)
    lapply(boundsnps.pred, function(b) {
        vote[b[[1]]] <<- vote[b[[1]]] + 1
    })
    
    if(verbose) {
        cat(paste0('max vote:', max(vote), '\n'))
    }
    
    if(max(vote)==0) {
        if(verbose) {
            cat('Exiting; no new bound SNPs found.\n')
        }
        return() ## exit iteration, no more bound SNPs found
    }
    
    vote[vote > 0] <- 1
    mv <- 1 ## at least 1 vote
    cs <- 1
    bound.snps.cont <- rep(0, length(vote))
    names(bound.snps.cont) <- names(vote)
    for(i in 2:length(vote)) {
        if(vote[i] >= mv & vote[i] == vote[i-1]) {
            bound.snps.cont[i] <- cs
        } else {
            cs <- cs + 1
        }
    }
    tb <- table(bound.snps.cont)
    tbv <- as.vector(tb); names(tbv) <- names(tb)
    tbv <- tbv[-1] # get rid of 0
    
    ## all detected deletions have fewer than 5 SNPs...reached the end
    tbv[tbv < min.num.snps] <- NA
    tbv <- na.omit(tbv)
    if(length(tbv)==0) {
        if(verbose) {
            cat(paste0('Exiting; less than ', min.num.snps, ' new bound SNPs found.\n'))
        }
        return()
    }
    
    return(list("tbv" = tbv, "bound.snps.cont" = bound.snps.cont))
}


my_calcAlleleCnvProb <- function(r.sub=NULL, n.sc.sub=NULL, l.sub=NULL, n.bulk.sub=NULL, pe=0.1, mono=0.7, n.iter=1000, quiet=TRUE ) {
    
    ## Convert to multi-dimensions based on j
    I.j <- rep(1, nrow(r.sub))#unlist(lapply(genes2snps.dict, length))
    numGenes <- nrow(r.sub)
    numSnpsPerGene <- 1
    numCells <- ncol(r.sub)
    ## j, i, k
    r.array <- array(0, c(numGenes, numSnpsPerGene, numCells))
    for(i in seq_len(numGenes)) {
        r.array[i, , ] <- r.sub[i, ]
    }
    n.sc.array <- array(0, c(numGenes, numSnpsPerGene, numCells))
    for(i in seq_len(numGenes)) {
        n.sc.array[i, , ] <- n.sc.sub[i, ]
    }
    l.array <- array(0, c(numGenes, numSnpsPerGene))
    for(i in seq_len(numGenes)) {
        l.array[i, ] <- l.sub[i]
    }
    n.bulk.array <- array(0, c(numGenes, numSnpsPerGene))
    for(i in seq_len(numGenes)) {
        n.bulk.array[i, ] <- n.bulk.sub[i]
    }
    
    data <- list(
        'l' = l.array,
        'r' = r.array,
        'n.bulk' = n.bulk.array,
        'n.sc' = n.sc.array,
        'J' = length(I.j),  # how many genes
        'K' = ncol(r.sub),  # how many cells
        'I.j' = I.j,
        'pseudo' = pe,
        'mono' = mono)
    
    modelFile <- system.file("bug", "snpModel.bug", package = "HoneyBADGER")
    
    require(rjags)
    model <- rjags::jags.model(modelFile, data=data, n.chains=4, n.adapt=300, quiet=quiet)
    update(model, 300, progress.bar=ifelse(quiet,"none","text"))
    
    parameters <- 'S'
    samples <- rjags::coda.samples(model, parameters, n.iter=n.iter, progress.bar=ifelse(quiet,"none","text"))
    samples <- do.call(rbind, samples) # combine samples across chains
    pm <- do.call(cbind,lapply(seq_len(numCells),function(ci) {
        c(mean(samples[,paste("S[",ci,"]",sep="")]))
    }))
    
    return(pm)
}

calcCombCnvProb=function(r.sub, n.sc.sub, l.sub, n.bulk.sub, gexp.norm.sub, mvFit = NULL, m=0.15, region=NULL, pe=0.1, mono=0.7, n.iter=1000, quiet=FALSE, verbose=FALSE) {
    
    gexp <- gexp.norm.sub
    fits <- mvFit
    
    if(!is.null(region)) {
        gexp <- gexp[region, ]
        r.sub <- r.sub[region,]
        n.sc.sub <- n.sc.sub[region,]
        l.sub <- l.sub[region]
        n.bulk.sub <- n.bulk.sub[region]
    }
    
    ## smooth
    mu0 <- apply(gexp, 2, mean)
    ng <- nrow(gexp)
    sigma0 <- unlist(lapply(fits, function(fit) sqrt(10^predict(fit, newdata=data.frame(x=ng), interval="predict")[, 'fit'])))
    
    ## Convert to multi-dimensions based on j
    I.j <- rep(1, nrow(r.sub))
    numGenes <- nrow(r.sub)
    numSnpsPerGene <- max(I.j)
    numCells <- ncol(r.sub)
    ## j, i, k
    r.array <- array(0, c(numGenes, numSnpsPerGene, numCells))
    for(i in seq_len(numGenes)) {
        r.array[i, , ] <- r.sub[i, ]
    }
    n.sc.array <- array(0, c(numGenes, numSnpsPerGene, numCells))
    for(i in seq_len(numGenes)) {
        n.sc.array[i, , ] <- n.sc.sub[i, ]
    }
    l.array <- array(0, c(numGenes, numSnpsPerGene))
    for(i in seq_len(numGenes)) {
        l.array[i, ] <- l.sub[i]
    }
    n.bulk.array <- array(0, c(numGenes, numSnpsPerGene))
    for(i in seq_len(numGenes)) {
        n.bulk.array[i, ] <- n.bulk.sub[i]
    }
    
    
    data <- list(
        'l' = l.array,
        'r' = r.array,
        'n.bulk' = n.bulk.array,
        'n.sc' = n.sc.array,
        'J' = length(I.j),  # how many genes
        'K' = ncol(r.sub),  # how many cells
        'I.j' = I.j,
        'pseudo' = pe,
        'mono' = mono,
        'gexp' = gexp,
        'JJ' = nrow(gexp),
        'sigma0' = sigma0,
        'mag0' = m
    )
    
    modelFile <- system.file("bug", "combinedModel.bug", package = "HoneyBADGER")
    
    model <- rjags::jags.model(modelFile, data=data, n.chains=4, n.adapt=300, quiet=quiet)
    update(model, 300, progress.bar=ifelse(quiet,"none","text"))
    
    parameters <- c('S', 'dd')
    samples <- rjags::coda.samples(model, parameters, n.iter=300, progress.bar=ifelse(quiet,"none","text"))
    samples <- do.call(rbind, samples) # combine chains
    
    snpLike <- samples
    v <- colnames(snpLike)
    S <- snpLike[,grepl('S', v)]
    dd <- snpLike[,grepl('dd', v)]
    
    delcall <- apply(S*(1-dd), 2, mean)
    ampcall <- apply(S*dd, 2, mean)
    
    return(list('posterior probability of amplification'=ampcall,
                'posterior probability of deletion'=delcall)
    )
}