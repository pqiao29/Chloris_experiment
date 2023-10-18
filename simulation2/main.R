args <- commandArgs(trailingOnly = TRUE)
sim_runs <- as.numeric(args[[1]])
ncores <- as.numeric(args[[2]])
batch <- as.numeric(args[[3]])

import::from("infercnv_HMM.R", infercnv_HMM)
import::from("infercnv_clustering_randomtree.R", single_tumor_subclustering_recursive_random_smoothed_trees)
import::from("infercnv_emission_pars.R", infercnv_emission_pars)
import::from("infercnv_main.R", run_infercnv)
source("HoneyBADGER_functions.R")
import::from("HoneyBADGER.R", run_HoneyBADGER)
source("../simulation1/util.R") #state_align
library(Chloris)



library(doParallel)
registerDoParallel(ncores) 
import::from(foreach, "%dopar%", foreach)
import::from("simulation.R", compare)

set.seed(1)
res_fin <- foreach (seed = sample(1:1000000, sim_runs, replace = F)) %dopar% {
  try(compare(seed))
}

saveRDS(res_fin, file = paste0("result/compare.rds"))
