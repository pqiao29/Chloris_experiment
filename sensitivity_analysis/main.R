args <- commandArgs(trailingOnly = TRUE)
min_cluster_size <- as.numeric(args[[1]])
cluster_shrink_tol <- as.numeric(args[[2]])
ncores <- as.numeric(args[[3]])
sim_runs <- as.numeric(args[[4]])

library(Chloris)

library(doParallel)
registerDoParallel(ncores) 
import::from(foreach, "%dopar%", foreach)
import::from("simulation.R", run)

set.seed(1)
res_fin <- foreach (seed = sample(1:1000000, sim_runs, replace = F)) %dopar% {
  try(run(min_cluster_size, cluster_shrink_tol, seed))
}

saveRDS(res_fin, file = paste0("result/w", min_cluster_size, "T", cluster_shrink_tol, ".rds"))
