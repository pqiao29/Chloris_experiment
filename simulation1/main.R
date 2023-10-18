args <- commandArgs(trailingOnly = TRUE)
RDR_var <- as.numeric(args[[1]])
K <- as.numeric(args[[2]])
ncores <- as.numeric(args[[3]])
sim_runs <- as.numeric(args[[4]])


library(doParallel)
registerDoParallel(ncores) 
library(Chloris)
import::from(foreach, "%dopar%", foreach)
import::from("simulation.R", run)
import::from("util.R", state_align)

set.seed(0)
res_fin <- foreach (seed = sample(1:1000000, sim_runs, replace = F)) %dopar% {
  try(run(RDR_var, K, seed, 200, 200))
}

saveRDS(res_fin, file = paste0("result/RDRvar", RDR_var*10, "_K", K, ".rds"))
