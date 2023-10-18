sink("sbatch.sh")

RDR_var <- seq(0.5, 5, 0.5)
K_range <- 4:7
ncores <- 20
nsims <- 600

for(RDR_var in RDR_var){
    for(K in K_range){
        cmd <- "\nsbatch"
        cmd <- paste0(cmd, " --job-name=var", RDR_var*10, "_K", K)
        cmd <- paste0(cmd, " --output=log/var", RDR_var*10, "_K", K, "-out")
        cmd <- paste0(cmd, " --error=log/var", RDR_var*10, "_K", K, "-error")
        cmd <- paste0(cmd, " run_sims.sh")
        cmd <- paste0(cmd, " ", RDR_var, " ", K, " ", ncores, " ", nsims)
        cat(cmd)
    }
}

sink()
