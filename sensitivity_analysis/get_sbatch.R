sink("sbatch.sh")

Omegas <- c(1:5, 10, 20, 30, 40, 50:60)
Ts <- c(10, 20, 50, 100)

for(w in Omegas){
  for(TT in Ts){
    cmd <- "\nsbatch"
    cmd <- paste0(cmd, " --job-name=w", w, "_T", TT)
    cmd <- paste0(cmd, " --output=log/w", w, "_T", TT, "-out")
    cmd <- paste0(cmd, " --error=log/w", w, "_T", TT, "-error")
    cmd <- paste0(cmd, " run_sims.sh")
    cmd <- paste0(cmd, " ", w, " ", TT)
    cat(cmd)
  }
}

sink()