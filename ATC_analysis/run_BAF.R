load("ATC_BAF.Rdata")
library(Chloris)

res <- Chloris(A = B, D = D, S = 2, init = "random")

saveRDS(res, file = "Chloris_BAFref.rds")