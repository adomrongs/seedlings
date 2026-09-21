# Pedigree inbreeding of the crossing block and of the new F1, and the
# resulting inbreeding rate. Runs after 5_StoreRecords.R, which is what puts
# the current parents in the pedigree file.

Aparents <- getA(ped.file, parents)
Ft0 <- getFt(Aparents, ploidy)
Ft1 <- getFtProgeny(Aparents, parents, F1)
dF1 <- 100 * (Ft1 - Ft0) / (1 - Ft0)
