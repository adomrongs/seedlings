# Write the two files COMA::read_data() needs: marker effects plus dosages of
# the candidate pool, and the pedigree kinship matrix of that same pool.

geno.file <- paste0('outputs/', scenario, '/geno_', rep, '.csv.gz')
k.file <- paste0('outputs/', scenario, '/k_', rep, '.csv.gz')
ped.file2 <- paste0('outputs/', scenario, '/ped_cand_', rep, '.csv')

# ----------- Pedigree kinship of the candidates -----------
candPed <- getPed(candidatePop)
candPed[candPed == 0] <- NA
combined_ped <- bind_rows(ped, candPed)
combined_ped <- combined_ped[!duplicated(combined_ped$id), ]
write.csv(combined_ped, file = ped.file2, row.names = FALSE)

Acand <- getA(ped.file2, candidatePop)
suppressMessages(fwrite(
  as.data.frame(Acand / ploidy),
  file = k.file,
  row.names = TRUE
))

# ----------- Marker effects and dosages -----------
geno <- data.frame(markers = colnames(Mtrain), add = ameff)
geno <- cbind(geno, t(Mcand))
suppressMessages(fwrite(geno, file = geno.file, row.names = FALSE))

cat("    COMA files written for", candidatePop@nInd, "candidates\n")
