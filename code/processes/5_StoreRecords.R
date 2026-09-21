# Store the phenotypic records that feed the training population and keep the
# pedigree of everything that can ever become a parent or a candidate.

ped.file <- paste0('outputs/', scenario, '/ped_', rep, '.csv')
if (!dir.exists(dirname(ped.file))) {
  dir.create(dirname(ped.file), recursive = TRUE)
}

# ----------- Phenotypic records of S3 and S4 -----------
newRecords <- rbind(
  data.frame(id = S3@id, year = year, stage = 'S3', pheno = S3@pheno[, 1]),
  data.frame(id = S4@id, year = year, stage = 'S4', pheno = S4@pheno[, 1])
)

# ----------- Pedigree of the current parents -----------
# Founders come out of getPed() with mother and father set to 0; they have to
# become NA before the rows are stacked, because doing it afterwards would mean
# subsetting a data frame that already contains NAs.
newPed <- getPed(parents)
newPed[newPed == 0] <- NA

if (year == 1) {
  records <- newRecords
  tpPop <- c(S3, S4)
  ped <- newPed
} else {
  records <- rbind(records, newRecords)
  tpPop <- c(tpPop, S3, S4)
  ped <- bind_rows(ped, newPed)
}

# ----------- Keep only the last nRecords years -----------
records <- records[records$year > (year - nRecords), ]
tpPop <- tpPop[!duplicated(tpPop@id)]
tpPop <- tpPop[tpPop@id %in% records$id]

# ----------- Pedigree of the whole simulation -----------
ped <- ped[!duplicated(ped$id), ]
write.csv(ped, file = ped.file, row.names = FALSE)

cat("    records:", nrow(records), "| training individuals:", tpPop@nInd, "\n")
