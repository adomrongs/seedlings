# Fill the breeding pipeline so that every stage holds a different cohort.
# After the loop: F1 comes from cohort 5, S1 from cohort 4, S2 from cohort 3,
# S3 from cohort 2 and S4 from cohort 1, which is the steady state that
# 4_AdvanceYear.R expects.

Pfill <- runif(5)

for (cohort in 1:5) {
  F1 <- randCross(parents, nCrosses = nCrosses, nProgeny = nProgeny)

  if (cohort < 5) {
    S1 <- setPheno(F1, h2 = h2S1, p = Pfill[cohort])
  }

  if (cohort < 4) {
    S2 <- selectInd(S1, nInd = nS2, use = 'pheno')
    S2 <- setPheno(S2, h2 = h2S2, p = Pfill[cohort])
  }

  if (cohort < 3) {
    S3 <- selectInd(S2, nInd = nS3, use = 'pheno')
    S3 <- setPheno(S3, h2 = h2S3, p = Pfill[cohort])
  }

  if (cohort < 2) {
    S4 <- selectInd(S3, nInd = nS4, use = 'pheno')
    S4 <- setPheno(S4, h2 = h2S4, p = Pfill[cohort])
  }
}

cat("Pipeline ready\n")
