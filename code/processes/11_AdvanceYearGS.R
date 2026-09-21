# Advance the program one year with genomic selection from S2 to S3.
# Every other selection stays phenotypic. F1 comes from the mating plan built
# in 8_SelectParentsTrunc.R or 10_RunOMA.R.

S4 <- selectInd(S3, nInd = nS4, use = 'pheno')
S4 <- setPheno(S4, h2 = h2S4, p = P[year])

S3 <- selectInd(S2, nInd = nS3, use = 'ebv')
S3 <- setPheno(S3, h2 = h2S3, p = P[year])

S2 <- selectInd(S1, nInd = nS2, use = 'pheno')
S2 <- setPheno(S2, h2 = h2S2, p = P[year])

S1 <- setPheno(F1, h2 = h2S1, p = P[year])

F1 <- F1_tmp
