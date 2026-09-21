# Advance the program one year under phenotypic selection at every stage.
# Stages are advanced from the oldest to the youngest so that each object is
# still holding last year's cohort when it is used.

accuracy_S2 <- cor(S2@gv, S2@pheno)

S4 <- selectInd(S3, nInd = nS4, use = 'pheno')
S4 <- setPheno(S4, h2 = h2S4, p = P[year])

S3 <- selectInd(S2, nInd = nS3, use = 'pheno')
S3 <- setPheno(S3, h2 = h2S3, p = P[year])

S2 <- selectInd(S1, nInd = nS2, use = 'pheno')
S2 <- setPheno(S2, h2 = h2S2, p = P[year])

S1 <- setPheno(F1, h2 = h2S1, p = P[year])

F1 <- randCross(parents, nCrosses = nCrosses, nProgeny = nProgeny)
