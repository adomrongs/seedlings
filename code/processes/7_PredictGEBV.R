# Fit a GBLUP model on the training population (S3 + S4 records of the last
# nRecords years) and predict GEBVs for every genotyped candidate.
#
# Marker effects are back-solved from the training GEBVs so that (i) candidates
# that are not in the model can be predicted and (ii) COMA can be fed directly.

# The candidate pool is the current S2, S3 and S4 plus the parents of last
# year's mating plan. Those parents are still available for crossing, and most
# of them are no longer in any current stage: a parent taken from S2 last year
# is only in S3 now if it advanced, and is otherwise gone. Leaving them out
# would forbid reusing a parent across years and would remove the overlapping
# generations that the inbreeding constraint has to arbitrate.
candidatePop <- c(parents, S2, S3, S4)
candidatePop <- candidatePop[!duplicated(candidatePop@id)]

# ----------- Marker matrix of the training population -----------
Mtrain <- pullSnpGeno(tpPop)
freq <- colMeans(Mtrain) / ploidy
keep <- freq > 0.01 & freq < 0.99
Mtrain <- Mtrain[, keep, drop = FALSE]
freq <- freq[keep]
Wtrain <- sweep(Mtrain, 2, ploidy * freq)

# ----------- Genomic relationship matrix -----------
G <- AGHmatrix::Gmatrix(Mtrain, method = 'VanRaden', ploidy = ploidy, maf = 0)
G <- G + diag(1e-05, nrow(G))

# ----------- Fit the model -----------
model_df <- data.frame(
  gid = factor(records$id, levels = rownames(G)),
  year = factor(records$year),
  stage = factor(records$stage),
  pheno = records$pheno
)

model <- lme4breeding::lmebreed(
  pheno ~ year + stage + (1 | gid),
  data = model_df,
  relmat = list(gid = G)
)

# ----------- GEBVs of the training population -----------
gebvTrain <- ranef(model)$gid[, 1]
names(gebvTrain) <- rownames(ranef(model)$gid)
gebvTrain <- gebvTrain[rownames(G)]

# ----------- Back-solve marker effects -----------
ameff <- crossprod(
  Wtrain,
  solve(tcrossprod(Wtrain) + diag(1e-06, nrow(Wtrain)), gebvTrain)
)
ameff <- as.vector(ameff)

# ----------- Predict the candidates -----------
Mcand <- pullSnpGeno(candidatePop)[, colnames(Mtrain), drop = FALSE]
Wcand <- sweep(Mcand, 2, ploidy * freq)
candidatePop@ebv <- matrix(Wcand %*% ameff, ncol = 1)

S2 <- setEBV(S2, candidatePop)
S3 <- setEBV(S3, candidatePop)
S4 <- setEBV(S4, candidatePop)

nCandidates <- candidatePop@nInd
accuracy_S2 <- cor(S2@gv, S2@ebv)
accuracy_cand <- cor(candidatePop@gv, candidatePop@ebv)

cat("    GEBVs ready | accuracy in S2:", round(accuracy_S2, 3), "\n")
