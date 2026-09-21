# Fit a GBLUP model on the training population (S3 + S4 records of the last
# nRecords years) and predict GEBVs for every genotyped candidate.
#
# Marker effects are back-solved from the training GEBVs so that (i) candidates
# that are not in the model can be predicted and (ii) COMA can be fed directly.

# Everything genotyped gets a GEBV, whatever the scenario. This is wider than
# the mating candidate pool on purpose: the S2 -> S3 selection is genomic in
# every future scenario, so S2 has to be predicted even in COMA_S4 and
# COMA_S3S4, where S2 is not a candidate for crossing.
allGenotyped <- c(parents, S2, S3, S4)
allGenotyped <- allGenotyped[!duplicated(allGenotyped@id)]

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

# ----------- Predict everything genotyped -----------
Mall <- pullSnpGeno(allGenotyped)[, colnames(Mtrain), drop = FALSE]
Wall <- sweep(Mall, 2, ploidy * freq)
allGenotyped@ebv <- matrix(Wall %*% ameff, ncol = 1)

parents <- setEBV(parents, allGenotyped)
S2 <- setEBV(S2, allGenotyped)
S3 <- setEBV(S3, allGenotyped)
S4 <- setEBV(S4, allGenotyped)

# ----------- Mating candidate pool -----------
# Each scenario declares it in poolStages. 'parents' is the crossing block of
# last year's mating plan: still available for crossing and, in most years, no
# longer present in any current stage, so it is what keeps generations
# overlapping. COMA_onlyS2 leaves it out on purpose, as the seedlings-only
# control. Built after setEBV so the pool carries the GEBVs.
stagePops <- list(parents = parents, S2 = S2, S3 = S3, S4 = S4)
candidatePop <- do.call(c, unname(stagePops[poolStages]))
candidatePop <- candidatePop[!duplicated(candidatePop@id)]

Mcand <- pullSnpGeno(candidatePop)[, colnames(Mtrain), drop = FALSE]

nCandidates <- candidatePop@nInd
accuracy_S2 <- cor(S2@gv, S2@ebv)
accuracy_cand <- cor(candidatePop@gv, candidatePop@ebv)

cat(
  "    GEBVs ready | accuracy in S2:",
  round(accuracy_S2, 3),
  "| candidates:",
  nCandidates,
  "\n"
)
