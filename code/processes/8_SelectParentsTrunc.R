# Control GS scenario: the crossing block is the top nParents candidates by
# GEBV, taken from the whole genotyped pool (S2 + S3 + S4), and crosses are
# random among them.

old <- parents
parents <- selectInd(candidatePop, nInd = nParents, use = 'ebv')
nKept <- sum(old@id %in% parents@id)

F1_tmp <- randCross(parents, nCrosses = nCrosses, nProgeny = nProgeny)

# ----------- Share of the mating plan coming from seedlings -----------
matings <- data.frame(
  parent1 = F1_tmp@mother,
  parent2 = F1_tmp@father,
  stringsAsFactors = FALSE
)
matings <- matings[!duplicated(matings), ]

share <- seedlingShare(
  parentIds = parents@id,
  matings = matings,
  s2ids = S2@id,
  progeny = rep(nProgeny, nrow(matings))
)

nParentsUsed <- parents@nInd
nCrossesUsed <- nrow(matings)
time <- NA_real_
