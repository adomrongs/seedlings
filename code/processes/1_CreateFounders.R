# Create the founder population, the trait and the crossing block.

founderpop <- runMacs2(
  nInd = nFounders,
  nChr = nChr,
  segSites = segSites,
  ploidy = ploidy,
  bp = bp,
  genLen = genLen,
  mutRate = mutRate,
  histNe = NULL,
  histGen = NULL
)

SP <- SimParam$new(founderpop)

# ----------- SNP chip -----------
SP$restrSegSites(minQtlPerChr = nQtlPerChr, minSnpPerChr = nSnpPerChr)
SP$addSnpChip(nSnpPerChr = nSnpPerChr)

# ----------- Additive trait -----------
SP$addTraitA(nQtlPerChr = nQtlPerChr, mean = meanG0, var = varG0)
SP$setTrackPed(TRUE)

# ----------- Crossing block -----------
parents <- newPop(founderpop)
parents <- setPheno(parents, h2 = h2S4)

# ----------- Year effects, shared by the whole simulation -----------
P <- runif(nBurnin + nFuture)

cat("Founders ready\n")
