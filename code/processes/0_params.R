# Global parameters and helper functions for the seedlings simulation
# Sourced at the top of every scenario.

# ----------- Simulation parameters -----------
nBurnin <- 20 # years of phenotypic-only burnin
nFuture <- 30 # years of the future scenario
nRecords <- 8 # years of S3 + S4 records kept in the training population

# ----------- Genome parameters (Covarrubias et al. 2026, diploid outbred) -----------
nFounders <- 40
nChr <- 10
segSites <- 1500 # segregating sites per chromosome
nQtlPerChr <- 100
nSnpPerChr <- 1000
ploidy <- 2
bp <- 8e08 # base pairs per chromosome
genLen <- 1.43 # Morgans per chromosome
mutRate <- 2e-09

# ----------- Trait parameters -----------
meanG0 <- 0
varG0 <- 1

# ----------- Breeding program parameters -----------
nParents <- 40 # size of the crossing block
nCrosses <- 100
nProgeny <- 200 # -> 20000 seedlings per year
nS2 <- 1000
nS3 <- 200
nS4 <- 20
nRecycled <- 10 # S4 clones entering the crossing block every year (burnin / PS)

h2S1 <- 0.1
h2S2 <- 0.2
h2S3 <- 0.4
h2S4 <- 0.8

# ----------- COMA parameters -----------
dF <- 0.01 # target inbreeding rate per year
totalProgeny <- nCrosses * nProgeny
minProgeny <- 200

# ----------- Helper functions -----------

## getA: numerator relationship matrix of a population, from a pedigree on disk.
getA <- function(ped.file, pop) {
  tmp <- polyBreedR::get_pedigree(pop@id, ped.file, na.string = "NA")
  polyBreedR::A_mat(tmp, ploidy = pop@ploidy)[pop@id, pop@id]
}

## getFt: mean pedigree inbreeding of a population.
getFt <- function(A, ploidy) {
  (mean(diag(A)) - 1) / (ploidy - 1)
}

## getFtProgeny: mean inbreeding of a progeny population, computed as the mean
## kinship of the parent pairs that produced it. Avoids building a relationship
## matrix for the 20000 seedlings.
getFtProgeny <- function(A, parentPop, progeny) {
  idx <- cbind(
    match(progeny@mother, parentPop@id),
    match(progeny@father, parentPop@id)
  )
  mean(A[idx]) / parentPop@ploidy
}

## setEBV: copy the ebv of an individual from the candidate pool into a stage.
setEBV <- function(pop, candidatePop) {
  pop@ebv <- matrix(
    candidatePop@ebv[match(pop@id, candidatePop@id)],
    ncol = 1
  )
  pop
}

## seedlingShare: how much of a mating plan comes from unphenotyped seedlings (S2).
## matings needs columns parent1, parent2; progeny is the number (or weight) of
## offspring allocated to each mating.
seedlingShare <- function(parentIds, matings, s2ids, progeny) {
  w <- progeny / sum(progeny)
  p1 <- matings$parent1 %in% s2ids
  p2 <- matings$parent2 %in% s2ids
  list(
    propParents = mean(parentIds %in% s2ids),
    propCrosses = mean(p1 | p2),
    propProgeny = sum(w * (p1 + p2) / 2)
  )
}

## newOutput: empty results table, one row per year.
newOutput <- function(rep, scenario, years) {
  data.frame(
    rep = as.numeric(rep),
    scenario = scenario,
    year = years,
    meanG_S1 = NA_real_,
    varG_S1 = NA_real_,
    genicVarA_S1 = NA_real_,
    accuracy_S2 = NA_real_,
    accuracy_cand = NA_real_,
    nCandidates = NA_real_,
    nParentsUsed = NA_real_,
    nCrossesUsed = NA_real_,
    nKept = NA_real_,
    Ft0 = NA_real_,
    Ft1 = NA_real_,
    dF1 = NA_real_,
    propParents_S2 = NA_real_,
    propCrosses_S2 = NA_real_,
    propProgeny_S2 = NA_real_,
    time = NA_real_
  )
}
