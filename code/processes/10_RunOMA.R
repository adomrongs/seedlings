# Optimal mate allocation over the whole genotyped pool, at a fixed inbreeding
# rate dF. The mating plan is then turned into an actual F1 with sim_mate().

object <- COMA::read_data(
  geno.file = geno.file,
  kinship.file = k.file,
  ploidy = ploidy,
  matings = 'all',
  standardize = FALSE
)

parents_df <- data.frame(object$parents['id'], min = 0, max = 1)
matings_df <- data.frame(object$matings, min = 0, max = 1)

# read_data() already drops reciprocals: for hermaphrodites it keeps only
# female >= male, so (A,B) and (B,A) never coexist. What it does keep are the
# selfs, and the filter below is what removes them. They have to go: the trait
# is purely additive, so a self carries no inbreeding depression in its merit
# and only the kinship constraint argues against it. Leaving them in would also
# inflate the seedling share, since a self of an S2 counts as two S2 slots.
matings_df <- matings_df[matings_df$parent1 != matings_df$parent2, ]

# ----------- Run OMA -----------
start <- Sys.time()
oma_object <- COMA::oma(
  dF = dF,
  parents = parents_df,
  matings = matings_df,
  ploidy = ploidy,
  K = object$K,
  base = 'current',
  dF.adapt = list(step = 0.005, max = 0.1)
)
time <- as.numeric(difftime(Sys.time(), start, units = 'mins'))

# ----------- Share of the mating plan coming from seedlings -----------
plan <- oma_object$om[oma_object$om$value > 0, ]
share <- seedlingShare(
  parentIds = unique(c(plan$parent1, plan$parent2)),
  matings = plan,
  s2ids = S2@id,
  progeny = plan$value
)

# ----------- Turn the plan into progeny -----------
oma_plan <- SimPlus::sim_mate(
  pop = candidatePop,
  SP = SP,
  matings = oma_object$om,
  total.progeny = totalProgeny,
  min.progeny = minProgeny
)

old <- parents
parents <- candidatePop[
  candidatePop@id %in%
    unique(c(oma_plan$matings$parent1, oma_plan$matings$parent2))
]
nKept <- sum(old@id %in% parents@id)
nParentsUsed <- parents@nInd
nCrossesUsed <- nrow(oma_plan$matings)
F1_tmp <- oma_plan$progeny

cat(
  "    OMA done in",
  round(time, 1),
  "min |",
  nCrossesUsed,
  "crosses |",
  round(100 * share$propProgeny, 1),
  "% seedlings\n"
)
