# Burnin: 20 years of the breeding program under phenotypic selection only.
# The final state is saved so that every future scenario starts from the same
# point for a given replicate.

rm(list = ls())
library(AlphaSimR)
library(dplyr)
library(polyBreedR)
source('code/processes/0_params.R')

# ----------- Replicate and outputs -----------
rep <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"))
scenario <- 'burnin'
set.seed(rep)

image <- paste0('outputs/', scenario, '/rdata/', scenario, '_', rep, '.Rdata')
if (!dir.exists(dirname(image))) dir.create(dirname(image), recursive = TRUE)

output <- newOutput(rep, scenario, 1:nBurnin)

# ----------- Run -----------
cat("Running", scenario, "rep", rep, "\n")

source('code/processes/1_CreateFounders.R')
source('code/processes/2_FillPipeline.R')

for (year in 1:nBurnin) {
  cat("  year", year, "\n")

  source('code/processes/3_UpdateParents.R')
  source('code/processes/4_AdvanceYear.R')
  source('code/processes/5_StoreRecords.R')
  source('code/processes/6_Inbreeding.R')

  output$meanG_S1[year] <- meanG(S1)
  output$varG_S1[year] <- varG(S1)
  output$genicVarA_S1[year] <- genicVarA(S1)
  output$accuracy_S2[year] <- accuracy_S2
  output$nParentsUsed[year] <- parents@nInd
  output$nCrossesUsed[year] <- nCrosses
  output$nKept[year] <- nKept
  output$Ft0[year] <- Ft0
  output$Ft1[year] <- Ft1
  output$dF1[year] <- dF1
}

# ----------- Save -----------
save.image(file = image)
cat("Burnin done\n")
