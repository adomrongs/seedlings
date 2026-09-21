# PS control: the burnin program carried on unchanged for nFuture more years.
# Phenotypic selection at every stage, crossing block recycled from S4.
# This is the baseline the two GS scenarios are compared against.

rm(list = ls())
library(AlphaSimR)
library(dplyr)
library(polyBreedR)
source('code/processes/0_params.R')

# ----------- Replicate and outputs -----------
rep <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"))
load(paste0('outputs/burnin/rdata/burnin_', rep, '.Rdata'))
source('code/processes/0_params.R')
scenario <- 'PS'

csv <- paste0('outputs/', scenario, '/csv/', scenario, '_', rep, '.csv')
if (!dir.exists(dirname(csv))) dir.create(dirname(csv), recursive = TRUE)

output_future <- newOutput(rep, scenario, (nBurnin + 1):(nBurnin + nFuture))

# ----------- Run -----------
cat("Running", scenario, "rep", rep, "\n")

for (year in (nBurnin + 1):(nBurnin + nFuture)) {
  cat("  year", year, "\n")
  i <- year - nBurnin

  source('code/processes/3_UpdateParents.R')
  source('code/processes/4_AdvanceYear.R')
  source('code/processes/5_StoreRecords.R')
  source('code/processes/6_Inbreeding.R')

  output_future$meanG_S1[i] <- meanG(S1)
  output_future$varG_S1[i] <- varG(S1)
  output_future$genicVarA_S1[i] <- genicVarA(S1)
  output_future$accuracy_S2[i] <- accuracy_S2
  output_future$nParentsUsed[i] <- parents@nInd
  output_future$nCrossesUsed[i] <- nCrosses
  output_future$nKept[i] <- nKept
  output_future$Ft0[i] <- Ft0
  output_future$Ft1[i] <- Ft1
  output_future$dF1[i] <- dF1
}

# ----------- Save -----------
final_output <- bind_rows(output, output_future)
final_output$scenario <- scenario
write.csv(final_output, file = csv, row.names = FALSE)
cat(scenario, "done\n")
