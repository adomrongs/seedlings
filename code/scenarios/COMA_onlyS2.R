# COMA con la candidate pool formada por solo seedlings. Sin padres y sin clones evaluados: no hay solapamiento
# de generaciones y cada ano el crossing block se renueva por completo.
# Es el control 'unicamente seedlings' del PDF.
# Todo lo demas (burnin, training population, GS de S2 a S3, dF) es identico
# entre los cinco escenarios COMA.

rm(list = ls())
library(AlphaSimR)
library(dplyr)
library(data.table)
library(AGHmatrix)
library(lme4breeding)
library(polyBreedR)
library(COMA)
library(SimPlus)
source('code/processes/0_params.R')

# ----------- Replicate and outputs -----------
rep <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"))
load(paste0('outputs/burnin/rdata/burnin_', rep, '.Rdata'))
source('code/processes/0_params.R')
scenario <- 'COMA_onlyS2'
poolStages <- c('S2')

csv <- paste0('outputs/', scenario, '/csv/', scenario, '_', rep, '.csv')
if (!dir.exists(dirname(csv))) dir.create(dirname(csv), recursive = TRUE)

output_future <- newOutput(rep, scenario, (nBurnin + 1):(nBurnin + nFuture))

# ----------- Run -----------
cat("Running", scenario, "rep", rep, "\n")

for (year in (nBurnin + 1):(nBurnin + nFuture)) {
  cat("  year", year, "\n")
  i <- year - nBurnin

  source('code/processes/7_PredictGEBV.R')
  source('code/processes/9_COMAFiles.R')
  source('code/processes/10_RunOMA.R')
  source('code/processes/11_AdvanceYearGS.R')
  source('code/processes/5_StoreRecords.R')
  source('code/processes/6_Inbreeding.R')

  output_future$meanG_S1[i] <- meanG(S1)
  output_future$varG_S1[i] <- varG(S1)
  output_future$genicVarA_S1[i] <- genicVarA(S1)
  output_future$accuracy_S2[i] <- accuracy_S2
  output_future$accuracy_cand[i] <- accuracy_cand
  output_future$nCandidates[i] <- nCandidates
  output_future$nParentsUsed[i] <- nParentsUsed
  output_future$nCrossesUsed[i] <- nCrossesUsed
  output_future$nKept[i] <- nKept
  output_future$Ft0[i] <- Ft0
  output_future$Ft1[i] <- Ft1
  output_future$dF1[i] <- dF1
  output_future$propParents_S2[i] <- share$propParents
  output_future$propCrosses_S2[i] <- share$propCrosses
  output_future$propProgeny_S2[i] <- share$propProgeny
  output_future$time[i] <- time
}

# ----------- Save -----------
final_output <- bind_rows(output, output_future)
final_output$scenario <- scenario
write.csv(final_output, file = csv, row.names = FALSE)
cat(scenario, "done\n")
