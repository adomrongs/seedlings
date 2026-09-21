# One-off environment setup. Run it once, from the project root, with the
# seedlings.Rproj already open. It is not sourced by any scenario.
#
# Why bare = TRUE: renv::init() by default scans the code, finds every
# library() call and tries to install all of them from CRAN. Half of this
# project lives on GitHub, that scan fails, and renv leaves the project in a
# half-initialised state. Starting bare and installing explicitly avoids it.

# ----------- 1. Empty renv environment -----------
install.packages('renv')
renv::init(bare = TRUE)

# ----------- 2. CRAN packages -----------
renv::install(c(
  'AlphaSimR',
  'dplyr',
  'data.table',
  'AGHmatrix',
  'lme4'
))

# ----------- 3. GitHub packages -----------
# Installing them through renv::install() with the user/repo form is what makes
# renv record the remote in renv.lock, so renv::restore() can rebuild the same
# library on the cluster.
renv::install('covaruber/lme4breeding')
renv::install('jendelman/polyBreedR')
renv::install('jendelman/COMA')
# TODO: SimPlus is only used by 10_RunOMA.R (sim_mate). Replace the line below
# with the correct remote before running the COMA scenario.
# renv::install('<user>/SimPlus')

# ----------- 4. Freeze -----------
renv::snapshot()

# ----------- 5. On the cluster -----------
# Clone the repo, open R at the project root and run:
#   renv::restore()
# If GitHub rate-limits the build, set a token first:
#   Sys.setenv(GITHUB_PAT = '...')
# Sharing one cache across jobs saves a lot of time; put this in ~/.Renviron:
#   RENV_PATHS_CACHE=/path/to/shared/renv/cache
