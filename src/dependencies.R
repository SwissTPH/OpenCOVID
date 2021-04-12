###########################################################
# DEPENDENCIES
#
# Deal with all package dependencies in one place.
#
###########################################################

# ---- R version check ----

# NOTE: Several packages seem to have back-compatibility issues past 3.6.0, so stick with this to be on the safe side

# Get details of R version currently running
version_info = R.Version()

# Construct version number from list details
version_num = paste0(version_info$major, ".",  version_info$minor)

# Throw an error if this version is anything other than 3.6.0
if (version_num != "3.6.0")
  stop("This software is stable with R version 3.6.0 (currently running ", version_info$version.string, ")")

# ---- Source all files ----

# NOTE: This isn't necessarily all R files in the directory

source("myRfunctions.R")
source("options.R")
source("directories.R")
source("load_input.R")
source("load_data.R")
source("parameters.R")
source("model.R")
source("emulator.R")
source("likelihood.R")
source("calibration.R")
source("analysis.R")
source("scenarios.R")
source("strategies.R")
source("simulate.R")
source("postprocess.R")
source("results.R")
source("plotting.R")
source("unit_tests.R")

# ---- Define packages ----

# Complete list of all R packages required for this project
packages = c("tidyverse",      # Includes ggplot2, dplyr, tidyr, and others (https://www.tidyverse.org/packages/)
             "tgp",            # Latin hypercube sampler
             "hetGP",          # Gaussian Process model and acquisition functions
             "philentropy",    # Distance measures
             "smooth",         # Simple moving average model
             "RMAWGEN",        # Monthly to daily weather interpolation
             "forecast",       # Linear regression for time series
             "imputeTS",       # Imputation for time series
             "stats",          # Statistical calculations and random number generation
             "matrixStats",    # Matrix row and column operations
             "tidygraph",      # Network functionality
             "ggraph",         # Network functionality
             "socialmixr",     # Age structured contact matrixes from POLYMOD
             "widyr",          # Compiling network properties
             # "Rfast",        # Numerous C++ speed functions [RD having compiler complications]
             "wrswoR",         # Fast weighted integer sampling without replacement
             "data.table",     # Next generation dataframes
             # "RPostgreSQL",  # SQL database connections [AS having installation issues]
             "rio",            # Data loading functionality
             "readxl",         # Data loading functionality
             "lubridate",      # Data formatting functionality
             "naniar",         # Data formatting functionality
             "wrapr",          # Convenience functions (eg qc)
             # "operators",    # Convenience operators (eg %+=%) [Probably overkill]
             "coda",           # Plotting functionality
             "gridExtra",      # Plotting functionality
             "ggpubr",         # Plotting functionality
             "cowplot",        # Plotting functionality
             "scales",         # Plotting functionality
             "pals")           # Colour palettes

# List of all packages available from github
gh_packages = c("jameshay218/lazymcmc")  # MCMC algorithm

# ---- Install and/or load packages with pacman ----

message("* Installing required packages")

# Check whether pacman itself has been installed
pacman_installed = "pacman" %in% rownames(installed.packages())

# If not, install it
if (!pacman_installed) 
  install.packages("pacman")

# Load pacman
require(pacman) 

# Load all required packages, installing them if required
pacman::p_load(char = packages)

# Same for github packages
pacman::p_load_gh(gh_packages)

# ---- Redefine or unmask particular functions ----

# Redefine Rfast functions as throwing c++ compilier errors
colsums = base::colSums
rowsums = base::rowSums

# Unmask certain functions otherwise overwritten
union  = dplyr::union
select = dplyr::select
filter = dplyr::filter
groups = tidygraph::groups
as.data.frame = base::as.data.frame

