###########################################################
# DEPENDENCIES
#
# Deal with all package dependencies in one place.
#
###########################################################

# ---- R version check ----

# R versions for which this project has been tested and is stable
stable_versions = "4.1.0"

# R versions for which this project is stable (as a string)
stable_str = paste(stable_versions, collapse = ", ")

# Get details of R version currently running
version_info = R.Version()

# Construct version number from list details
version_num = paste0(version_info$major, ".",  version_info$minor)

# Throw an error if this R version is unsuitable
if (!version_num %in% stable_versions)
  stop("This software is stable with R version(s): ", stable_str, 
       " (currently running ", version_num, ")")

# ---- Source all files ----

# NOTE: This isn't necessarily all R files in the directory

source("auxiliary.R")
source("options.R")
source("directories.R")
source("parse_input.R")
source("model.R")
source("networks.R")
source("calibration.R")
source("load_data.R")
source("emulator.R")
source("asd.R")
source("scenarios.R")
source("cluster_jobs.R")
source("postprocess.R")
source("results.R")
source("plotting.R")
source("unit_tests.R")

# Also load my_results if it exists
if (file.exists("my_results.R"))
  source("my_results.R")

# ---- Define packages ----

# Complete list of all R packages required for this project
packages = c("tidyverse",      # Includes ggplot2, dplyr, tidyr, and others (https://www.tidyverse.org/packages/)
             "tgp",            # Latin hypercube sampler
             "hetGP",          # Gaussian Process model and acquisition functions
             "philentropy",    # Distance measures
             "smooth",         # Simple moving average model
             "splines",        # Spline fitting models
             "forecast",       # Linear regression for time series
             "imputeTS",       # Imputation for time series
             "akima",          # Bivariate interpolation
             "stats",          # Statistical calculations and random number generation
             "matrixStats",    # Matrix row and column operations
             # "tidygraph",      # Network functionality
             # "ggraph",         # Network functionality
             "widyr",          # Compile network properties
             "socialmixr",     # Age structured contact matrixes from POLYMOD
             "EpiEstim",       # R_effective calculated from incidence and serial interval
             "wrswoR",         # Fast weighted integer sampling without replacement
             "data.table",     # Next generation dataframes
             "useful",         # General helper functions (eg compare.list)
             "rlist",          # List-related helper functions (eg list.remove)
             "gsubfn",         # Output multiple variables from functions
             "httr",           # Read data from API endpoint
             "jsonlite",       # Convert data to/from json format
             "rio",            # Data loading functionality
             "readxl",         # Data loading functionality
             "yaml",           # Data loading functionality
             "lubridate",      # Data formatting functionality
             "naniar",         # Data formatting functionality
             "wrapr",          # Convenience functions (eg qc)
             "coda",           # Plotting functionality
             "gridExtra",      # Plotting functionality
             "ggnewscale",     # Plotting functionality
             "ggpubr",         # Plotting functionality
             "cowplot",        # Plotting functionality
             "scales",         # Plotting functionality
             "ggh4x",          # Plotting functionality (flexible faceting)
             "ggtext",         # Plotting functionality (use markdown in labels)
             "pals",           # Colour palettes
             "colorspace",     # Colour palettes
             "GGally",         # Network plotting
             "igraph",         # Network plotting
             "RColorBrewer",   # Network plotting colour palettes
             "progress")       # Cool progress bar

# List of all packages only available from github
gh_packages = c("eliocamp/tagger")

# ---- Install and/or load packages with pacman ----

message("* Installing required packages")

# Check whether pacman itself has been installed
pacman_installed = "pacman" %in% rownames(installed.packages())

# If not, install it
if (!pacman_installed) 
  install.packages("pacman")

# Load pacman
library(pacman) 

# Load all required packages, installing them if required
pacman::p_load(char = packages)

# Same for github packages
pacman::p_load_gh(gh_packages)

# ---- Redefine or unmask particular functions ----

# Unmask certain functions otherwise overwritten
union   = dplyr::union
select  = dplyr::select
filter  = dplyr::filter
rename  = dplyr::rename
recode  = dplyr::recode
predict = stats::predict
# groups  = tidygraph::groups
as.data.frame = base::as.data.frame

