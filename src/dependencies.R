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

# ---- Source files ----

# NOTE: This isn't necessarily all R files in the directory

# Load R files and functions
source("auxiliary.R")
source("options.R")
source("directories.R")
source("parse_input.R")
source("array.R")
source("model.R")
source("prognosis.R")
source("networks.R")
source("calibration.R")
source("load_data.R")
source("emulator.R")
source("asd.R")
source("scenarios.R")
source("uncertainty.R")
source("cluster_jobs.R")
source("postprocess.R")
source("results.R")
source("plotting.R")
source("manuscript.R")
source("unit_tests.R")

# Also load my_results if it exists
if (file.exists("my_results.R"))
  source("my_results.R")

# List of python files to source
py_source = c("auxiliary.py")

# ---- Define packages ----

# Complete list of all R packages required for this project
packages = c("tidyverse",      # Includes ggplot2, dplyr, tidyr (www.tidyverse.org/packages/)
             "magrittr",       # Additional pipe operators, such as %<>%
             "data.table",     # Next generation dataframes
             "useful",         # General helper functions (eg compare.list)
             "rlist",          # List-related helper functions (eg list.remove)
             "gsubfn",         # Output multiple variables from functions
             "wrapr",          # Convenience functions (eg qc)
             "stats",          # Statistical calculations and random number generation
             "matrixStats",    # Matrix row and column operations
             "tgp",            # Latin hypercube sampler
             "hetGP",          # Gaussian Process model and acquisition functions
             "philentropy",    # Distance measures
             "smooth",         # Simple moving average model
             "splines",        # Spline fitting models
             "forecast",       # Linear regression for time series
             "imputeTS",       # Imputation for time series
             "akima",          # Bivariate interpolation
             "widyr",          # Compile network properties
             "socialmixr",     # Age structured contact matrixes from POLYMOD
             "EpiEstim",       # Re calculated from incidence and serial interval
             "wrswoR",         # Fast weighted integer sampling without replacement
             "phonenumber",    # Letter to number conversion
             "progress",       # Progress bar
             "tictoc",         # Code timer
             "httr",           # Read data from API endpoint
             "jsonlite",       # Convert data to/from json format
             "reticulate",     # Execute python functions in R
             "rio",            # Data loading functionality
             "readxl",         # Data loading functionality
             "yaml",           # Data loading functionality
             "lubridate",      # Data formatting functionality
             "naniar",         # Data formatting functionality
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
             "RColorBrewer",   # Colour palettes
             "GGally",         # Network plotting
             "igraph")         # Network plotting

# List of all packages only available from github
gh_packages = c("eliocamp/tagger")

# List of python packages (and associated modules) to be installed
py_packages = c("pyyaml::yaml")

# ---- Install and/or load R packages with pacman ----

message("* Installing required R packages")

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

# ---- Python packages and files ----

# Specify conda environment
# use_condaenv("r-reticulate")

# Skip this step if no python packages needed
if (length(py_packages) > 0) {
  
  message("* Installing required python packages")
  
  # Split package and module names
  py_modules = str_split(py_packages, "::", simplify = TRUE)
  py_modules = setNames(py_modules[, 2], py_modules[, 1])
  
  # Loop through modules
  for (i in seq_along(py_modules)) {
    
    # This module from this package
    py_module  = py_modules[i] %>% unname()
    py_package = py_modules[i] %>% names()
    
    # Only install package if module is not available
    if (!py_module_available(py_module))
      reticulate::conda_install("r-reticulate", py_package)
    
    # Import the module
    reticulate::import(py_module)
  }
  
  # Now ready to source all python files 
  for (py_file in py_source)
    reticulate::source_python(py_file)
}

# ---- Redefine or unmask particular functions ----

# Unmask certain functions otherwise overwritten
select  = dplyr::select
filter  = dplyr::filter
rename  = dplyr::rename
recode  = dplyr::recode
count   = dplyr::count
union   = dplyr::union
predict = stats::predict

