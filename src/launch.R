###########################################################
# LAUNCH
#
# Main launch function for OpenCOVID, a dynamic individual-based 
# model of SARS-CoV-2 transmission and COVID-19 disease.
#
# Usage:
#  Via Rstudio: Open opencovid.Rproj and 'source' this file (without echo)
#  Via shell command line: sh bash_launch.sh
#
# The OpenCOVID model is maintained by members of the Disease
# Modelling Unit at the Swiss Tropical and Public Health Institute.
#
###########################################################

# Clear global environment
rm(list = ls())

# Set working directory to sourced file
if (interactive()) setwd(getSrcDirectory(function() {}))

# Load all required packages and functions
source("dependencies.R")

# Tidy up
if (interactive()) clf()  # Close figures
if (interactive()) clc()  # Clear console

message("Running OpenCOVID v4.2 alpha")

# Set options (see options.R)
o = set_options(do_step = 0)

# Step 0) Test run a single simulation
run_model_test(o)  # See unit_tests.R

# Step 1) Calibrate model
run_calibration(o)  # See calibration.R

# Step 2) Run all scenarios
run_scenarios(o)  # See scenarios.R

# Step 3) Operate on array scenarios
run_arrays(o)  # See array.R

# Step 4) Plot results
run_results(o)  # See results.R

# Finish up
message("* Finished!")

