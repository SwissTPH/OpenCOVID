###########################################################
# LAUNCH
#
# Main launch function for OpenCOVID, a dynamic individual-based 
# model of SARS-CoV-2 transmission and COVID-19 disease.
#
# Usage:
#  Interactively via Rstudio: 'source' this file (without echo)
#  From shell command line: sh bash_launch.sh
#
# The OpenCOVID model is maintained by members of the Disease
# Modelling Unit at Swiss TPH.
#
###########################################################

# Clear global environment
rm(list = ls())

# Set working directory to sourced file
if (interactive()) setwd(getSrcDirectory(function() {}))

source("dependencies.R")

# Tidy up
if (interactive()) clf()  # Close figures
if (interactive()) clc()  # Clear console

message("Running OpenCOVID v2.0")

# Set options (see options.R)
o = set_options(do_step = 0)

# Step 0) Test run a single simulation
run_model_test(o)  # See unit_tests.R

# Step 1) Fit effective reproduction number
run_calibration(o)  # See calibration.R

# Step 2) Run all scenarios
run_scenarios(o)  # See scenarios.R

# Step 3) Plot results
run_results(o)  # See results.R

# Finish up
message("* Finished!")

