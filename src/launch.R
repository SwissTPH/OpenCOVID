###########################################################
# LAUNCH
#
# Main launch function for a dynamic COVID-19 individual-based 
# transmission model developed at Swiss TPH.
#
# Easiest use is to 'source' this file without echo.
# Alternativly, call 'sh bash_launch.sh' from shell command line.
# 
# Stable with R3.6.0.
#
###########################################################
  
# Clear global environment
rm(list = ls())

# Set working directory to sourced file
if (interactive()) setwd(getSrcDirectory(function() {}))

source("dependencies.R")

# Tidy u
if (interactive()) clf()  # Close figures 
if (interactive()) clc()  # Clear console

message("Running COVID-19 IBM transmission model")

# Set options (see options.R)
o = set_options(do_step = 1)

# Step NA) Test run a single simulation
run_model_test(o)  # See unit_tests.R

# Step 0) Visualise cantonal data
run_data_exploration(o)  # See load_data.R

# Step 1) Train a model emulator
run_emulator(o)  # See emulator.R

# Step 2) Calibrate model
run_calibration(o)  # See calibration.R

# Step 3) Run analysis scenarios and strategies
run_analysis(o)  # See analysis.R

# Step 4) Plot results
run_results(o)  # See results.R

# Finish up
message("* Finished!")

