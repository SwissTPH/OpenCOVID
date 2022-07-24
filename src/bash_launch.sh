#!/bin/bash

############################################################
# BASH LAUNCH
#
# Launch COVID-19 IBM pipeline.
#
# Command line usage:
#   sh bash_launch.sh
#
############################################################

# Load R
module purge
ml R/3.6.0-foss-2018b

# Call main launch script
Rscript launch.R

