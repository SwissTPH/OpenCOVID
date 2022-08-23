#!/bin/bash

############################################################
# BASH COPY PROJECT
#
# Upload or download OpenCOVID simulation output to/from the
# shared GROUP directory on scicore.
#
# Command line usage:
#   sh bash_copy.sh upload <analysis_name>  # For uploading files
#   sh bash_copy.sh download <analysis_name>  # For downloading files
#
############################################################

# Load R
module purge
ml R/4.1.0-foss-2018b

# Extract inputs
copy_direction=$1
analysis_name=$2

# Call copy project script with any provided arguments
Rscript copy_project.R $copy_direction $analysis_name

