#!/bin/bash

#SBATCH --job-name=simulate
#SBATCH --account=penny
#SBATCH --cpus-per-task=1
#SBATCH --output=scicore_output.txt

############################################################
# BASH SUBMIT
#
# Submit array of COVID-19 IBM simulations.
#
# Arguments:
#  SIM_TYPE: Type of job we're simulating - see simulate.R for options
#  LOG_FILE: Text log file that stores names of completed scenarios
#
# Various sbatch options can be set in options.R (see 'cluster settings' 
# section). For example, cluster partition, job memory, and slurm queue.
#
# Written by A.J.Shattock
############################################################

# Load R
module purge
ml R/3.6.0-foss-2018b

# Extract inputs
sim_type=$1
log_file=$2

# Extract array task ID
task_id=$(expr ${SLURM_ARRAY_TASK_ID})

echo "Simulating model for parameter set $task_id"

# Run R script which calls simulation function
Rscript submit.R $sim_type $task_id

# Write scenario name to log file
echo "$task_id" >> $log_file
