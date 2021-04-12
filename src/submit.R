############################################################
# SUBMIT
#
# Submit model simulations from a cluster array.
#
############################################################

source("dependencies.R")

# Extract input from bash file
args = commandArgs(trailingOnly = TRUE)

# Name arguments when provided from bash
if (length(args) > 0) {
  sim_type = as.character(args[1])
  task_id  = as.numeric(args[2])
  
} else {  # Otherwise define some defaults - used for testing and debugging
  sim_type = "calibration::0"
  task_id  = 1
}

# Call main function and catch any errors
tryCatch(
  expr = do_simulate(sim_type, task_id),  # See simulate.R
  
  # Error handler
  error = function(e) {
    
    # Concatenate (ideally useful) error message
    err_message = paste0(" ! ", e$message, " (array ID: ", task_id, ")")
    
    # Append this error message to an error log file
    write(err_message, file = "scicore_error.txt", append = TRUE)
  },
  
  # Close up
  finally = message("Closing cluster job")
)

