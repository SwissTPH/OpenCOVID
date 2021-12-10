###########################################################
# SET DIRECTORIES
#
# Set and get directories in one place in the name of consistency and ease.
# Creates any directories that do not currently exist.
#
# OUTPUTS:
#	- A list of relevant directories (within o$pth) which can be referenced
#   elsewhere.
#
###########################################################

# ---------------------------------------------------------
# Define paths for project inputs and outputs
# ---------------------------------------------------------
set_dirs = function(o, force_analysis_name = NULL) {
  
  # Initiate file path lists
  pth = out = list()
  
  # We've already moved to code directory
  pth$code = getwd()
  
  # Path to cluster log files
  out$log = file.path(pth$code, "log")
  
  # ---- Input files ----
  
  # Parent path of all input files
  pth$input  = file.path(pth$code, "input")
  pth$config = file.path(pth$code, "config")
  pth$data   = file.path(pth$code, "data")
  
  # Paths to specific configuration files
  pth$states     = file.path(pth$config, "state_flows.xlsx")
  pth$metrics    = file.path(pth$config, "model_metrics.xlsx")
  pth$my_options = file.path(pth$config, "my_options.csv")
  
  # ---- Sub folder location (analysis name) ----
  
  # If user has a 'my options' file, load it
  if (file.exists(pth$my_options)) {
    
    # Check whether analysis name has been defined
    overwrite = filter(read.csv(pth$my_options), option == "analysis_name")
    
    # If so, overwrite value defined in options.R
    if (length(overwrite$value) > 0)
      o$analysis_name = overwrite$value
  }
  
  # Force overwrite this if desired
  if (!is.null(force_analysis_name))
    o$analysis_name = force_analysis_name
  
  # ---- Output directories ----
  
  # Parent path of all output files
  pth_output = file.path(pth$code, "output")
  
  # Path to test run files
  out$testing = file.path(pth_output, "0_testing")
  
  # Path to calibration files
  out$fitting = file.path(pth_output, "1_calibration", o$analysis_name)
  
  # Paths to scenario files
  out$scenarios   = file.path(pth_output, "2_scenarios", o$analysis_name)
  out$simulations = file.path(out$scenarios, "simulations")
  out$arrays      = file.path(out$scenarios, "array_info")
  
  # Path to figures and other outputs
  out$figures = file.path(pth_output, "3_figures", o$analysis_name)
  
  # ---- Create directory structure ----
  
  # Make all output directories
  make_out_dirs(out)
  
  # Append paths to o list
  o = append_dirs(o, pth, out)
  
  return(o)
}

# ---------------------------------------------------------
# Make all output directories if they do not already exist
# ---------------------------------------------------------
make_out_dirs = function(out) {
  
  # Extract all path names in list
  pth_names = names(out)
  
  # Loop through these path names
  for (pth_name in pth_names) {
    this_pth = out[[pth_name]]
    
    # If it does not already exist, create it
    if (!dir.exists(this_pth) & !grepl("\\*", this_pth))
      dir.create(this_pth, recursive = TRUE)
  }
}

# ---------------------------------------------------------
# Concatenate separators and append directories to o list
# ---------------------------------------------------------
append_dirs = function(o, pth, out) {
  
  # Extract all path names in list
  pth_names = names(out)
  
  # Loop through these path names
  for (pth_name in pth_names) {
    this_pth = out[[pth_name]]
    
    # We use * to denote a partial path
    if (grepl("\\*", this_pth)) {
      out[[pth_name]] = substr(this_pth, 1, nchar(this_pth) - 1)
      
    } else {  # Otherwise add a file separator to end of output paths
      out[[pth_name]] = paste0(this_pth, .Platform$file.sep)
    }
  }
  
  # Concatenate lists
  o$pth = c(pth, out)
  
  return(o)
}

