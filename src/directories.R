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
set_dirs = function(o) {
  
  # Initiate file path lists
  pth = out = list()
  
  # We've already moved to code directory
  pth$code = getwd()
  
  # Path to cluster log files and data cache
  out$log   = file.path(pth$code, "log")
  out$cache = file.path(pth$code, "cache")
  
  # ---- Input and configuration files ----
  
  # Parent path of all input files
  pth$input  = file.path(pth$code, "input")
  pth$config = file.path(pth$code, "config")
  
  # Paths to model configuration files
  pth$states    = file.path(pth$config, "model_states.yaml")
  pth$metrics   = file.path(pth$config, "model_metrics.yaml")
  pth$variables = file.path(pth$config, "model_variables.yaml")
  
  # ---- Set analysis name ----
  
  # User's custom options (may or may not exist)
  pth$my_options = file.path(pth$config, "my_options.yaml")
  
  # If user has a 'my options' file, load it
  if (file.exists(pth$my_options)) {
    
    # Check whether analysis name has been defined
    overwrite_name = read_yaml(pth$my_options)$analysis_name
    
    # If so, overwrite value defined in options.R
    if (!is.null(overwrite_name))
      o$analysis_name = overwrite_name
  }
  
  # ---- Model parameter files ----
  
  # Default model parameters
  pth$params_default = file.path(pth$config, "default.yaml")
  
  # User-specified model parameters for this analysis
  pth$params_user = file.path(pth$input, paste0(o$analysis_name, ".yaml"))
  
  # ---- Output directories ----
  
  # Parent path of all output files
  pth_output = file.path(pth$code, "output")
  
  # Path to test run files
  out$testing = file.path(pth_output, "0_testing")
  
  # Path to calibration files
  out$fitting     = file.path(pth_output, "1_calibration", o$analysis_name)
  out$fit_samples = file.path(out$fitting, "fit_samples")
  
  # Paths to scenario files
  pth_scenarios   = file.path(pth_output, "2_scenarios", o$analysis_name)
  out$scenarios   = file.path(pth_scenarios, "scenarios")
  out$simulations = file.path(pth_scenarios, "simulations")
  out$uncertainty = file.path(pth_scenarios, "uncertainty")
  
  # Path to predictor models for LHC scenarios
  pth_arrays     = file.path(pth_output, "3_arrays", o$analysis_name)
  out$array_info = file.path(pth_arrays, "array_info")
  out$parents    = file.path(pth_arrays, "grid_parents")  # TODO: Pre-summarise grid parents
  out$endpoints  = file.path(pth_arrays, "lhc_endpoints")
  
  # Path to figures and other output results
  out$results = file.path(pth_output, "4_results", o$analysis_name)
  
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

