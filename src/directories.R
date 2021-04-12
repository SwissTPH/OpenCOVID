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
  
  # Parent path of all input files
  pth$config = file.path(pth$code, "config")
  pth$data   = file.path(pth$code, "data")
  
  # ---- Input files ----
  
  # Paths to specific configuration files
  pth$states     = file.path(pth$config, "state_flows.xlsx")
  pth$metrics    = file.path(pth$config, "model_metrics.xlsx")
  pth$calib_file = file.path(pth$config, "model_parameters.xlsx")
  pth$calib_mult = file.path(pth$config, "calibration_weights.csv")
  pth$response   = file.path(pth$config, "response_efficacy.csv")
  pth$cantons    = file.path(pth$config, "canton_codes.csv")
  pth$stations   = file.path(pth$config, "weather_stations.txt")
  pth$my_options = file.path(pth$config, "my_options.csv")
  
  # Paths to specific data files
  pth$capacity = file.path(pth$data, "data_icu_capacity.csv")
  pth$seroprev = file.path(pth$data, "data_seroprevalence.csv")
  pth$risk_pop = file.path(pth$data, "data_risk_groups.csv")
  pth$disease  = file.path(pth$data, "disease_probabilities.csv")
  pth$variants = file.path(pth$data, "variant_properties.csv")
  pth$vaccines = file.path(pth$data, "vaccine_properties.csv")
  
  # ---- Analysis and calibration names ---
  
  # If user has a 'my options' file, load it
  if (file.exists(pth$my_options)) {
    
    # Check whether analysis or calibration names have been defined
    overwrite_names = read.csv(pth$my_options) %>%
      filter(option %in% qc(analysis_name, calibration_name))
    
    # If so, overwrite any names as defined in options.R
    for (i in seq_len(nrow(overwrite_names)))
      o[[overwrite_names[i, ]$option]] = overwrite_names[i, ]$value
  }
  
  # If calibration_name is yet to be defined, set it as analysis name
  if (is.na(o$calibration_name))
    o$calibration_name = o$analysis_name
  
  # Also store full set of canton names - useful when loading data
  o$all_switzerland = read.csv(pth$cantons)$code    # With national level
  o$all_cantons = setdiff(o$all_switzerland, "CH")  # Without national level
  
  # ---- Output directories ----
  
  # Parent path of all output files
  pth_output = file.path(pth$code, "output")
  
  # Paths to emulator directories
  out$testing = file.path(pth_output, "0_testing")
  
  # Paths to emulator directories
  out$emulator = file.path(pth_output, "1_emulator", o$calibration_name)
  out$samples  = file.path(out$emulator, "samples*")
  
  # Paths to calibration directories
  out$fitting  = file.path(pth_output, "2_calibration", o$calibration_name)
  out$chains   = file.path(out$fitting, "chains")
  out$archive  = file.path(out$fitting, "archive")
  
  # Paths to analysis directories
  out$analysis = file.path(pth_output, "3_analyses", o$analysis_name)
  out$sims     = file.path(out$analysis, "simulations")
  out$scenario = file.path(out$analysis, "scenarios")
  out$strategy = file.path(out$analysis, "strategies")
  
  # Path to figures (all outputs now directed to this one folder)
  out$figures  = file.path(pth_output, "4_figures", o$analysis_name)
  
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

