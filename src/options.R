###########################################################
# OPTIONS
#
# Set key options for all things model related.
#
# Any numerical, logical, or character value or vector defined
# in the main set_options() function can be overridden through
# the non-mandatory \config\my_options.csv file. Note that
# this my_options file is git ignored (indeed, that is the very
# point of such a file), so each user will need to create one.
#
###########################################################

# ---------------------------------------------------------
# Set model options and assumptions
# ---------------------------------------------------------
set_options = function(do_step = NA, quiet = FALSE) {
  
  if (!quiet) message("* Setting options")
  
  # Reset R's most annoying default options
  options(stringsAsFactors = FALSE, scipen = 999, dplyr.summarise.inform = FALSE)
  
  # Initiate options list
  o = list(do_step = do_step)
  
  # Detect user - required for checking cluster jobs
  o$user = Sys.info()[["user"]]
  
  # Name of analysis to run (cannot contain period symbol)
  o$analysis_name = "demo"
  
  # Set analysis name and create output directory system
  o = set_dirs(o)  # See directories.R
  
  # ---- Data references ----
  
  # API endpoint for national-level Oxford Stringency Index data
  o$osi_api = "https://covidtrackerapi.bsg.ox.ac.uk/api/v2/stringency/date-range/"
  
  # ---- Calibration settings ----
  
  # Take R_eff as the mean across these days
  o$fit_days = 2 : 4
  
  # Number of ASD iterations
  o$fit_iters_max = 50
  
  # Number of ASD iterations between progress plots
  o$fit_iters_plot = 5
  
  # Selection of model parameters that can be changed with the need for re-fitting
  o$fit_changeable_items = c("n_days", 
                             "model_metrics",
                             "contacts_scaler",
                             "variants_novel", 
                             "vaccine_rollout", 
                             "diagnosis_delay",
                             "testing", 
                             "isolation")
  
  # ---- Uncertainty settings ----
  
  # Which parameter set to use for 'best estimate' projection
  #
  # OPTIONS:
  #  "median" := Median of uncertainty simulations (stochastic and parameter uncertainty)
  #    "mean" := Mean of uncertainty simulations (stochastic and parameter uncertainty)
  o$best_estimate_simulation = "mean"
  
  # Number of parameters sets to sample from MCMC posteriors
  o$n_parameter_samples = 0  # Not including best parameter set
  
  # Number of seeds to run for each scenario (including baseline)
  o$n_seeds_analysis = 10
  
  # Flag to impute vaules for scenarios that did not run (mean across all other seeds)
  o$impute_failed_jobs = TRUE
  
  # Quantiles for credibility intervals
  o$quantiles = c(0.025, 0.975)
  
  # Force quantile summary even if no new simulations have been run
  o$force_summarise = FALSE
  
  # ---- Cluster settings ----
  
  # Flag for overwriting any existing simulations
  o$overwrite_simulations = FALSE
  
  # Check YAML file consistency for any existing simulations before skipping
  o$check_yaml_consistency = TRUE
  
  # Choose cluster partition for parallel jobs
  o$cluster_partition = "scicore" # OPTIONS: "covid19" or "scicore"
  
  # If running on scicore partition, you also need to define a job time 
  # 
  # NOTES:
  #  1) This should be an OVER-ESTIMATE of time needed per job - this is important as jobs over this will fail
  #  2) For more details, see: https://wiki.biozentrum.unibas.ch/display/scicore/4.+Queues+and+partitions
  o$job_time  = "00:30:00"  # Use "HH:MM:SS" or "D-HH:MM:SS" format
  o$job_queue = "30min"  # Queue to use (check out scicore wiki if unfamiliar)
  
  # Memory to allocate for each cluster job
  #
  # NOTES: 
  #  1) General principle - the higher this is set, the less resources you get
  #  2) Anything upto (and incuding) 3.7GB will give maximum resources (~500 cores at a time)
  #  3) For populations sizes of over 1m, you'll likely want something around 8GB
  o$job_memory = "8GB"
  
  # Set an upper limit for jobs that can be run at any one time
  o$job_limit = 600
  
  # Define names for cluster log and error files
  o$log_file = "scicore_log.txt"
  o$err_file = "scicore_error.txt"
  
  # Flag to remove cluster logs when jobs have successfully completed
  o$rm_cluster_log = TRUE  # Set to FALSE to bebug any cluster errors
  
  # Action to take if user is already running cluster jobs
  o$cluster_conflict_action = "error"  # Set to 'none' to turn off
  
  # ---- Plotting settings ----
  
  # Lower bound of age groups for plotting - bounded above by maximum age
  o$plot_ages = c(0, 18, 50, 65, 75, 85)
  
  # Plot a maximum number of scenarios on temporal plots
  o$max_scenarios = 25
  
  # Colour packages and palettes for scenarios, metrics, and cantons (see colour_scheme in myRfunctions.R)
  o$palette_scenario = "pals::cols25" # "brewer::dark2"
  o$palette_metric   = "base::rainbow"
  
  # Colour packages and palettes for groupings (see colour_scheme in myRfunctions.R)
  o$palette_age     = "brewer::set2"
  o$palette_variant = "brewer::set1"
  o$palette_vaccine_group = "brewer::set1"
  
  # Define some nice properties for baseline metric plots
  o$baseline_name   = "Baseline scenario"
  o$baseline_colour = "grey50"  # Light grey
  
  # Grey colour for current date dashed line
  o$data_colour = "#555555"  # Dark grey
  o$dash_colour = "#808080"  # Even darker grey

  # Saved figure size
  o$save_width  = 14
  o$save_height = 10
  
  # Image format for saving figure
  # 
  # NOTE: Use a character vector to save with multiple formats at once
  o$figure_format = "png" # Classic options: "png", "pdf", or "svg"
  
  # ---- Plotting flags ----
  
  # Turn figures on or off
  o$plot_baseline    = TRUE  # Standard baseline figures
  o$plot_cumulative  = TRUE  # Plot cumulative outcomes
  o$plot_scenarios   = TRUE  # Plot alternative (non-array) scenarios
  o$plot_arrays      = TRUE  # Plot array scenario bundles
  o$plot_heatmaps    = TRUE  # Plot heat maps for multidimension arrays
  o$plot_assumptions = TRUE  # Model structure and assumptions figures
  
  # Flags for custom figures
  o$plot_manuscript  = TRUE  # Plot figures for Omicron manuscript (December 2021)
  o$plot_custom      = TRUE  # Run my_results.R (if it exists)
  
  # ---- Advanced functionality ----
  
  # Override options set in my_options file
  o = override_options(o, quiet = quiet)
  
  # Display analysis details
  if (!quiet) message(" - Analysis name: ", o$analysis_name)
  
  return(o)
}

# ---------------------------------------------------------
# Override options set in my_options file
# ---------------------------------------------------------
override_options = function(o, quiet = FALSE) {
  
  # If user has a 'my options' file, load it
  if (file.exists(o$pth$my_options)) {
    my_options = read.csv(o$pth$my_options)
    
    # Continue if we have options to override
    if (nrow(my_options) > 0) {
      
      if (!quiet) message(" - Overriding options using config/my_options.csv file")
      
      # Throw an error if there are entries that are not well defined options
      unrecognised = my_options$option[!my_options$option %in% names(o)]
      if (length(unrecognised) > 0)
        stop("Unrecognised entries in 'my options' file: ", paste(unrecognised, collapse = ", "))
      
      # Throw an error if there are multiple entries for any one option
      duplicates = my_options$option[duplicated(my_options$option)]
      if (length(duplicates) > 0)
        stop("Duplicate entries in 'my options' file: ", paste(duplicates, collapse = ", "))
      
      # Variable class of the options we wish to overwrite
      class_conversion = paste0("as.", lapply(o[my_options$option], class))
      
      # Iterate through the options and overwrite with value of correct class
      for (i in seq_len(nrow(my_options))) {
        
        # Format entry - primary job is to seperate out multiple values
        format_value = base::trimws(str_split(my_options$value[i], ",")[[1]])
        
        # Overwrite the converted option as defined in 'my options' file
        o[[my_options$option[i]]] = get(class_conversion[i])(format_value)
      }
    }
  }
  
  return(o)
}

