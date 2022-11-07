###########################################################
# CALIBRATION
#
# Fit parameters to either 1) a set of epi metrics over time, 
# or 2) Re at the 'initial' point for forward projections.
#
# Uses adaptive sampling to improve likelihood of locating 
# global minimum of error between model output and fitting
# target.
#
###########################################################

# ---------------------------------------------------------
# Parent function for model calibration process
# ---------------------------------------------------------
run_calibration = function(o) {
  
  # Only continue if specified by do_step
  if (!is.element(1, o$do_step)) return()
  
  message("* Calibrating model")
  
  # ---- Load input and data ----
  
  # Initiate fit list and perform a few checks on input yaml
  fit = setup_calibration(o)
  
  # Load data - see load_data.R
  fit = load_data(o, fit, fit$synthetic)
  
  # ---- Main adaptive sampling loop ----
  
  # Iterate through adaptive sampling rounds (including initial 'r0' step)
  for (r_val in 0 : fit$input$adaptive_sampling$rounds) {
    
    message(" - Adaptive sampling round ", r_val)
    
    # If we want this process to be reproducible, set a seed on each round
    if (o$fit_reproducible)
      set.seed(r_val)
    
    # Sample parameter sets for this sampling round
    param_ids = sample_parameters(o, fit, r_val)
    
    # Break out of adaptive sampling if no more decent samples to simulate
    if (length(param_ids) == 0)
      break
    
    # Update round index string
    r_idx = paste0("r", r_val)
    
    # Simulate these parameter samples
    simulate_parameters(o, r_idx, param_ids)
    
    # Calculate quality of fit
    quality_of_fit(o, fit, r_idx)
    
    # Select all parameter sets that we'll use to train the emulator
    select_samples(o, r_val)
    
    # Train model emulator (see emulator.R)
    train_emulator(o, fit, r_idx)
    
    # Find global minimum of emulated space (see emulator.R)
    search_emulator(o, fit, r_idx)
    
    message("  > Plotting performance")
    
    # Performance plots (see plotting.R)
    plot_best_samples(o, "Best simulated samples",   r_idx)
    plot_emulator(o,     "Emulator performance",     r_idx)
    plot_optimisation(o, "Optimisation performance", r_idx)
  }
  
  # Also plot calibration weight assumptions (only need do this once)
  plot_calibration_weights(o, "Calibration weights")
  
  # Save final result
  save_calibration(o, r_idx)
}

# ---------------------------------------------------------
# Initiate fit list and perform a few checks on input yaml
# ---------------------------------------------------------
setup_calibration = function(o) {
  
  # Load parsed parameters from yaml file
  p = parse_yaml(o, scenario = "baseline")$parsed
  
  # Shorthand for calibration options
  opts = p$calibration_options
  
  # ---- Determine type of calibration to perform ----
  
  message(" - Calibration method: ", p$calibration_type)
  
  # Use logical flags rather than string compare for calibration_type
  if (p$calibration_type == "Re_user")  do = list(.use_data = FALSE, .use_Re = TRUE)
  if (p$calibration_type == "Re_data")  do = list(.use_data = TRUE,  .use_Re = TRUE)
  if (p$calibration_type == "epi_data") do = list(.use_data = TRUE,  .use_Re = FALSE)
  
  # Also use a flag for whether we'll be working with synthetic data
  do$.use_synthetic = do$.use_data && opts$data_synthetic
  
  # ---- Sanity checks on input yaml ----
  
  # Metrics to be fitted over some time frame
  fit_metrics = get_fit_metrics(p)
  fit_days    = get_fit_days(p)
  
  # Required metrics not currently being reported
  missing_metrics = fit_metrics[!fit_metrics %in% p$metrics$df$metric]
  
  # Throw error if required metrics are missing
  if (length(missing_metrics) > 0)
    stop("For calibration type '", p$calibration_type, 
         "' you must report the following metrics: ", 
         paste(missing_metrics, collapse = ", "))
  
  # Throw error if n_days from yaml file is not at least max of fit_days
  if (p$n_days < fit_days)
    stop("Model calibration requires n_days >= ", fit_days)
  
  # Objective function scaling doesn't work well when fitting Re
  # if (do$.use_Re && p$objective_scaling != "none")
  #   stop("It is strongly recommended to use objective_scaling: 'none' when fitting Re")
  
  # ---- Construct list for fitting details ----
  
  # Initiate fit list (with easy access to calibration options)
  fit = c(do, list(input = p, opts = opts))
  
  # Extract parameters to calibrate
  fit_df = list2dt(fit$input$calibration_parameters) %>%
    mutate(lower = as.numeric(lower), 
           upper = as.numeric(upper))
  
  # Shorthand for fitted parameters
  fit$params = fit_df$param
  
  # Bounds of all parametrs to be fitted, as matrix
  fit$bounds = as.matrix(fit_df[, .(lower, upper)])
  
  # ---- Checks on synthetic data values ----
  
  # No need if not using synthetic data
  if (!do$.use_synthetic) 
    fit$synthetic = NA
  
  # Check whether this is necessary
  if (do$.use_synthetic) {
    
    # Add column if this doesn't already exist
    if (!"synth" %in% names(fit_df))
      fit_df$synth = NA
    
    # Ensure value is within bounds, use mean if doesn't exist, and convert to list
    fit$synthetic = fit_df %>%
      mutate(synth = pmax(synth, lower), 
             synth = pmin(synth, upper), 
             synth = ifelse(is.na(synth), lower+(upper-lower)/2, synth)) %>%
      pull(synth) %>%
      setNames(fit_df$param) %>%
      as.list()
  }
  
  # ---- Population scaler ----
  
  # Set value in fit list
  fit$pop_scaler = opts$country_pop / p$population_size
  
  return(fit)
}

# ---------------------------------------------------------
# Save final file with both 'emulated' and 'simulated' best params
# ---------------------------------------------------------
save_calibration = function(o, rx_idx) {
  
  # Load result of final adaptive sampling round
  rx_fit = try_load(o$pth$fitting, paste0(rx_idx, "_result"))
  
  # Most of the items can be copied directly
  fit = rx_fit[qc(input, params, bounds, data, target, synthetic, pop_scaler)]
  
  # ---- Best emulated parameter set ----
  
  # Take straight from optimisation output
  fit$best$emulated = rx_fit$result
  
  # ---- Best simulated parameter set ----
  
  # Load all simulated samples
  all_samples = try_load(o$pth$fitting, "rx_samples")
  
  # ID of parameter set with lowest objective value (mean across seeds)
  best_id = all_samples %>%
    mutate(paramset_id = get_paramset_id(param_id)) %>%
    group_by(paramset_id) %>%
    summarise(obj_value = mean(obj_value)) %>%
    slice_min(obj_value, n = 1) %>%
    pull(paramset_id)
  
  # Extract parameter set associated with this ID
  best_simulated = all_samples %>%
    mutate(paramset_id = get_paramset_id(param_id)) %>%
    filter(paramset_id == best_id) %>%
    select(all_of(fit$params)) %>%
    unique()
  
  # Store as list for consistency 
  fit$best$simulated = as.list(best_simulated)
  
  # Save this final result
  saveRDS(fit, file = paste0(o$pth$fitting, "fit_result.rds"))
}

# ---------------------------------------------------------
# Load final fit file and select 'best' parameter set
# ---------------------------------------------------------
load_calibration = function(o, ...) {
  
  # Load model fitting result - throw an error if it doesn't exist
  err_msg = "Cannot find a fitting file for this analysis - have you run step 1?"
  fit_result = try_load(o$pth$fitting, "fit_result", msg = err_msg, ...)
  
  # Leave only the 'best' parameter set to avoid confusion
  fit_result$best = fit_result$best[[o$best_param_set]]
  
  return(fit_result)
}

# ---------------------------------------------------------
# Generate a set of parameter sets to simulate
# ---------------------------------------------------------
sample_parameters = function(o, fit, r_val) {
  
  message("  > Sampling parameters")
  
  # Shorthand for emulator and adaptive sampling options
  em_opts = fit$input$emulator
  as_opts = fit$input$adaptive_sampling
  
  # First iteration is simple: sample a load of points
  if (r_val == 0)
    samples = tgp::lhs(em_opts$init_samples, fit$bounds)
  
  # On subsequent iterations, evaluate emulator and determine expected improvement
  if (r_val > 0) {
    
    # Load emulator produced in the previous round
    emulator_file = paste0("r", r_val - 1, "_emulator")
    emulator = try_load(o$pth$fitting, emulator_file)
    
    # We'll evaluate the emulator at a load of points and take the best ones
    eval_samples = tgp::lhs(as_opts$acquisition_points, fit$bounds) %>%
      as_named_dt(fit$params)
    
    # Evalute the points, calculate expected improvement, and take the most likely
    samples = eval_samples %>%
      normalise_params(fit$params, fit$bounds) %>% 
      as.matrix() %>% 
      crit_EI(emulator) %>% 
      as_named_dt("expected_improvement") %>% 
      cbind(eval_samples) %>%
      arrange(-expected_improvement) %>%
      filter(expected_improvement > as_opts$min_value_accepted) %>%
      slice(1 : min(n(), as_opts$max_samples_accepted)) %>%
      select(-expected_improvement) %>%
      as.matrix()
    
    # NOTE: We previously removed points that were very close, but have now scrapped that
  }
  
  # Number of samples to simulate (each one for multiple seeds)
  n_samples = nrow(samples)
  
  # If none have been found, provide trivial output
  if (n_samples == 0) {
    param_ids = character()
    
    # Warn the user that we'll be breaking out (min_value_accepted may be too high)
    message("  ! No suitable parameter sets identified - breaking out of adaptive sampling")
  }
  
  # Only continue if we have new points worth simulating  
  if (n_samples > 0) {
    
    # Expand for each seed and prepare for simulation
    paramset_df = samples %>%
      as_named_dt(fit$params) %>%
      mutate(param_idx = 1 : n_samples, 
             round     = r_val) %>%
      expand_grid(seed = 1 : em_opts$seeds) %>%
      mutate(param_id = get_param_id(round, param_idx, seed)) %>%
      select(param_id, round, seed, all_of(fit$params)) %>%
      setDT()
    
    # Extract parameter IDs
    param_ids = paramset_df$param_id
    
    # Save simulation df for reference on cluster
    saveRDS(paramset_df, file = paste0(o$pth$fitting, "r", r_val, "_paramsets.rds"))
  }
  
  return(param_ids)
}

# ---------------------------------------------------------
# Simulate newly-generated parameter sets (uses the cluster)
# ---------------------------------------------------------
simulate_parameters = function(o, r_idx, param_ids) {
  
  # Full path to all output files to be produced
  output_files = paste0(o$pth$fit_samples, param_ids, ".rds")
  
  # If we don't want to overwrite, remove reference to files that already exist
  if (o$overwrite_samples == FALSE)
    param_ids = param_ids[!file.exists(output_files)]
  
  # TODO: See if we can reuse check_existing() from scenarios.R here
  
  # Number of simulations to run (one per cluster job)
  n_simulations = length(param_ids)
  
  # ---- Simulate model on the cluster ----
  
  message("  > Simulating ", thou_sep(n_simulations), " samples")
  
  # Skip this process if nothing to run
  if (n_simulations > 0) {

    # Append reference to round index in job description string
    job_type = paste0("fitting::", r_idx)

    # Submit all jobs to the cluster (see auxiliary.R)
    submit_cluster_jobs(o, n_simulations, "bash_submit.sh", job_type)

    # Throw an error if any cluster jobs failed (see auxiliary.R)
    err_tol = floor(n_simulations * o$sample_err_tol)
    stop_if_errors(o$pth$log, o$err_file, err_tol = err_tol)

    # Remove all log files if desired (generally a good idea unless debugging)
    if (o$rm_cluster_log) unlink(paste0(o$pth$log, "*"), force = TRUE)
  }
  
  # ---- Load and summarise output ----
  
  message("  > Concatenating simulation outcomes")
  
  # Load all files into one long datatable and combine with samples
  output_list = lapply(output_files, function(x) try(readRDS(x), silent = TRUE))
  output_df   = rbindlist(output_list[unlist(lapply(output_list, is.data.table))])
  
  # Overwrite samples file for plotting purposes
  saveRDS(output_df, file = paste0(o$pth$fitting, r_idx, "_output.rds"))
}

# ---------------------------------------------------------
# Calculate how well model output matches the fitting target
# ---------------------------------------------------------
quality_of_fit = function(o, fit, r_idx, do_plot = FALSE) {
  
  message("  > Calculating quality of fit")
  
  # Load simulated parameter sets and associated model output
  paramset_df = try_load(o$pth$fitting, paste0(r_idx, "_paramsets"))
  output_df   = try_load(o$pth$fitting, paste0(r_idx, "_output"))
  
  # Model output may need to be scaled to get to per 100k people
  scaler   = 1e5 / fit$input$population_size
  scale_df = fit$input$metrics$df %>% 
    select(metric, scale) %>%
    mutate(scaler = ifelse(scale, scaler, 1))
  
  # Apply scaler to model output and remove burn in phase
  output_df %<>%
    left_join(scale_df, by = "metric") %>%
    mutate(value = value * scaler) %>%
    select(-scale, -scaler) %>%
    filter(date > fit$opts$data_burn_in) 
  
  # Change daily data to different time period if necessary
  if (fit$opts$data_period != "day") {
    
    # Formatting the dates
    all_dates = seq(fit$opts$data_start, fit$opts$data_end, by = "day")
    dates_df  = data.table(date = all_dates, 
                           day  = 1 : length(all_dates)) %>%
      filter(day > fit$opts$data_burn_in)
    
    # Transform the data
    output_df %<>%
      mutate(date = as_date(date, origin = fit$opts$data_start), 
             date = ceiling_date(date, fit$opts$data_period)) %>%
      group_by(param_id, round, seed, metric, date) %>%
      summarise_at("value", sum, na.rm = TRUE) %>%
      ungroup() %>%
      inner_join(dates_df, by = "date") %>%
      select(param_id, round, seed, metric, date = day, value) %>%
      setDT()
  }
  
  # When fitting to Re
  if (fit$.use_Re) {
    
    # Take the mean over dates and calculate squared difference between model and data
    obj_value_df = output_df %>%
      group_by(param_id, round, seed) %>%
      summarise(value = mean(value)) %>%
      ungroup() %>%
      mutate(obj_value = (fit$target - value) ^ 2, 
             obj_value = log(1 + obj_value * 100) + 1) %>%
      select(-value) %>%
      setDT()
  }
  
  # When fitting to epi data
  if (!fit$.use_Re) {
    
    # Cumulatively sum data and append weights
    data_df = fit$data %>%
      group_by(metric) %>%
      mutate(value = cumsum(value)) %>%
      ungroup() %>%
      format_weights(fit$input) %>%
      setDT()

    # Also cumulatively sum model output
    model_df = output_df %>%
      group_by(param_id, round, seed, metric) %>%
      mutate(value = cumsum(value)) %>%
      ungroup() %>%
      setDT()
    
    # Sum the log weighted squared difference between model and data (all metrics, all time points)
    obj_value_df = model_df %>%
      inner_join(data_df, by = c("metric", "date")) %>%
      group_by(metric) %>%
      mutate(value  = value  / max(target), 
             target = target / max(target)) %>%
      ungroup() %>%
      mutate(err = weight * (target - value) ^ 2) %>%
      group_by(param_id, round, seed) %>%
      summarise(obj_value = log(sum(err) + 1)) %>% 
      ungroup() %>%
      setDT()
  }
  
  # Join parameter sets to obj values and normalise error
  samples_df = paramset_df %>%
    inner_join(obj_value_df, by = c("round", "param_id", "seed")) %>%
    mutate(obj_scale_none   = obj_value, 
           obj_scale_weak   = 1 - (min(obj_value) / obj_value),
           obj_scale_strong = 1 - (min(obj_value) / obj_value) ^ 2) %>%
    select(all_of(names(paramset_df)), 
           obj_value = !!paste0("obj_scale_", fit$input$objective_scaling))
  
  # Save simulation df for reference on cluster
  saveRDS(samples_df, file = paste0(o$pth$fitting, r_idx, "_samples.rds"))
}

# ---------------------------------------------------------
# Select all parameter sets that we'll use to train an emulator
# ---------------------------------------------------------
select_samples = function(o, r_val) {
  
  # NOTE: At the moment this is a pretty trivial function, but could
  #       be used to filter out unlikely regions of parameter space
  
  # Initiate list of samples
  samples_list = list()
  
  # Iterate through all already simulated rounds
  for (i in 0 : r_val) {
    
    # Load samples from each round
    round_path    = paste0("r", i, "_samples")
    round_samples = try_load(o$pth$fitting, round_path)
    
    # Store in list as we iterate
    samples_list[[i + 1]] = round_samples
  }
  
  # Compile all samples into a datatable
  all_samples = rbindlist(samples_list)
  
  # Save all these samples with a slightly different naming convention
  saveRDS(all_samples, file = paste0(o$pth$fitting, "rx_samples.rds"))
}

# ---------------------------------------------------------
# Extract and select calibration weightings for metrics, time, and peaks
# ---------------------------------------------------------
format_weights = function(data, model_input, all = FALSE) {
  
  # Shorthand for calibration weights list
  w = model_input$calibration_weights
  
  # Weights for each metric
  w_metric = as.data.table(w$metric) %>% 
    pivot_longer(cols = everything(), 
                 names_to  = "metric", 
                 values_to = "w_metric")
  
  # All possible time and peak weights
  weight_all_df = data %>%
    filter(!is.na(value)) %>%
    select(date, metric, target = value) %>%
    group_by(metric) %>%
    # Punish errors at more recent dates...
    mutate(w_time_none   = 1, 
           w_time_weak   = seq(from = 0.5, to = 1, length.out = n()), 
           w_time_strong = seq(from = 0,   to = 1, length.out = n()), 
           w_time_exponential = rev(dexp((1 : n()) / n(), rate = 5)) / 
                                      max(dexp((1 : n()) / n(), rate = 5))) %>%
    # Punish errors at peaks...
    mutate(w_peak_none   = 1,
           w_peak_weak   = (target / max(target, na.rm = TRUE)) ^ 0.5, 
           w_peak_strong = (target / max(target, na.rm = TRUE)) ^ 2) %>%
    # Join with metric weights...
    ungroup() %>%
    left_join(w_metric, by = "metric") %>%
    setDT()
  
  # We may want this full dataframe for plotting purposes
  if (all == TRUE)
    return(weight_all_df)
  
  # Select the types of weights we want, discard the rest, and multiply for overall weight
  weight_df = weight_all_df %>%
    select(date, metric, target, w_metric, 
           w_time = !!paste0("w_time_", w$time),
           w_peak = !!paste0("w_peak_", w$peak)) %>%
    mutate(weight = w_metric * w_time * w_peak) %>%
    select(-starts_with("w_"))
  
  return(weight_df)
}

# ---------------------------------------------------------
# Create consistent IDs for parameter sets
# ---------------------------------------------------------
get_param_id = function(round, param, seed) {
  
  # Pad values for file name consistency
  param_id = paste0("r", sprintf("%02i", round),
                    "p", sprintf("%04i", param),
                    "s", sprintf("%02i", seed))
  
  return(param_id)
}

# ---------------------------------------------------------
# Reduce param_id for reference to round and paramset only
# ---------------------------------------------------------
get_paramset_id = function(param_id) {
  
  # Strip away reference to seed number
  paramset_id = str_remove(param_id, "s[0-9]+$")
  
  return(paramset_id)
}

# ---------------------------------------------------------
# Get days we need to simulate for for fitting purposes
# ---------------------------------------------------------
get_fit_days = function(p) {
  
  # If fitting to Re
  if (p$calibration_type != "epi_data")
    fit_days = p$Re_fit$start_day + 
      p$Re_fit$n_days + 
      p$Re_window - 1
  
  # If fitting to epi data a longer period is likely needed
  if (p$calibration_type == "epi_data")
    fit_days = p$calibration_options$data_days + 
      p$calibration_options$data_burn_in
  
  return(fit_days)
}

# ---------------------------------------------------------
# Which metrics are to be reported for fitting purposes
# ---------------------------------------------------------
get_fit_metrics = function(p) {
  
  # If fitting to Re: only one metric needed
  if (p$calibration_type != "epi_data")
    fit_metrics = "Re"
  
  # If fitting to epi data
  if (p$calibration_type == "epi_data") {
    
    # All non-trivial epi metrics required
    fit_metrics = p$calibration_weights$metric
    fit_metrics = names(fit_metrics[fit_metrics > 0])
  }
  
  return(fit_metrics) 
}

# ---------------------------------------------------------
# Alter key parameters in yaml file for fitting simulation
# ---------------------------------------------------------
apply_fit = function(y, fit_list) {
  
  # Only needed if fitting OR using fitted parameters
  if (!is.null(fit_list)) {
    
    # Check if we are fitting here - indicated by .perform_fit
    if (isTRUE(fit_list$.perform_fit)) {
      
      # In this case, we only need to run the bare essentials...
      
      # Bare essential metrics
      fit_metrics = get_fit_metrics(y)
      
      # Function for removing any disaggregation of metrics 
      reduce_fn = function(x) {if (!is.null(x$by)) {x$by = "none"}; return(x)}
      
      # Extract only the metrics of interest and removing any disaggregation
      y$model_metrics = lapply(y$model_metrics[fit_metrics], reduce_fn)
      
      # Only need to simulate for a shortened period
      fit_list$n_days = get_fit_days(y)
      
      # Now safely remove flag, it's done it's job
      fit_list$.perform_fit = NULL
    }
    
    # ---- Alter parameters ----
    
    # Loop through parameters to alter
    for (param in names(fit_list)) {
      
      # Check whether parameter we want to overwrite already exists
      eval_str("param_valid = is.numeric(y$", param, ")")
      
      # Check parameter is valid - throw error if not
      if (!param_valid)
        stop("Attempting to alter invalid parameter '", param, "' during model fitting")
      
      # Re-define paramter value (using eval to enable change of listed items)
      eval_str("y$", param, " = fit_list[[param]]")
    }
  }
  
  return(y)
}

