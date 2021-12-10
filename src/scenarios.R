###########################################################
# SCENARIOS
#
# Simulates all scenarios defined in yaml file for this 
# analysis, including baseline.
#
###########################################################

# ---------------------------------------------------------
# Parent function for running scenarios
# ---------------------------------------------------------
run_scenarios = function(o) {
  
  # Only continue if specified by do_step
  if (!is.element(2, o$do_step)) return()
  
  message("* Running all scenarios")
  
  # Check calibration file is suitable
  check_fit(o)
  
  # ---- Set up full set of simulations ----
  
  # All scenarios to simulate
  all_scenarios = names(parse_yaml(o, "*read*", read_array = TRUE))
  
  # We'll run each simulation for different seeds
  seed_id = str_pad(1 : o$n_seeds_analysis, 4, pad = "0")
  
  # ... and potentially different parameter sets (sampled from fitted posteriors)
  param_id = str_pad(0 : o$n_parameter_samples, 4, pad = "0")
  
  # Full factorial set of cantons, parameter sets, seeds, and simulations
  all_simulations = expand_grid(param_id = param_id, 
                                seed     = seed_id, 
                                scenario = all_scenarios)
  
  # Use all of this to create unique simulation identifiers
  sim_df = all_simulations %>% 
    mutate(sim_id = unite(all_simulations, "x", sep = "_")$x, 
           sim_id = paste0("s", sim_id), 
           seed = as.numeric(seed)) %>%
    as.data.table()
  
  # Save to file for use on the cluster
  saveRDS(sim_df, file = paste0(o$pth$simulations, "all_simulations.rds"))
  
  # ---- Simulate on the cluster ----
  
  # We may not want/need to run them all
  n_simulations = check_existing(o, sim_df)

  # Do we need to simulate at all?
  if (n_simulations > 0) {

    # Submit all jobs to the cluster (see myRfunctions.R)
    submit_cluster_jobs(o, n_simulations, "bash_submit.sh", "scenarios")

    # Throw an error if any cluster jobs failed (see myRfunctions.R)
    stop_if_errors(o$pth$log, o$err_file, err_tol = 1)

    # Remove all log files if desired (generally a good idea unless debugging)
    if (o$rm_cluster_log) unlink(paste0(o$pth$log, "*"), force = TRUE)
  }
  
  # ---- Summarise and aggregate results ----
  
  # If no new simulations we may be able to skip this step
  if (o$force_summarise == FALSE && n_simulations == 0) {

    # Summarised file of each scenario
    summarised_files = paste0(o$pth$scenarios, all_scenarios, ".rds")

    # If all exist and we haven't re-simulated, we're done here
    if (all(file.exists(summarised_files)))
      return()
  }
  
  message(" - Quantifying parameter and stochastic uncertainty")
  
  # Re-save full set of simulations for use on the cluster
  #
  # NOTE: Needed as original may have been overwritten in check_existing function
  saveRDS(sim_df, file = paste0(o$pth$simulations, "all_simulations.rds"))

  # Number of unique scenarios... we'll summarise each on the cluster
  n_jobs = length(unique(sim_df$scenario))

  # Submit all summarising jobs to the cluster (see myRfunctions.R)
  submit_cluster_jobs(o, n_jobs, "bash_submit.sh", "summarise")

  # Throw an error if any cluster jobs failed (see myRfunctions.R)
  stop_if_errors(o$pth$log, o$err_file)
  
  # Missing simulations that have been imputed
  files_exist  = str_remove(list.files(o$pth$simulations), ".rds")
  missing_sims = filter(sim_df, !sim_id %in% files_exist)$sim_id
  
  if (length(missing_sims) > 0)
    warning("Imputed values for ", length(missing_sims), " missing simulation(s):\n", 
            paste(missing_sims, collapse = "\n"))
  
  # Remove all log files if desired (generally a good idea unless debugging)
  if (o$rm_cluster_log) unlink(paste0(o$pth$log, "*"), force = TRUE)
}

# ---------------------------------------------------------
# Check that calibration is still suitable
# ---------------------------------------------------------
check_fit = function(o) {
  
  # Load model fitting result - throw an error if it doesn't exist
  err_msg = "Cannot find a fitting file for this analysis - have you run step 1?"
  fit_result = try_load(o$pth$fitting, "fit_result", msg = err_msg)
  
  # Load current baseline - we need to check input hasn't changed
  baseline = parse_yaml(o, scenario = "baseline")
  
  # Remove items that it is ok to change (see o$fit_changeable_items)
  fit_input = list.remove(fit_result$fit_input, o$fit_changeable_items)
  baseline_input = list.remove(baseline$parsed, o$fit_changeable_items)

  # Are these two lists identical - if yes we are done here
  check_flag = setequal(baseline_input, fit_input)
  
  # If yes, we are done here
  if (check_flag == TRUE)
    return()
  
  # If not, try to find the issue and report back to user...
  
  # Create a standard error message - user will need to rerun step 1
  rerun_msg = "Your calibration file contains inconsistencies and must be rerun\n\n"
  
  # First check for missing items
  missing_items = symdiff(names(baseline_input), names(fit_input))
  
  # Throw an error if any parameters missing in either file
  if (length(missing_items) > 0)
    stop(rerun_msg, "! Parameters missing from calibration or baseline file: ", 
         paste(missing_items, collapse = ", "))
  
  # We'll then search for inconsistent parameters
  inconsistent_items = character()
  
  # Loop through all non-changeable parameters
  for (param in names(fit_input)) {
    
    # Check for consistency between baseline and fitted file
    check_flag = setequal(fit_input[[param]], baseline_input[[param]])
    
    # If any inconsistency is identified, store parameter name
    if (check_flag == FALSE)
      inconsistent_items = c(inconsistent_items, param)
  }
  
  # Throw an error if any parameters are inconsistent
  if (length(inconsistent_items) > 0)
    stop(rerun_msg, "! Inconsistent parameters in calibration and baseline files: ", 
         paste(inconsistent_items, collapse = ", "))
  
  # If we've got this far, we've identified a problem, but not sire where it is
  stop(rerun_msg, "! Unknown error")
}

# ---------------------------------------------------------
# Check if we can avoid rerunning any existing simulations
# ---------------------------------------------------------
check_existing = function(o, sim_df) {
  
  # Start by assuming we'll need to rerun everything
  run_df  = sim_df
  n_total = nrow(run_df)
  
  message(" - Total number of simulations: ", thou_sep(n_total))
  
  # Simulations that already exist (ie from previous analyses)
  files_exist = str_remove(list.files(o$pth$simulations), ".rds")
  sims_exist  = filter(sim_df, sim_id %in% files_exist)
  
  # Do we want to avoid overwriting any previously run simulations?
  if (o$overwrite_simulations == FALSE && nrow(sims_exist) > 0) {
    
    # We won't overwrite - only rerun any missing simulations any with inconsistent yamls
    if (o$check_yaml_consistency == TRUE) {
      
      message("  > Checking consistency of YAML files")
      
      # Load fitted parameters - we'll need to rerun if these have changed
      fit_result = try_load(o$pth$fitting, "fit_result")$fit_output
      
      # Assess each scenario that already exists
      for (scenario_exist in unique(sims_exist$scenario)) {
        
        # Only need to look at one simulation - load it and extract input list
        check_sim   = filter(sims_exist, scenario == scenario_exist)[1, ]
        check_input = try_load(o$pth$simulations, check_sim$sim_id)$input
        
        # How the input looks now - direct from the current yaml
        test_input = parse_yaml(o, scenario_exist, read_array = TRUE)$parsed
        
        # Apply calibrated parameters - important that these are the same
        test_input[names(fit_result)] = fit_result
        
        # Compare all items in these two lists - aside from calibrated parameters
        check_flag = setequal(check_input, test_input)
        
        # If lists are identical, we can safely skip
        if (check_flag == TRUE) {
          
          # IDs of sims that have already been run - might not be all
          skip_sims = run_df %>%
            filter(scenario == scenario_exist) %>%
            pull(sim_id) %>%
            intersect(files_exist)
          
          # Remove these from run_df 
          run_df = filter(run_df, !sim_id %in% skip_sims)
        }
      }
    }
    
    # Alternatively, don't bother checking yaml consistency - won't overwrite regardless
    if (o$check_yaml_consistency == FALSE)
      run_df = filter(run_df, !sim_id %in% files_exist)
    
    # Save to file for use on the cluster
    saveRDS(run_df, file = paste0(o$pth$simulations, "all_simulations.rds"))
  }
  
  # Updated number of simulations that need to be run / re-run
  n_run = nrow(run_df)
  
  # Report this back to the user
  message("  > Skipping: ", thou_sep(n_total - n_run))
  message("  > Simulating: ", thou_sep(n_run))
  
  return(n_run)
}

