############################################################
# CLUSTER JOBS
#
# Perform a job from a cluster array submission.
#
############################################################

# ---------------------------------------------------------
# Perform cluster job of choice (see job_type input)
# ---------------------------------------------------------
run_cluster_job = function(o, job_type, task_id) {
  
  # ---- Run fitting samples -----
  
  # Run fitting samples used to train model emulator
  if ((grepl("fitting", job_type))) {
    
    # Parse (re)sampling number from input
    r_idx = unlist(str_split(job_type, "::"))[2]
    r_val = as.numeric(str_remove(r_idx, "r"))
    
    # Load parameter set samples
    sim_df = try_load(o$pth$fitting, paste0(r_idx, "_paramsets"))
    
    # Select parameter set associated with this task ID (all seeds)
    param_df = sim_df[task_id, ]
    param_id = param_df$param_id
    
    # Parameter values in list format (for input into model)
    param_list = param_df %>%
      select(-round, -param_id, -seed) %>%
      unique() %>%
      as.list()
    
    # Append flag that we want to perform fit
    fit_list = list.append(param_list, .perform_fit = TRUE)
    
    # If uncertainty defined, take the average over the distribution(s)
    uncert_list = sample_average(o)  # See uncertainty.R
    
    message(" - Running model")
    
    # Simulate model with this parameter set
    result = model(o, "baseline", 
                   seed    = param_df$seed, 
                   fit     = fit_list, 
                   uncert  = uncert_list, 
                   verbose = "none")
    
    # Extract metrics of interest over time points of interest
    output_df = result$output %>%
      filter(!is.na(value)) %>%
      select(seed, metric, date, value) %>%
      mutate(param_id = !!param_id, 
             round    = r_val, 
             .before  = 1)
    
    # Save sample output as an RDS file
    saveRDS(output_df, paste0(o$pth$fit_samples, param_id, ".rds"))
  }
  
  # ---- Run simulations -----
  
  # Run all user-defined scenarios
  if (job_type == "scenarios") {
    
    # Load full set of simulations to run
    sim_df = try_load(o$pth$simulations, "all_simulations")
    
    # Select analysis to run based on job ID
    this_sim = sim_df[task_id, ]
    
    # Check whether there are multiple parameter sets to sample
    if (length(unique(sim_df$param_set)) > 1) {
      
      # If yes, load uncertainty datatable
      uncert_df = try_load(o$pth$uncertainty, "uncertainty")
      
      # Select values relevant for this simulation
      this_uncert = uncert_df[param_set == this_sim$param_set, ]
      
      # Convert to list
      uncert_list = this_uncert$value %>% 
        setNames(this_uncert$param) %>%
        as.list()
      
    } else {  # Use trivial list if no parameter uncertainty
      uncert_list = NULL
    }
    
    # Loaded fitted model parameters
    fit_list = load_calibration(o)$best
    
    # Simulate the model for this scenario and this seed (see model.R)
    result = model(o, this_sim$scenario, 
                   seed    = this_sim$seed_num, 
                   fit     = fit_list,
                   uncert  = uncert_list, 
                   verbose = "none") 
    
    # Save model inputs and outputs as an RDS file
    saveRDS(result, paste0(o$pth$simulations, this_sim$sim_id, ".rds"))
  }
  
  # ---- Summarise simulations ----
  
  # Summarise all simulations for a scenarios
  if (job_type == "summarise") {
    
    # Load full set of simulations to summarise
    sim_df = try_load(o$pth$simulations, "all_simulations")
    
    # Scenario to summarise
    scenario_name = unique(sim_df$scenario)[task_id]
    
    # All simulations to summarise over
    sim_ids = sim_df %>% 
      filter(scenario == scenario_name) %>%
      pull(sim_id)
    
    # Preallocate list for model outcomes
    output_list = time_list = list()
    
    # Loop through simulations and load model output
    for (sim_id in sim_ids) {
      sim_result = try_load(o$pth$simulations, sim_id, throw_error = !o$impute_failed_jobs)
      
      # Store model output and time taken for each simulation
      output_list[[sim_id]] = sim_result$output
      time_list[[sim_id]]   = sim_result$time_taken
      
      # Check if simulation was succesfully loaded
      if (!is.null(sim_result)) {
        
        message("Simulation '", sim_id, "' succesfully loaded")
        
        # Store simulation inputs (in case last sim in set failed and gives NULL)
        sim = list(yaml    = sim_result$yaml,
                   input   = sim_result$input, 
                   network = sim_result$network)
      }
    }
    
    # Throw an error if no results found
    if (length(output_list) == 0)
      stop("No results available for scenario '", scenario_name, "'")
    
    # Otherwise any missing values will be imputed
    impute_sims = setdiff(sim_ids, names(output_list))
    
    # Warn the user if this is required
    if (length(impute_sims) > 0)
      warning("Imputing values for missing simulations: ", 
              paste(impute_sims, collapse = ", "))
    
    # Convert time taken of each sim into seconds (lubridate periods don't summarise well)
    time_vec = unlist(lapply(time_list, period_to_seconds))
    
    # Mean simulation time 
    time_mean = seconds_to_period(round(mean(time_vec)))
    
    # Min and max simulation times
    time_min = seconds_to_period(min(time_vec))
    time_max = seconds_to_period(max(time_vec))
    
    # Compile this info into a string
    time_str = paste0(time_mean, " (", time_min, " - ", time_max, ")")
    
    # Initiate analysis list - we'll add summarised model output to this
    result = list(analysis_name = o$analysis_name,
                  scenario_name = scenario_name,
                  time_stamp = format(Sys.time(), "%Y%m%d_%H%M"), 
                  time_taken = time_str, 
                  time_raw   = time_vec, 
                  yaml       = sim$yaml,
                  input      = sim$input, 
                  network    = sim$network)
    
    # Bind model output into single datatable
    raw_output = rbindlist(output_list)[!is.na(value), ]
    
    # Aggregate any groupings (for appropriate metrics)
    raw_df = aggregate_results(result$input, raw_output)
    
    # Save this raw (ie unsummarised) form of model output
    saveRDS(raw_df, file = paste0(o$pth$scenarios, scenario_name, "_raw.rds"))
    
    # Aggregative and summarise raw model output (see postprocess.R)
    result = process_results(result, raw_output)
    
    # Store summarised result
    saveRDS(result, file = paste0(o$pth$scenarios, scenario_name, ".rds"))
  }
}

