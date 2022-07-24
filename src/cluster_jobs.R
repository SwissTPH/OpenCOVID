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
  if (grepl("fitting", job_type)) {
    
    # Load full factorial set of cantons, seeds, and simulation names
    sim_df = readRDS(paste0(o$pth$fit_samples, "all_samples.rds"))
    
    # Select sims associated with this param ID
    param_df = filter(sim_df, param_id == task_id)
    
    # Preallocate output dataframe
    output_df = param_df %>%
      select(param_id, seed) %>%
      mutate(r_eff = NA)
    
    # Parameter values in list format (for input into model)
    param_list = param_df %>%
      select(-param_id, -seed) %>%
      unique() %>%
      as.list()
    
    # Append number of days to simulate for (not many needed for fitting Reff)
    fit_list = list.append(param_list, n_days = max(o$fit_days))
    
    # Loop through seeds
    for (seed in param_df$seed) {
      
      message("  > Seed ", seed, " of ", o$emulator_seeds)
      
      # TODO: Put model inside trycatch - it's ok if we get occasional error
      
      # Simulate model
      result = model(o, "baseline", 
                     seed = seed, 
                     fit  = fit_list, 
                     verbose = "none")
      
      # Extract mean R_eff over time points of interest
      r_eff_df = filter(result$output, date %in% o$fit_days)
      
      # Store in output dataframe
      output_df$r_eff[seed] = mean(r_eff_df$value)
    }
    
    # Save sample outputs as an RDS file
    saveRDS(output_df, paste0(o$pth$fit_samples, "sample_", task_id, ".rds"))
  }
  
  # ---- Run simulations -----
  
  # Run all user-defined scenarios
  if (grepl("scenarios", job_type)) {
    
    # Load full factorial set of cantons, seeds, and simulation names
    sim_df = readRDS(paste0(o$pth$simulations, "all_simulations.rds"))
    
    # Select analysis to run based on job ID
    this_sim = sim_df[task_id, ]
    
    # Loaded fitted model parameters
    fit_result = readRDS(paste0(o$pth$fitting, "fit_result.rds"))
    
    # Simulate the model for this scenario and this seed (see model.R)
    result = model(o, this_sim$scenario, 
                   seed = this_sim$seed, 
                   fit  = fit_result$result,
                   verbose = "date")
    
    # Save model inputs and outputs as an RDS file
    saveRDS(result, paste0(o$pth$simulations, this_sim$sim_id, ".rds"))
  }
  
  # ---- Summarise simulations ----
  
  # Summarise all simulations for a scenarios
  if (grepl("summarise", job_type)) {
    
    # Load full set of simulations
    sim_df = readRDS(paste0(o$pth$simulations, "all_simulations.rds"))
    
    # Scenario to summarise
    scenario_name = unique(sim_df$scenario)[task_id]
    
    # All simulations to summarise over
    sim_ids = sim_df %>% 
      filter(scenario == scenario_name) %>%
      pull(sim_id)
    
    # Preallocate list for model outcomes
    output_list = list()
    
    # Loop through simulations and load model output
    for (sim_id in sim_ids) {
      sim_result = try_load(o$pth$simulations, sim_id, throw_error = !o$impute_failed_jobs)
      
      # Convert to datatable and append sample number
      output_list[[sim_id]] = sim_result$output
      
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
    
    # Otherwise values will be imputed
    impute_sims = setdiff(sim_ids, names(output_list))
    
    # Warn the user if this is required
    if (length(impute_sims) > 0)
      warning("Imputing values for missing simulations: ", 
              paste(impute_sims, collapse = ", "))
    
    # Initiate analysis list - we'll add summarised model output to this
    result = list(analysis_name = o$analysis_name,
                  scenario_name = scenario_name,
                  time_stamp = format(Sys.time(), "%Y%m%d_%H%M"), 
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

