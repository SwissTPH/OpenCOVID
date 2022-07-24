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
                   fit  = fit_result$fit_output,
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
    results_list = list()
    
    # Loop through simulations and load model output
    for (sim_id in sim_ids) {
      sim_result = try_load(o$pth$simulations, sim_id, throw_error = !o$impute_failed_jobs)
      
      # Convert to datatable and append sample number
      results_list[[sim_id]] = sim_result$output
      
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
    if (length(results_list) == 0)
      stop("No results available for scenario '", scenario_name, "'")
    
    # Otherwise values will be imputed
    impute_sims = setdiff(sim_ids, names(results_list))
    
    # Warn the user if this is required
    if (length(results_list) < 0)
      warning("Imputing values for missing simulations: ", 
              paste(impute_sims, collapse = ", "))
    
    # Initiate analysis list
    result = list(analysis_name = o$analysis_name,
                  scenario_name = scenario_name,
                  time_stamp = format(Sys.time(), "%Y%m%d_%H%M"), 
                  yaml       = sim$yaml,
                  input      = sim$input, 
                  network    = sim$network)
    
    # Create key summary statistics across seeds
    result$output = rbindlist(results_list) %>% 
      filter(!is.na(value)) %>% 
      group_by(date, metric, grouping, group, scenario) %>% 
      summarise(mean   = mean(value),
                median = quantile(value, 0.5),
                lower  = quantile(value, o$quantiles[1]),
                upper  = quantile(value, o$quantiles[2])) %>% 
      as.data.table()
    
    message("Scenario summary complete")
    
    # Store summarised result
    saveRDS(result, file = paste0(o$pth$scenarios, scenario_name, ".rds"))
  }
}

