###########################################################
# ANALYSIS
#
# Simulates fitted baseline analysis and any defined scenarios
# and/or strategies.
#
###########################################################

# ---------------------------------------------------------
# Parent function for running desired type of analysis
# ---------------------------------------------------------
run_analysis = function(o) {
  
  # Only continue if specified by do_step
  if (!is.element(3, o$do_step)) return()
  
  message("* Running analyses")
  
  # ---- Set up set of simulations to run ----
  
  # All scenarios and strategies to run
  #
  # NOTE: These are run for all selected cantons
  all_simulations = get_simulations(o)
  
  # Sample and store parameter sets for selected cantons
  param_id = sample_parameters(o)
  
  # We'll run each simulation for different seeds
  #
  # NOTE: Seed ID is not necessarily the seed number, that depends on baseline$best_seed
  seed_id = str_pad(1 : o$n_seeds_analysis, 4, pad = "0")
  
  # Full factorial set of cantons, parameter sets, seeds, and simulations
  analysis_df = expand_grid(canton   = o$cantons, 
                            param_id = param_id,
                            seed_id  = seed_id, 
                            sim_name = all_simulations$sim_name)
  
  # Use all of this to create unique simulation identifiers
  analysis_df$sim_id = unite(analysis_df, "x", sep = "_")$x
  
  # Convert to data.table and append strategy names
  analysis_df = data.table(analysis_df) %>%
    left_join(all_simulations, by = "sim_name")
  
  # Save to file for use on the cluster
  saveRDS(analysis_df, file = paste0(o$pth$sims, "all_simulations.rds"))
  
  # ---- Simulate on the cluster ----
  
  # Determine which jobs actually need to run
  n_simulations = which_jobs(o, analysis_df)
  
  # No need to run and check if zero
  if (n_simulations > 0) {
    
    # Submit all jobs to the cluster (see myRfunctions.R)
    submit_cluster_jobs(o, n_simulations, "bash_submit.sh", "analysis")
    
    # Throw an error if any cluster jobs failed (see myRfunctions.R)
    stop_if_errors(o$pth$code, o$err_file)
  }
  
  # ---- Summarise and aggregate results ----
  
  # Summarise and store final results
  summarise_results(o, analysis_df)
}

# ---------------------------------------------------------
# Generate lit of all scenarios and strategies to run
# ---------------------------------------------------------
get_simulations = function(o) {
  
  # Unpack names of all scenarios and strategies
  all_scenarios  = unname(unlist(o$scenarios))
  all_strategies = setdiff(o$strategies, "none")
  
  # Initiate set with scenarios (including baseline)
  all_simulations = data.table(sim_name = c("baseline", all_scenarios), 
                               strategy = NA)
  
  # Strategies contain multiple scenarios
  for (strategy in all_strategies) { 
    
    # Throw an error if not strategy not recognised recognised
    if (!exists(strategy, mode = "function"))
      stop("Strategy '", strategy, "' not well defined")
    
    # Run the strategy by calling the named function
    details = get(strategy)()
    
    # Full factorial of strategies
    combinations = merge(data.frame(roll_out = names(details$roll_out)), 
                         data.frame(npi_exit = names(details$npi_exit)))
    
    # Construct scenario ID from strategy name and elements
    combinations$id = paste0(strategy, ".", unite(combinations, "x", sep = "_")$x)
    
    # TODO: Use defaults if any names have not been defined...
    
    # Extract descriptive names from the strategy definitions
    roll_out_names = unlist(lapply(details$roll_out, "[[", "name"))
    npi_exit_names = unlist(lapply(details$npi_exit, "[[", "name"))
    
    # Crunch these together for a (reasonably) concise set of scenario names
    combinations$name = paste0(roll_out_names[combinations$roll_out], ", ", 
                               npi_exit_names[combinations$npi_exit])
    
    # Remove reference to names in nested details list
    details$roll_out = lapply(details$roll_out, function(x) x[names(x) != "name"])
    details$npi_exit = lapply(details$npi_exit, function(x) x[names(x) != "name"])
    
    # Format all of this into one list that can be read when running on the cluster
    this_strategy = list(strategy     = strategy,
                         combinations = combinations, 
                         details      = details)
    
    # Save to file for use on the cluster
    saveRDS(this_strategy, file = paste0(o$pth$strategy, "strategy.", strategy, ".rds"))
    
    # Append all combinations to simulations vector
    all_simulations = rbind(all_simulations, data.table(sim_name = combinations$id, 
                                                        strategy = strategy))
  }
  
  return(all_simulations)
}

# ---------------------------------------------------------
# Sample parameter sets for selected cantons
# ---------------------------------------------------------
sample_parameters = function(o) {
  
  # Construct vector of parameter IDs
  param_id = str_pad(0 : o$n_parameter_samples, 4, pad = "0")
  
  # Loop through selected cantons
  for (canton in o$cantons) {
    
    # Read MCMC calibration output file
    err_msg = paste0("Cannot find calibration file for canton ", canton, 
                     " - either run a calibration or use 'calibration_name' option")
    fit_result = try_load(o$pth$fitting, paste0(canton, "_calibration"), msg = err_msg) 
    
    # Store data used to create fit, along with some fitting info and options
    data_list = list(fit_data    = fit_result$data, 
                     fit_info    = fit_result[qc(fit_name, fit_cantons, time_stamp, param_table)], 
                     fit_options = fit_result$options)
    
    # Save to file for use on the cluster
    saveRDS(data_list, file = paste0(o$pth$analysis, canton, "_data.rds"))
    
    # Extract MCMC chain 
    mcmc_chain = fit_result$chain
    
    # Initiate list of parameter details for this canton
    param_list = list(data = fit_result$data)
    
    # Extract best estimate parameter set based on value of best_parameter_set
    if (o$best_parameter_set == "simulated") {
      
      # Best of all actually simulated parameter sets
      param_list$best_params = fit_result$best_simulated
      param_list$best_seed   = fit_result$best_sim_seed  # Seed which gave the best fit
    }
    
    # Extract best estimate parameter set based on value of best_parameter_set
    if (o$best_parameter_set == "emulated") {
      
      # Best parameter set identified by MCMC of emulated space
      param_list$best_params = get_best_pars(mcmc_chain)
      param_list$best_seed   = 1  # Any seed will do as none have yet been simulated
    }
    
    # Throw an error if there were any issues - likely invalid value of best_parameter_set
    if (is.null(param_list$best_params)) 
      stop("Best parameter set option '", o$best_parameter_set, "' not recognised")
    
    # Generate a load of indices to sample from MCMC chain
    samples_idx = sample(nrow(mcmc_chain), o$n_parameter_samples)
    
    # Sample the chain, bind with best params, and label with param IDs
    param_list$param_sets = mcmc_chain[samples_idx, ] %>% 
      select(-sampno, -lnlike) %>% 
      rbind(as.data.table(t(param_list$best_params))) %>% 
      mutate(param_id = rev(param_id), .before = 1) %>% 
      arrange(param_id)
    
    # Save to file for use on the cluster
    saveRDS(param_list, file = paste0(o$pth$analysis, canton, "_parameters.rds"))
  }
  
  return(param_id)
}

# ---------------------------------------------------------
# Determine which jobs actually need to run
# ---------------------------------------------------------
which_jobs = function(o, analysis_df) {
  
  # Simulations that already exist (ie from previous analyses)
  sims_exist = str_remove(list.files(o$pth$sims), ".rds")
  sims_exist = sims_exist[sims_exist %in% analysis_df$sim_id]
  
  # Total number of simulations to run
  n_jobs = nrow(analysis_df)   # Including any already simulated
  n_done = length(sims_exist)  # Number already simulated
  
  message(" - Running ", thou_sep(n_jobs), " model simulations")
  
  # Have any of these simulations already been run?
  if (n_done > 0) {
    
    # We will overwrite - explain to user
    if (o$overwrite_analysis == TRUE)
      message("  > Overwriting ", thou_sep(n_done), " previously completed simulations")
    
    # We won't overwrite - only rerun any missing simulations
    if (n_done > 0 && o$overwrite_analysis == FALSE) {
      message("  > Skipping ", thou_sep(n_done), " previously completed simulations")
      
      # Subset analyses to only those we haven't yet simulated
      analysis_df = filter(analysis_df, !sim_id %in% sims_exist)
      
      # Save to file for use on the cluster
      saveRDS(analysis_df, file = paste0(o$pth$sims, "all_simulations.rds"))
    }
  }
  
  # Total number of simulations we now need to run 
  n_simulations = nrow(analysis_df)
  
  return(n_simulations)
}

# ---------------------------------------------------------
# Summarise simulations and store final results
# ---------------------------------------------------------
summarise_results = function(o, analysis_df) {
  
  message(" - Quantifying parameter and stochastic uncertainty")
  
  # Number of unique scenarios/strategies... we'll summarise each on the cluster
  n_jobs = nrow(unique(select(analysis_df, canton, sim_name)))
  
  # Save to file for use on the cluster
  #
  # NOTE: We do this again here as original may have been overwritten in which_jobs
  saveRDS(analysis_df, file = paste0(o$pth$sims, "all_simulations.rds"))
  
  # We'll need some extra juice if a large enough number of simulations to summarise
  o$job_memory = o$force_memory
  
  # Submit all summarising jobs to the cluster (see myRfunctions.R)
  submit_cluster_jobs(o, n_jobs, "bash_submit.sh", "summarise")
  
  # Throw an error if any cluster jobs failed (see myRfunctions.R)
  stop_if_errors(o$pth$code, o$err_file)
  
  # ---- Aggregate strategies ----
  
  # We'll use time stamps to differentiate analyses
  time_stamp = format(Sys.time(), "%Y%m%d_%H%M")
  
  # Carry over key aspects from the fitting process, such as data and relevant options
  fit_inherit = c("fit_data", "fit_info", "fit_options")  # From data_list, loaded below
  
  # Loop through all cantons
  for (this_canton in o$cantons) {
    
    # Load info from fitting file - includes data and options used to generate fit
    data_list = readRDS(paste0(o$pth$analysis, this_canton, "_data.rds"))
    
    # Loop through all strategies we are analysing 
    all_strategies = unique(na.omit(unique(analysis_df$strategy)))
    for (this_strategy in all_strategies) {
      
      # Load details of this strategy from file
      strategy_info = readRDS(paste0(o$pth$strategy, "strategy.", this_strategy, ".rds"))
      
      # Store key details in a file which will contain all outcomes from all strategy combinations
      a_strategy = list(analysis_name = o$analysis_name, 
                        time_stamp    = time_stamp, 
                        combinations  = strategy_info$combinations, 
                        details       = strategy_info$details)
      
      # Loop through each of the simulation combinations and load result
      sim_list = list()
      for (this_sim_id in a_strategy$combinations$id) {
        
        # Extract model outcomes and append sim ID
        a_sim = readRDS(paste0(o$pth$scenario, this_canton, "_", this_sim_id, ".rds"))
        sim_list[[this_sim_id]] = a_sim$model_output
      }
      
      # Continuously concatenate all model outcomes
      a_strategy$model_output = rbindlist(sim_list)
      
      # Inherit options from this analysis run and the fitting process
      a_strategy$options = store_options(o, "analysis")
      a_strategy[fit_inherit] = data_list[fit_inherit]
      
      # Save as a compiled file - this makes life much easier when plotting
      saveRDS(a_strategy, file = paste0(o$pth$strategy, this_canton, "_", this_strategy, ".rds"))
    }
  }
}

