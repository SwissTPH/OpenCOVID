############################################################
# SIMULATE
#
# Simulate model from a cluster array or by calling directly.
#
############################################################

# ---------------------------------------------------------
# Main function: perform simulation of choice (see sim_type input)
# ---------------------------------------------------------
do_simulate = function(sim_type, task_id, o = NULL) {
  
  # Reload simulation options (see options.R)
  if (is.null(o)) o = set_options()
  
  # ---- Run calibration simulations -----
  
  # Calibration step: A load of LHC-sampled parameter sets
  if (grepl("calibration", sim_type)) {
    
    # Parse (re)sampling number from input
    sample_num = unlist(str_split(sim_type, "::"))[2]
    sample_pth = paste0(o$pth$samples, sample_num, file_sep())
    
    # Load full parameter set dataframe from file
    param_sets = readRDS(paste0(sample_pth, "samples", sample_num, ".rds"))
    
    # Filter out those already simulated in previous sampling iteration
    param_samples = param_sets$samples %>% 
      filter(sample_id == as.numeric(sample_num)) 
    
    # Select the parameter set defined by task_id
    param_values = param_samples %>% 
      select(-sim_id, -sample_id, -param_id, -seed, -likelihood) %>% 
      dplyr::slice(task_id)
    
    # Re-transform this parameter set to real scale
    param_set = unlist(param_values) * param_sets$scale + param_sets$center
    
    # Select sim ID defined by task_id
    this_sim_id = param_samples$sim_id[task_id]
    
    # Does this file already exist?
    file_name   = paste0(sample_pth, this_sim_id, ".rds")
    file_exists = file.exists(file_name)
    
    # We may just want to recalculate for existing simulations 
    #
    # NOTE: This is helpful if experimenting with different weightings
    if (file_exists == TRUE && o$overwrite_fit == FALSE && o$recalculate_fit == TRUE) {
      
      # Load model outcomes generated from previous simulations
      original_result = readRDS(file_name)
      
      # Parameter sets should be the same - if not then a new set of simulations likely required
      if (!identical(original_result$param_set, param_set))
        stop("Attempting to recalculate likelihood, but parameters sets appear inconsistent")
      
      # Recalculate likelihood using this 
      result = likelihood_model(o, param_set = param_set, quiet = TRUE, 
                                m = original_result$model_output)
      
    } else {
      
      # Also select seed number defined by task_id
      seed_num = param_samples$seed[task_id]
      
      # Run model and calculate likelihood given the data
      result = likelihood_model(o, param_set = param_set, seed_num = seed_num)
    }
    
    # Save single value as an RDS file (unless theres been an issue)
    #
    # NOTE: Surely better to insert this into a relational database
    if (!is.na(result$likelihood))
      saveRDS(result, file_name)
  }
  
  # ---- Run analysis simulations -----
  
  # Analysis step: Best fitting parameter set and samples from posteriors
  if (grepl("analysis", sim_type)) {
    
    # Load full factorial set of cantons, seeds, and simulation names
    analysis_df = readRDS(paste0(o$pth$sims, "all_simulations.rds"))
    
    # Select analysis to run based on task ID
    this_analysis = analysis_df[task_id, ]
    
    # Load parameter details associated with this canton
    param_list = readRDS(paste0(o$pth$analysis, this_analysis$canton, "_parameters.rds"))
    
    # Extract parameter set to use for this simulation
    this_params = param_list$param_sets %>%
      filter(param_id == this_analysis$param_id) %>%
      select(-param_id) %>% unlist()
    
    # Seed number to run this set of parameters with
    #
    # NOTE: We want seed_id = 1 to provide seed number that produced best fit
    this_seed = param_list$best_seed + as.numeric(this_analysis$seed_id) - 1
    
    # Load info from fitting file - includes data and options used to generate fit
    data_list = readRDS(paste0(o$pth$analysis, this_analysis$canton, "_data.rds"))
    
    # Override subset of o with options from fitting file (see store_options() function)
    o[names(data_list$fit_options)] = data_list$fit_options
    
    # Set up model inputs using fitting data
    p = get_parameters(o, data_list$fit_data, this_params)  # See parameters.R
    
    # If not a strategy, we simply alter future model parameters via alter_model_params
    if (is.na(this_analysis$strategy)) {
      p = alter_model_params(o, p, this_analysis$sim_name)  # See scenarios.R
      
    } else {
      
      # Load strategy details from file
      this_strategy = readRDS(paste0(o$pth$strategy, "strategy.", this_analysis$strategy, ".rds"))
      
      # Select the combination we are interesting in here
      this_combination = this_strategy$combinations %>% 
        filter(id == this_analysis$sim_name)
      
      # Throw an error if we have anything other than one row here
      if (nrow(this_combination) != 1)
        stop("Strategy ID '", this_analysis$sim_name, "' not recognised")
      
      # Alter these parameters based on strategy (see strategies.R)
      p = parse_strategy(o, p, this_strategy$details, this_combination)
    }
    
    # Simulate the model
    m = model(o, p, seed = this_seed, verbose = "date")  # See model.R
    
    # Prepare output to be stored
    sim_result = list(canton   = this_analysis$canton, 
                      sim_name = this_analysis$sim_name, 
                      params   = this_params,
                      seed     = this_seed, 
                      input    = p, 
                      output   = m)
    
    # Save model outcomes as an RDS file
    saveRDS(sim_result, paste0(o$pth$sims, this_analysis$sim_id, ".rds"))
  }
  
  # ---- Summarise analysis simulations ----
  
  if (grepl("summarise", sim_type)) {
    
    # Load full factorial set of cantons, seeds, and simulation names
    analysis_df = readRDS(paste0(o$pth$sims, "all_simulations.rds"))
    
    # We'll summarise across parameter sets and seeds, so interested in cantons and simulation names
    summarise_df = unique(select(analysis_df, canton, sim_name))
    
    # Select analyses to summarise based on task ID
    this_canton   = summarise_df$canton[task_id]
    this_sim_name = summarise_df$sim_name[task_id]
    
    # All simulations to summarise over
    this_sim_df = analysis_df %>%
      filter(canton   == this_canton,
             sim_name == this_sim_name)
    
    # Initiate analysis list
    a = list(analysis_name = o$analysis_name,
             sim_name      = this_sim_name,
             time_stamp    = format(Sys.time(), "%Y%m%d_%H%M"))
    
    # Preallocate list for model outcomes
    m_list = list()
    
    # TODO: sim_result$model_input does not change over seeds... append it here?
    
    # Loop through outputs associated with this simulation
    for (i in 1 : nrow(this_sim_df)) {
      sim_id = this_sim_df$sim_id[i]
      
      # Load model outcomes
      sim_result = readRDS(paste0(o$pth$sims, sim_id, ".rds"))
      
      # Convert to datatable and append sample number
      m_list[[i]] = sim_result$output %>%
        mutate(sim_number = i)
    }
    
    # Create key summary statistics across seeds
    a$model_output = rbindlist(m_list) %>% 
      group_by(date, metric, grouping, group) %>% 
      summarise(mean   = mean(value),
                median = quantile(value, 0.5),
                lower  = quantile(value, o$quantiles[1]),
                upper  = quantile(value, o$quantiles[2])) %>% 
      mutate(scenario = this_sim_name) %>% 
      as.data.table()
    
    # Calculate and append likelihood score for this canton for baseline scenario
    if (this_sim_name == "baseline") {
      
      # TODO: Likelihood value associated with this parameter set (see likleihood.R)
      # a$likelihood = likelihood_model(o, m = a$best, data_from_cache = FALSE)
    }
    
    # Load info from fitting file - includes data and options used to generate fit
    data_list = readRDS(paste0(o$pth$analysis, this_canton, "_data.rds"))
    fit_inherit = c("fit_data", "fit_info", "fit_options")
    
    # Store a few options used to generate the simulations and inherit from fitting process
    a$options = store_options(o, "analysis")
    a[fit_inherit] = data_list[fit_inherit]
    
    # Store summarised result
    saveRDS(a, file = paste0(o$pth$scenario, this_canton, "_", this_sim_name, ".rds"))
  }
  
  # ---- Run test simulations: test stochasticity -----
  
  # Testing stochasticity: Run the same parameter set for a load of seeds
  if (grepl("test_stochasticity", sim_type)) {
    
    # Load full parameter set dataframe from file
    param_sets = readRDS(paste0(o$pth$testing, "stochasticity_samples.rds"))
    
    # Load recently cached data (see load_data.R)
    d = load_data(o, from_cache = TRUE)
    
    # Generate model parameters using prior means for non-fixed parameters
    p = get_parameters(o, d[[o$cantons]], calibration_params = NULL)  # See parameters.R
    
    # Extract simulation ID and seed number
    sim_id   = param_sets$sim_id[task_id]
    seed_num = param_sets$seed[task_id]
    
    # Run model (see model.R)
    m = model(o, p, seed = seed_num)
    
    # Save model outcomes as an RDS file
    saveRDS(m, paste0(o$pth$testing, sim_id, ".rds"))
  }
  
  # ---- Run test simulations: simulations with 'high' likelihood -----
  
  # Testing likelihood: Run model for sims with high likelihood
  if (grepl("test_likelihood", sim_type)) {
    
    # Load full parameter set dataframe from file
    param_sets = readRDS(paste0(o$pth$testing, "likelihood_samples.rds"))
    
    # Select the parameter set defined by task_id
    param_values = param_sets$samples %>% 
      select(-sim_id, -seed, -likelihood) %>% 
      dplyr::slice(task_id) %>% 
      unlist()
    
    # Re-transform this parameter set to real scale
    param_set = param_values * param_sets$scale + param_sets$center
    
    # Load recently cached data (see load_data.R)
    d = load_data(o, from_cache = TRUE)
    
    # Generate model parameters using this parameter set
    p = get_parameters(o, d[[o$cantons]], param_set)  # See parameters.R
    
    # Extract simulation ID and seed number
    sim_id   = param_sets$samples$sim_id[task_id]
    seed_num = param_sets$samples$seed[task_id]
    
    # Run model (see model.R)
    m = model(o, p, seed = seed_num)
    
    # Save model outcomes as an RDS file
    saveRDS(m, paste0(o$pth$testing, sim_id, ".rds"))
  }
}

