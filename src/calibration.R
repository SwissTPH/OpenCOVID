###########################################################
# CALIBRATION
#
# Fit selection of model parameters to epidemioligical data. 
# Operates on a trained emulator rather than directly calling
# the model and calculating likelihoods.
#
# See here for a simple example of the MCMC algorithm, lazyMCMC:
# https://jameshay218.github.io/lazymcmc/inst/doc/sir_example.html
#
###########################################################

# ---------------------------------------------------------
# Parent function to fit model to epidemiological data
# ---------------------------------------------------------
run_calibration = function(o) {
  
  # Only continue if specified by do_step
  if (!is.element(2, o$do_step)) return()
  
  message("* Calibrating model")
  
  # ---- Load model emulator ----
  
  # Load the fitted emulator - throw an error if it doesn't exist explaining step 1 needs to run first
  err_msg = "You must train a model emulator (step = 1) before calibrating (step = 2)"
  emulator = try_load(o$pth$emulator, "model_emulator", msg = err_msg) 
  
  # Extract data used to train the emulator
  d = emulator$data

  # Append model to o list for use in MCMC likelihood function
  o$emulator = emulator[c("model", "center", "scale", "calibration_df")]
  
  # Calibration dataframe - transform parameter bounds (do initial point when initialising MCMC)
  param_df = o$calibration_df %>% # emulator$calibration_df %>%
    mutate(steps = o$mcmc_init_steps,  # Also add step sizes before feeding into lazyMCMC
           lower_bound = (lower_bound - o$emulator$center) / o$emulator$scale, 
           upper_bound = (upper_bound - o$emulator$center) / o$emulator$scale, 
           prior_sd = prior_sd / o$emulator$scale)
  
  # ---- Load all simulated parameter sets ----
  
  # NOTE: We also store this in canton fit as we may use this as the 'best estimate' 
  #       simulation and use the posteriors to generate calibrated parameter uncertainty
  
  # Load all parameter sets we've actually simulated
  all_simulated = try_load(o$pth$emulator, "all_samples")
  
  # Get best parameter set (see o$best_parameter_by)
  best_sample = get_best_simulated(o, all_simulated$samples, n = 1)
  
  # Transform back to orginal scale
  best_simulated = best_sample$params * all_simulated$scale + all_simulated$center
  
  # ---- Display size of this fitting job ----
  
  # Total number of parameters
  n_params = nrow(param_df)
  
  message(" - Calibrating ", n_params, " parameters in total:")
  
  # Number of parameters global and canton-specific
  n_global = sum(param_df$scope == "global")
  n_hyper  = sum(param_df$scope == "hyper")
  n_local  = n_params - n_global - n_hyper
  
  message("  > ", n_global, " global parameters")
  message("  > ", n_local,  " canton-specific parameters")
  message("  > ", n_hyper,  " hyper parameters")
  
  # ---- Perform MCMC multiple times ----
  
  # TODO: Get this running on the cluster
  
  # Run the fitting proccess for mcmc_n_chains number of chains
  for (i in 1 : o$mcmc_n_chains)
    run_mcmc_chain(o, param_df, n_params, best_simulated, chain_id = i)
  
  # Preallocate for all chains
  all_chains = NULL
  
  # Loop through files we've just generated
  for (i in 1 : o$mcmc_n_chains) {
    
    # Extract MCMC chain
    this_chain = try_load(o$pth$chains, paste0("chain_", i), 
                          msg = "MCMC chain not found") 
    
    # Continuously bind to all others
    all_chains = rbind(all_chains, this_chain)
  }
  
  # Save chains and other bits in an RDS file
  saveRDS(all_chains, paste0(o$pth$fitting, "all_chains.rds"))
  
  # ---- Store canton specific results ----
  
  # NOTE: We split up by canton incase any cantons need to be rerun
  
  # We'll use time stamps to differentiate analyses
  time_stamp = format(Sys.time(), "%Y%m%d_%H%M")
  
  # Loop through cantons
  for (canton in o$cantons) {
    
    # Extract parameter names for this canton
    canton_idx    = emulator$calibration_df$scope %in% c("global", canton)
    canton_params = emulator$calibration_df$names[canton_idx]
    
    # Extract best simulated parameter values for this canton
    canton_best_simulated = setNames(best_simulated[canton_idx], canton_params)
    
    # Subset the chain columns
    canton_chain = all_chains %>% 
      select(sampno, one_of(canton_params), lnlike)
    
    # File paths to store this canton result - these are what are used in run_analysis
    canton_file    = paste0(o$pth$fitting, canton, "_calibration.rds")
    canton_archive = paste0(o$pth$archive, canton, "_calibration_", time_stamp, ".rds")
    
    # Append chain and the input options used to create this fit
    fit = list(fit_name       = o$calibration_name, 
               fit_cantons    = o$cantons, 
               time_stamp     = time_stamp, 
               chain          = canton_chain, 
               best_simulated = canton_best_simulated, 
               best_sim_seed  = best_sample$seed, 
               param_table    = o$parameter_table, 
               data           = d[[canton]], 
               known_params   = d$known_params, 
               options        = store_options(o, "calibration"))
    
    # Save the fit object
    saveRDS(fit, canton_file)
    
    # Also store this in archive (with time stamp)
    file.copy(from = canton_file, to = canton_archive)
  }
  
  # ---- Produce 'quick' calibration plot ----
  
  # Load file of best fitting simulation
  best_pth  = paste0(o$pth$samples, best_sample$sample_id, file_sep())
  best_file = readRDS(paste0(best_pth, best_sample$sim_id, ".rds"))

  # Collate results - use trivial bounds as we have only a single simulation
  model_output = best_file$model_output %>%
    select(date, metric, grouping, group,
           mean  = value, median = value,
           lower = value, upper  = value) %>%
    mutate(scenario = "model_test")

  # Append to an analysis file along with data
  a = list(fit_data     = fit$data, 
           model_input  = best_file$model_input, 
           model_output = model_output)

  # Names of metrics we calibrate to - can only plot these as that's all we store when calibrating
  calibration_metrics = names(o$calibration_multiplier[-1])

  # Plot a small selection of grouped metrics to quickly visualise calibration
  fig_name = c("Quick calibration plot ", canton)
  plot_temporal(o, canton, fig_name, plot_file = a,
                likelihood   = best_file$likelihood,
                plot_metrics = calibration_metrics,
                plot_to      = max(o$dates_data))
}

# ---------------------------------------------------------
# Call MCMC algorithm for given chain ID
# ---------------------------------------------------------
run_mcmc_chain = function(o, param_df, n_params, best_simulated, chain_id = NA) {
  
  message(" - Generating chain ", chain_id)
  
  # Set the random number seed generator if desired
  if (o$calibration_reproducible == TRUE)
    set.seed(chain_id)
  
  # File path stem for everything that'll be saved for this chain
  save_path = paste0(o$pth$chains, "chain_", chain_id)
  
  # LazyMCMC workaround to be able to provide the likelihood function with additional arguments
  data = o[qc(emulator, prior2data_weights)]
  
  # ---- Initial point ----
  
  # For the first chain use an initial condition we expect to get to the global minimum
  if (chain_id == 1) {
    
    # TODO: Perhaps best default would be emulator best estimate rather than prior means 
    # param_df$values = (as.vector(unname(best_simulated)) - o$emulator$center) / o$emulator$scale

    # Transform parameter prior means using standardisation from emulator samples
    param_df$values = (param_df$values - o$emulator$center) / o$emulator$scale

  } else {  # For other chains generate a random initial point in parameter space
    param_df$values = runif(nrow(param_df), min = param_df$lower_bound, max = param_df$upper_bound)
  }
  
  # ---- Univaritae fit ----
  
  # Sensible adaptive period - just use half for univariate
  uvr_adaptive_period = o$mcmc_opt_stages * o$mcmc_opt_freq / 2
  
  # Set a few options for the lazy MCMC algorithm
  uvr_options = c("iterations" = o$mcmc_iters / 2, # Post warm-up iterations - just use half for univariate
                  "popt"       = 0.44,
                  "opt_freq"   = o$mcmc_opt_freq, 
                  "thin"       = 1,
                  "save_block" = 100, 
                  "adaptive_period" = uvr_adaptive_period)
  
  message("  > Running univariate MCMC")
  
  # Run the algorithm (output saved to file)
  run_MCMC(parTab   = as.data.frame(param_df),
           data     = data,
           mcmcPars = uvr_options,
           mvrPars  = NULL,
           filename = save_path,
           PRIOR_FUNC = prior_function,
           CREATE_POSTERIOR_FUNC = likelihood_function)
  
  # Load the resulting univariate MCMC calibration output file
  uvr_file  = paste0(save_path, "_univariate_chain.csv")
  uvr_chain = read.csv(uvr_file) %>% 
    filter(sampno > max(uvr_adaptive_period, o$min_cut))
  
  # ---- Multivariate fit ----
  
  # No need to run multivariate if fitting only one parameter
  if (n_params == 1) {
    mvr_chain = uvr_chain
    
  } else {
    
    # Take as starting value the optimum found before
    param_df$values = get_best_pars(uvr_chain)
    
    # Avoid rounding errors and ensure we've stayed within bounds
    param_df$values = pmin(param_df$values, param_df$upper_bound)
    param_df$values = pmax(param_df$values, param_df$lower_bound)
    
    # Check how parameter values from univariate chain are correlated
    cov_mat = cov(uvr_chain[, param_df$names])  # Need this for the multivariate fit
    
    # Set additional parameters for the multivariate search, in which a
    # covariate matrix is learned from the fit, but initialised with cov_mat
    mvr_input = list(cov_mat, 2.38 / sqrt(n_params), w = 0.8)
    
    # Sensible adaptive period
    mvr_adaptive_period = o$mcmc_opt_stages * o$mcmc_opt_freq
    
    # Set a few options for the lazy MCMC algorithm
    mvr_options = c("iterations" = o$mcmc_iters,  # Post warm-up iterations
                    "popt"       = 0.234,  # Optimal acceptance rate is different for multivariate proposals
                    "opt_freq"   = o$mcmc_opt_freq, 
                    "thin"       = round(max(1, o$mcmc_iters / 1000)),  # 1000 post draws
                    "save_block" = 100, 
                    "adaptive_period" = mvr_adaptive_period)
    
    message("  > Running multivariate MCMC")
    
    # Now run again, this time multivariate fit
    run_MCMC(parTab   = as.data.frame(param_df),
             data     = data,
             mcmcPars = mvr_options,
             mvrPars  = mvr_input,
             filename = save_path,
             PRIOR_FUNC = prior_function,
             CREATE_POSTERIOR_FUNC = likelihood_function,
             OPT_TUNING = 0.2)
    
    # Load the resulting multivariate MCMC calibration output file
    mvr_file  = paste0(save_path, "_multivariate_chain.csv")
    mvr_chain = read.csv(mvr_file) %>% 
      filter(sampno > max(mvr_adaptive_period, o$min_cut))
  }
  
  # Transform all values back to real scale (trivial transformation for sampno and lnlike)
  real_chain = t(t(as.matrix(mvr_chain)) * 
                   c(1, o$emulator$scale, 1) + 
                   c(0, o$emulator$center, 0))
  
  # Append column to identify the different chains and reset sample numbers
  chain_df = as.data.frame(real_chain) %>%
    mutate(sampno   = 1 : nrow(real_chain), 
           chain_id = chain_id)
  
  # Save this to file
  saveRDS(chain_df, paste0(o$pth$chains, "chain_", chain_id, ".rds"))
}

# ---------------------------------------------------------
# Likelihood of parameters given the data using model emulator
# ---------------------------------------------------------
likelihood_function = function(parTab, data, PRIOR_FUNC) {
  
  # ---------------------------------------------------------
  # Likelihood function needs to be in a closure environment
  # ---------------------------------------------------------
  likelihood_emulator = function(param_set) {
    
    # Use model emulator to predict likelihood for this parameter set
    likelihood = stats::predict(x = matrix(param_set, nrow = 1), 
                                object = data$emulator$model)$mean
    
    # Evaluate parameter priors and increment
    if (!is.null(prior_function))
      likelihood = likelihood + prior_function(parTab, param_set)
    
    return(-likelihood)
  }
}

# ---------------------------------------------------------
# Define parameter priors
# ---------------------------------------------------------
prior_function = function(param_df, param_set) {

  # Initiate 
  prior_vals = 0
  
  # ---- All priors inc. is_hyper excl. has_hyper ----
  
  # Loop through parameters we're calibrating
  for (i in 1 : nrow(param_df) ) {
    param = param_df[i, ]
    
    # Extract prior distribution function
    prior_val = dnorm(param_set[i], param$values, param$prior_sd, log = TRUE)
    
    # Diffuse prior distribution for this parameter
    prior_vals = prior_vals + prior_val * param$prior_weight
  }
  
  return(prior_vals)
}

# ---------------------------------------------------------
# Get best parameter set from those simulated
# ---------------------------------------------------------
get_best_simulated = function(o, samples, n = 1) {
  
  # Options for o$best_parameter_by:
  #  "mean"    := The best mean likelihood across all seeds - any sets with NAs are discarded
  #  "mean_na" := As above but NAs are just ignored and can still be selected
  #  "max"     := The very best seed simulation (not so robust but occasionally produces better results)
  
  # The best n of all seed simulations (in terms of associated likelihood value)
  if (o$best_parameter_by == "max") {
    best_param_id = unique(arrange(samples, likelihood)$param_id)[1 : n]
    
  } else { # Otherwise we'll take the mean across seeds...
    
    # Calculate mean across all seeds for each parameter set
    if (o$best_parameter_by == "mean") {
      mean_df = group_by(samples, param_id) %>%
        summarise(mean_likelihood = mean(likelihood))
    }
    
    # As above, but do not discard sims with NA (just don't count them)
    if (o$best_parameter_by == "mean_na") {
      mean_df = group_by(samples, param_id) %>%
        summarise(mean_likelihood = mean(likelihood, na.rm = TRUE))
    }
    
    # Extract IDs of the best n parameter sets
    best_param_id = slice_min(mean_df, mean_likelihood, n = n)$param_id
  }
  
  # Throw an error if we haven't picked up a best parameter set
  if (!exists("best_param_id"))
    stop("Error in extracting best simulated parameter set - check value of 'best_parameter_by'")

  # Take the best simulations from identified parameter sets
  best_params_df = samples %>%
    filter(param_id %in% best_param_id) %>%
    group_by(param_id) %>%
    slice_min(likelihood , n = 1) %>%
    as.data.table()
  
  # Extract parameter values for these best parameter set
  best_params = best_params_df %>%
    select(-sim_id, -sample_id, -param_id, -seed, -likelihood)
  
  # Store parameter sets and seeds in output
  best_sample = list(params    = as.matrix(best_params), 
                     seed      = best_params_df$seed, 
                     sample_id = best_params_df$sample_id , 
                     sim_id    = best_params_df$sim_id)
  
  return(best_sample)
}

