###########################################################
# CALIBRATION
#
# Fit parameters to achieve a given R_eff at the 'initial'
# point for forward projections.
#
###########################################################

# ---------------------------------------------------------
# Parent function to fit contacts to desired R_eff
# ---------------------------------------------------------
run_calibration = function(o) {
  
  # Only continue if specified by do_step
  if (!is.element(1, o$do_step)) return()
  
  message("* Calibrating model")
  
  # If we want this process to be reproducible, set a seed 
  if (o$fit_reproducible == TRUE)
    set.seed(1)
  
  # ---- Parameters to be fitted ----
  
  message(" - Sampling parameters")
  
  # Load parsed parameters from yaml file
  fit_input = parse_yaml(o, scenario = "baseline")$parsed
  
  # Extract parameters to calibrate
  fit_df = list2dt(fit_input$calibration) %>%
    mutate(lower = as.numeric(lower), 
           upper = as.numeric(upper))
  
  # Shorthand for fitted parameters
  fit_params = fit_df$param
  
  # Bounds of all parametrs to be fitted, as matrix
  fit_bounds = as.matrix(fit_df[, -1])
  
  # Sample and expand each parameter set for each seed
  samples_df = lhs(o$emulator_samples, fit_bounds) %>%
    as_named_dt(fit_params) %>%
    cbind(param_id = 1 : o$emulator_samples) %>%
    expand_grid(seed = 1 : o$emulator_seeds) %>%
    select(param_id, seed, all_of(fit_params)) %>%
    as.data.table()
  
  # Save simulation df for reference on cluster
  saveRDS(samples_df, file = paste0(o$pth$fit_samples, "all_samples.rds"))
  
  # ---- Simulate parameter samples ----
  
  message(" - Simulating samples")

  # Submit all jobs to the cluster (see auxiliary.R)
  submit_cluster_jobs(o, o$emulator_samples, "bash_submit.sh", "fitting")

  # Throw an error if any cluster jobs failed (see auxiliary.R)
  err_tol = floor(o$emulator_samples * o$sample_err_tol)
  stop_if_errors(o$pth$log, o$err_file, err_tol = err_tol)

  # Remove all log files if desired (generally a good idea unless debugging)
  if (o$rm_cluster_log) unlink(paste0(o$pth$log, "*"), force = TRUE)
  
  # Output files produced
  output_files = list.files(path    = o$pth$fit_samples,
                            pattern = "sample_[0-9]+.rds",  
                            full.names = TRUE)

  # Load these files into one long datatable and combine with samples
  samples_df = rbindlist(lapply(output_files, readRDS)) %>%
    right_join(samples_df, by = c("param_id", "seed"))
  
  # Overwrite samples file for plotting purposes
  saveRDS(samples_df, file = paste0(o$pth$fit_samples, "all_samples.rds"))
  
  # ---- Train model emulator ----
  
  message(" - Training model emulator")
  
  # First job is to normalise response variable
  norm_df = samples_df # %>%
    # mutate(r_eff = normalise_0to1(r_eff))
  
  # We'll then do the same for all predictor variables
  for (param in fit_params) {
    param_idx = which(fit_params %in% param)
    
    # Use parameter bounds rather than bounds of sampled points
    norm_df[[param]] = normalise_0to1(
      x = norm_df[[param]], 
      x_min = fit_bounds[[param_idx, 1]], 
      x_max = fit_bounds[[param_idx, 2]])
  }
  
  # Sample some parameter sets to holdback
  n_test  = round(o$emulator_samples * o$test_train_split)
  test_id = sort(sample(1 : o$emulator_samples, n_test))
  
  # Split samples into train and test - summarise test outcomes
  train_data = filter(norm_df, !param_id %in% test_id)
  test_data  = filter(norm_df, param_id %in% test_id) %>%
    group_by(param_id, contacts) %>%
    summarise(r_eff = mean(r_eff)) %>%
    as.data.table()
  
  # Seperate inputs and output in the training data
  train_variables = as.matrix(train_data[, ..fit_params])
  train_response  = as.matrix(train_data[, r_eff])
  
  # Prepare data for mleHetGP: find and remove duplicates
  prep_data = find_reps(X = train_variables, Z = train_response)
  gp_input = list(X0 = as.matrix(prep_data$X0), 
                  Z0 = as.matrix(prep_data$Z0),
                  mult = prep_data$mult)
  
  # Run GP model (mleHetGP from hetGP package)
  emulator = mleHetGP(X = gp_input, 
                      Z = prep_data$Z, 
                      covtype = o$gp_kernel,
                      maxit   = o$gp_max_iter)
  
  # ---- Test performance of emulator ----
  
  message(" - Assessing emulator performance")
  
  # Split test independent and dependent variables
  test_variables = as.matrix(test_data[, ..fit_params])
  test_response  = as.vector(test_data[, r_eff])
  
  # Evaluate test set using model emulator
  test_predict = predict(emulator, x = test_variables)$mean
  
  # Compile actual and predicted points for test set
  test_df = data.table(actual  = test_response, 
                       predict = test_predict, 
                       group   = "test")
  
  # Also predict outcomes for points in the train set - just for plotting
  train_predict = predict(emulator, x = train_variables)$mean
  
  # Concatenate actual and predicted for points in the train set
  gp_performance = 
    data.table(actual  = as.vector(train_response), 
               predict = train_predict, 
               group   = "train") %>%
    rbind(test_df)
  
  # Append full performance info to emualator list
  emulator$performance = gp_performance
  
  # Save emulator to file
  saveRDS(emulator, paste0(o$pth$fitting, "model_emulator.rds"))
  
  # ---- Optimise parametrs for R_eff ----
  
  message(" - Optimising for target R_eff")
  
  # Preallocate output for each optimisation run
  fit_output = list(x = matrix(ncol = length(fit_params),  
                               nrow = o$optim_runs), 
                    y = matrix(ncol = o$optim_runs, 
                               nrow = o$fit_iters_max + 1))
  
  # Provide full o list as additonal argument
  obj_args = list(emulator = emulator, 
                  bounds   = fit_bounds, 
                  target   = fit_input$r_eff)
  
  # Generate set of starting points
  x0_mat = lhs(o$optim_runs, fit_bounds)
  
  # Initiate a progress bar
  pb = start_progress_bar(o$optim_runs)
  
  # Repeat optimisation optim_runs number of times
  for (i in 1 : o$optim_runs) {
    
    # Run ASD algorithm to determine optimal contacts
    this_optim = asd(objective_fn, 
                    x0   = x0_mat[i, ],
                    args = obj_args,
                    lb   = fit_df$lower,
                    ub   = fit_df$upper,
                    max_iters  = o$fit_iters_max,
                    plot_iters = NULL, # o$fit_iters_plot, 
                    verbose = FALSE)
    
    # Store result
    fit_output$x[i, ] = this_optim$x
    fit_output$y[, i] = this_optim$y_vec
    
    # Update progress bar
    setTxtProgressBar(pb, i)
  }
  
  # Close progress bar
  close(pb)
  
  # Optimisation producing best optimal
  best_idx = which.min(colMins(fit_output$y))
  best_val = fit_output$x[best_idx]
  
  # Compile optimal result into named list
  fit_list = as.list(setNames(best_val, fit_params))

  # Store result, input, and optim output in list
  fit_result = list(params = fit_params, 
                    bounds = fit_bounds, 
                    result = fit_list, 
                    input  = fit_input,
                    output = fit_output)

  # Save all these details to file
  saveRDS(fit_result, file = paste0(o$pth$fitting, "fit_result.rds"))
  
  # ---- Performance plots ----
  
  message(" - Assessing calibration performance")
  
  # Plot emulator performance (see plotting.R)
  plot_emulator(o)
  
  # Plot optimisation performance (see plotting.R)
  plot_optimisation(o)
}

# ---------------------------------------------------------
# Objective function to minimise - modelled and desired R_eff difference
# ---------------------------------------------------------
objective_fn = function(x, args) {
  
  # Normalise x for use with emulator
  x_norm = normalise_0to1(
    x = x, 
    x_min = args$bounds[1], 
    x_max = args$bounds[2]) %>% 
    as.numeric()
  
  # Evaluate emulator at this point
  reff_value = predict(args$emulator, x = x_norm)$mean
  
  # Squared difference with target R_eff
  objective_val = (reff_value - args$target) ^ 2
  
  # Output in list format
  return(list(y = objective_val))
}

# ---------------------------------------------------------
# Alter key parameters in yaml file for fitting simulation
# ---------------------------------------------------------
fit_yaml = function(o, yaml, fit_list) {
  
  # Only needed if fitting OR using fitted parameters
  if (!is.null(fit_list)) {
    
    # Check if we are fitting here - number of days will be specified
    if (!is.null(fit_list$n_days)) {
      
      # Throw error if R_eff is not being reported - we need this to calibrate to
      if (!"R_effective" %in% yaml$parsed$metrics$df$metric)
        stop("You must turn on reporting of 'R_effective' for calibration")
      
      # Only interested in collecting one metric when fitting: R_effective
      yaml$parsed$metrics$df = filter(yaml$parsed$metrics$df, metric == "R_effective")
      yaml$parsed$metrics$groupings = yaml$parsed$metrics$groupings["R_effective"]
      
      # Extend the simulation period so we collect over time period of interest
      fit_list$n_days = fit_list$n_days + yaml$parsed$r_eff_window - 1
      
      # Sanity check that n_days from yaml file is at least max of fit_days
      if (yaml$parsed$n_days < fit_list$n_days)
        stop("Model calibration requires n_days >= ", fit_list$n_days)
    }
    
    # Loop through parameters to alter
    for (param in names(fit_list)) {
      
      # Check parameter is valid - throw error if not
      if (!param %in% names(yaml$parsed))
        stop("Attempting to alter invalid parameter '", param, "' during model fitting")
      
      # Re-define paramter value
      yaml$parsed[[param]] = fit_list[[param]]
    }
  }
  
  return(yaml)
}

