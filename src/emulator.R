###########################################################
# EMULATOR
#
# Train a model emulator on simulated samples, with a response
# variable of normalised error between model output and fitting
# target (which could be a single Re value or temporal epi data). 
#
# Uses a Guassian Proccess model to interpolate parameter space.
# An optimisation algorithm is then used to locate the global
# minimum of the emulated space. 
#
###########################################################

# TODO: The transformation functions aren't so hot - tidy up required

# ---------------------------------------------------------
# Train a model emulator on simulated samples
# ---------------------------------------------------------
train_emulator = function(o, fit, r_idx) {
  
  message("  > Training model emulator")
  
  # Shorthand of emulator options
  opts   = fit$input$emulator
  params = fit$params
  
  # Load samples for this and all preceeding sampling rounds
  samples_df = try_load(o$pth$fitting, "rx_samples") %>% 
    normalise_params(fit$params, fit$bounds) %>%  # Transform to unit cube
    mutate(paramset_id = get_paramset_id(param_id), .before = 1)
  
  # Extract sample IDs (param_id without the seed reference)
  paramset_ids = unique(samples_df$paramset_id)
  
  # Sample some parameter sets to holdback
  n_test  = round(length(paramset_ids) * opts$test_train_split)
  test_id = sort(sample(paramset_ids, size = n_test))
  
  # Split samples into train and test - summarise test outcomes
  train_data = samples_df[!paramset_id %in% test_id, ]
  test_data  = samples_df[paramset_id %in% test_id, ] %>%
    group_by_at(c("paramset_id", params)) %>%
    summarise(obj_value = mean(obj_value)) %>%
    ungroup() %>%
    setDT()
  
  # ---- Train model emulator ----
  
  # Seperate inputs and output in the training data
  train_variables = as.matrix(train_data[, ..params])
  train_response  = as.matrix(train_data[, obj_value])
  
  # Prepare data for mleHetGP: find and remove duplicates
  prep_data = find_reps(X = train_variables, Z = train_response)
  gp_input = list(X0 = as.matrix(prep_data$X0), 
                  Z0 = as.matrix(prep_data$Z0),
                  mult = prep_data$mult)
  
  # Run GP model (mleHetGP from hetGP package)
  emulator = mleHetGP(X = gp_input, 
                      Z = prep_data$Z, 
                      covtype = opts$gp_kernel,
                      maxit   = opts$gp_max_iter)
  
  # ---- Test performance of emulator ----
  
  # Split test independent and dependent variables
  test_variables = as.matrix(test_data[, ..params])
  test_response  = as.vector(test_data[, obj_value])
  
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
               group   = paste0("train::r", train_data$round)) %>%
    bind_rows(test_df)
  
  # Append full performance info to emualator list
  emulator$performance = gp_performance
  
  # Save emulator to file
  saveRDS(emulator, paste0(o$pth$fitting, r_idx, "_emulator.rds"))
}

# ---------------------------------------------------------
# Find global minimum of emulated space
# ---------------------------------------------------------
search_emulator = function(o, fit, r_idx) {
  
  message("  > Searching emulated space")
  
  # Shorthand of optimisation options
  opts = fit$input$optimisation
  
  # Load fitted emulator
  emulator = try_load(o$pth$fitting, paste0(r_idx, "_emulator")) 
  
  # ---- Optimise parameters for best quality of fit ----
  
  # Preallocate output for each optimisation run
  fit$optim = list(x = matrix(ncol = length(fit$params),  
                              nrow = opts$optim_runs), 
                   y = matrix(ncol = opts$optim_runs, 
                              nrow = opts$optim_iters + 1))
  
  # Provide list as additonal argument for optimisation
  obj_args = list(emulator = emulator, 
                  params   = fit$params, 
                  bounds   = fit$bounds)
  
  # Generate set of starting points
  x0_mat = lhs(opts$optim_runs, fit$bounds)
  
  # Initiate a progress bar
  pb = start_progress_bar(opts$optim_runs)
  
  # Repeat optimisation optim_runs number of times
  for (i in 1 : opts$optim_runs) {
    
    # Run ASD algorithm to determine optimal contacts
    this_optim = asd(evaluate_emulator, 
                     x0   = x0_mat[i, ],
                     args = obj_args,
                     lb   = fit$bounds[, 1],
                     ub   = fit$bounds[, 2],
                     max_iters  = opts$optim_iters,
                     plot_iters = NULL,
                     verbose    = FALSE)
    
    # Store result
    fit$optim$x[i, ] = this_optim$x
    fit$optim$y[, i] = this_optim$y_vec
    
    # Update progress bar
    setTxtProgressBar(pb, i)
  }
  
  # Close progress bar
  close(pb)
  
  # Optimisation producing best optimal
  best_idx = which.min(colMins(fit$optim$y))
  best_val = fit$optim$x[best_idx, ]
  
  # Compile optimal result into named list
  fit$result = as.list(setNames(best_val, fit$params))
  
  # Save fitted result
  saveRDS(fit, file = paste0(o$pth$fitting, r_idx, "_result.rds"))
}

# ---------------------------------------------------------
# Evaluate ML model for a given set of predictor points
# ---------------------------------------------------------
evaluate_emulator = function(x, args) {
  
  # Transform parameters ready for emulator evaluation
  x_norm = x %>% 
    setNames(args$params) %>% 
    normalise_params(args$params, args$bounds) %>% 
    as.numeric()
  
  # Evaluate emulator at this point
  obj_value = predict(args$emulator, x = x_norm)$mean
  
  # Output in list format
  return(list(y = obj_value))
}

# ---------------------------------------------------------
# Train parameter input -> model output
# ---------------------------------------------------------
train_predictor = function(o, data, info, opts) {
  
  # Shorthand for predictor variables
  params = info$vars$variable_id
  
  # Transform all variables by standardising
  #
  # NOTE: Avoiding normalising as no natural bounds for response variable
  list[transform_df, transform_details] = 
    standardise_params(data, params, "value")
  
  # Sample some scernarios to holdback
  n_test  = round(length(unique(data$scenario)) * opts$test_train_split)
  test_id = sort(sample(unique(data$scenario), size = n_test))
  
  # Split samples into train and test - summarise test outcomes
  train_data = transform_df[!scenario %in% test_id, ]
  test_data  = transform_df[scenario %in% test_id, ] %>%
    group_by_at(c("scenario", params)) %>%
    summarise(value = mean(value)) %>%
    ungroup() %>%
    setDT()
  
  # ---- Train outcome predictor ----
  
  # Seperate inputs and output in the training data
  train_variables = as.matrix(train_data[, ..params])
  train_response  = as.matrix(train_data[, value])
  
  # Prepare data for mleHetGP: find and remove duplicates
  prep_data = find_reps(X = train_variables, Z = train_response)
  gp_input = list(X0 = as.matrix(prep_data$X0), 
                  Z0 = as.matrix(prep_data$Z0),
                  mult = prep_data$mult)
  
  # Run GP model (mleHetGP from hetGP package)
  gp_model = mleHetGP(X = gp_input, 
                      Z = prep_data$Z, 
                      covtype = opts$gp_kernel,
                      maxit   = opts$gp_max_iter)
  
  # ---- Test performance of predictor ----
  
  # Split test independent and dependent variables
  test_variables = as.matrix(test_data[, ..params])
  test_response  = as.vector(test_data[, value])
  
  # Evaluate test set using model predictor
  test_predict = predict(gp_model, x = test_variables)$mean
  
  # Compile actual and predicted points for test set
  test_df = data.table(actual  = as.vector(test_response), 
                       predict = test_predict, 
                       group   = "test")
  
  # Also predict outcomes for points in the train set - just for plotting
  train_predict = predict(gp_model, x = train_variables)$mean
  
  # Concatenate actual and predicted for points in the train set
  gp_performance = 
    data.table(actual  = as.vector(train_response), 
               predict = train_predict, 
               group   = "train") %>%
    bind_rows(test_df)
  
  # Output data, transformation, model, and performance details
  predictor = list(data        = data, 
                   transform   = transform_details, 
                   model       = gp_model, 
                   performance = gp_performance)
  
  return(predictor)
}

# ---------------------------------------------------------
# Evaluate ML model for a given set of predictor points
# ---------------------------------------------------------
evaluate_predictor = function(data, predictor, p_int = 0.95) {
  
  # Loop through variables
  for (var in names(data)) {
    
    # Extract transformation details
    tf = predictor$transform[p == var, ]
    
    # Apply this transformation
    data[[var]] = (data[[var]] - tf$mean) / tf$sd
  }
  
  # Evaluate emulator at this point
  p = predict(predictor$model, x = as.matrix(data))
  
  # The mean represents the best estimate
  best = p$mean
  
  # Generate confidence (or prediction) intervals
  lb = qnorm(1-p_int, p$mean, sqrt(p$sd2)) # + p$nugs))
  ub = qnorm(p_int,   p$mean, sqrt(p$sd2)) # + p$nugs))
  
  # Extract output transformation details
  tf = predictor$transform[p == "value", ]
  
  # Apply inverse of transformation to backtransform
  result = data.table(
    best = (best * tf$sd) + tf$mean, 
    lb   = (lb   * tf$sd) + tf$mean, 
    ub   = (ub   * tf$sd) + tf$mean)
  
  return(result)
}

# ---------------------------------------------------------
# Transform parameters by normalising to a unit cube
# ---------------------------------------------------------
normalise_params = function(df, params, bounds) {
  
  # Normalise all predictor variables
  for (param in params) {
    param_idx = which(params %in% param)
    
    # Use parameter bounds rather than bounds of sampled points
    df[[param]] = normalise_0to1(
      x = df[[param]], 
      x_min = bounds[[param_idx, 1]], 
      x_max = bounds[[param_idx, 2]])
  }
  
  return(df)
}

# ---------------------------------------------------------
# Transform parameters to normal distribution about zero
# ---------------------------------------------------------
standardise_params = function(df, params, response) {
  
  # Initiate list to store details of the transformation
  s_list = list()
  
  # Loop through all columns (including response)
  for (param in c(params, response)) {
    
    # Extract the means and standard deviation
    p_mean = df[[param]] %>% mean()
    p_sd   = df[[param]] %>% sd()
    
    # Use these to standardise the data
    df[[param]] = (df[[param]] - p_mean) / p_sd
    
    # Store the mean and standard deviation
    s_list[[param]] = 
      data.table(p = param, mean = p_mean, sd = p_sd)
  }
  
  # Squish transformation details into single datatable
  s_df = rbindlist(s_list)
  
  return(list(df, s_df))
}

