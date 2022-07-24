############################################################
# EMULATOR
#
# Train an algorithm to interpolate the whole solution space
# using model simulations and associated likelihood.
#
# THE STEPS:
# - Identify parameters to be fitted
# - Sample a load of points between bounds
# - Run model with for all of these parameter sets (see simulate.R)
# - Calculate liklihood function for each parameter set (see likelihood.R)
# - Train an ML emulator on the sampled points
#
############################################################

# ---------------------------------------------------------
# Parent function to run model emulator with adaptive sampling
# ---------------------------------------------------------
run_emulator = function(o) {
  
  # Only continue if specified by do_step
  if (!is.element(1, o$do_step)) return()
  
  message("* Training model emulator")
  
  # ---- Load data ----
  
  # Load data directly from source and save in cache
  #
  # NOTE: Workaround as cluster jobs can't seem to directly load data from source
  d = load_data(o, from_cache = FALSE)  # See load_data.R
  
  # Calibration test - create synthetic data from user-defined parameters
  if (!is.null(o$known_params))
    d = create_synthetic_data(o, d)  # See load_data.R
  
  # Set the random number seed generator if desired
  #
  # NOTE: Done after creating synthetic data as we may want to set a seed there too
  if (o$emulator_reproducible == TRUE)
    set.seed(1)
  
  # ---- Display size of this job ----
  
  message(" - Working with ", nrow(o$calibration_df), " parameters in total:")
  
  # Number of parameters global and canton-specific
  n_global = sum(o$calibration_df$scope == "global")
  n_local  = nrow(o$calibration_df) - n_global
  
  message("  > ", n_global, " global parameters")
  message("  > ", n_local,  " canton-specific parameters")
  
  # ---- Adaptive sampling loop ----
  
  # Quick sanity queck on value of resample_iters
  if (!(is.numeric(o$resample_iters) & o$resample_iters >= 0))
    stop("'", o$resample_iters, "' is not a valid choice for resample_iters")
  
  # Initiate objects updated in adaptive sampling loop
  param_sets = emulator = NULL
  
  # Path to store and extract samples from changes on each iteration
  o$pth$samples_stem = o$pth$samples
  
  # Main (re)sampling loop - start at zero to represent initial run
  for (sample_num in 0 : o$resample_iters) {
    
    message(" - Sampling iteration ", sample_num + 1, 
            " (of ", o$resample_iters + 1, ")")
    
    # Update file path to relevant samples 
    o$pth$samples = paste0(o$pth$samples_stem, sample_num, file_sep())
    
    # If it does not already exist, create it
    if (!dir.exists(o$pth$samples)) 
      dir.create(o$pth$samples, recursive = TRUE)
    
    # Sample points from parameter space
    param_sets = get_samples(o, param_sets, sample_num, emulator)
    
    # Run the model with these samples and calculate liklihood
    param_sets = run_samples(o, param_sets, sample_num)
    
    # Train an ML model emulator on the simulated samples
    emulator = train_emulator(o, param_sets, sample_num)
    
    # Plot performance of the emulator
    plot_emulator(o, sample_num = sample_num)
    
    # if (stop_criteria(o, param_sets, sample_num) == TRUE)
    #   break
  }
  
  # ---- Produce emulator diagnostics ----
  
  # Store all samples combined in a database - for now just save as an RDS file
  saveRDS(param_sets, paste0(o$pth$emulator, "all_samples.rds"))
  
  # Store data in final emulator output
  emulator$data = d
  
  # Save the final fitted emulator to file
  saveRDS(emulator, file = paste0(o$pth$emulator, "model_emulator.rds"))
  
  # Plot emulator performance (train vs test)
  plot_emulator(o)
  
  # Plot emulator performance (resampling iterations)
  plot_emulator(o, by_resampling = TRUE)
  
  # Plot parameter-likelihood relationship
  plot_param_likelihood(o)
}

# ---------------------------------------------------------
# Sample points from parameter space - uses acquisution function on later calls
# ---------------------------------------------------------
get_samples = function(o, param_sets, sample_num, emulator) {
  
  # Number of points to sample is different for initial and adaptive sampliong calls
  if (sample_num == 0) n_samples = o$n_init_samples else n_samples = o$acquisition_points
  
  # Form parameter ranges into a matrix (n x 2)
  param_ranges = matrix(c(o$calibration_df$lower_bound, 
                          o$calibration_df$upper_bound), ncol = 2)
  
  # Latin hypercube sample across this space
  param_samples = tgp::lhs(n_samples, param_ranges)
  
  # Standardise values to (mean of 0 and standard deviation of 1)
  if (sample_num == 0) {
    use_samples = base::scale(param_samples)
    
  } else {  # On resamples we already have this scaling - apply it
    resamples = t(param_samples) - param_sets$center
    resamples = t(resamples / param_sets$scale)
    
    # Calculate expected improvement for each and select those with the highest
    best_resamples = as.data.frame(resamples) %>%
      mutate(expected_improvement = crit_EI(resamples, emulator$model)) %>%
      arrange(desc(expected_improvement)) %>%
      filter(expected_improvement > o$min_value_accepted) %>%
      slice(1 : min(n(), o$max_samples_accepted)) %>%
      select(-expected_improvement) %>%
      as.matrix()
    
    # Calculate pairwise distance between all selected points
    dist_resamples = quiet(philentropy::distance(best_resamples, 
                                                 method = o$distance_method))
    
    # Scale distance tolerated by number of dimensions (x' = x * 10^sqrt(n))
    #
    # NOTE: Currently experimenting with this - more reading needed
    dist_tolerance = o$distance_tolerance * 10 ^ sqrt(nrow(o$calibration_df)) 
    
    # Loop through all but the first and remove if too close
    for (i in seq_len(nrow(best_resamples))[-1]) {
      
      # Distance of this point to all better points
      dist_from_better = dist_resamples[i, 1 : (i - 1)]
      
      # If too close to any better, we'll remove this point 
      if (!all(dist_from_better > dist_tolerance, na.rm = TRUE)) {
        
        # Trivialise in sample set matrix
        best_resamples[i, ] = NA
        
        # Also trivialise in distance matrix so we don't compare to a point that's not there
        dist_resamples[i, ] = NA
        dist_resamples[, i] = NA
      }
    }
    
    # Remove parameter sets that have been trivialised for being too close to others
    use_idx     = rowSums(is.na(best_resamples)) != ncol(best_resamples)
    use_samples = as.matrix(best_resamples[use_idx, ])
  }
  
  # Repeat these samples n_seeds_emulator times
  all_samples = matrix(rep(use_samples, each = o$n_seeds_emulator), 
                       ncol = nrow(o$calibration_df))
  
  # Set column names as parameter names
  colnames(all_samples) = o$calibration_df$names
  
  # Parameter ID numbers start from 1 on first iteration
  param_nums = 1 : nrow(use_samples)
  
  # Continue to increment for later iterations
  if (sample_num > 0)
    param_nums = param_nums + max(param_sets$samples$param_id)
  
  # Seed numbers are more straightforward
  seed_nums = 1 : o$n_seeds_emulator
  
  # We'll create simulation IDs with a full factorial merge of parameter set and seed numbers
  p_id = paste0("p", str_pad(param_nums, 6, pad = "0"))
  s_id = paste0("s", str_pad(seed_nums,  4, pad = "0"))
  
  # To complete the ID, use the adaptive sampling index
  a_id = paste0("a", str_pad(sample_num, 3, pad = "0"))
  
  # Merge, bind, collapse, and format into strings
  param_set_ids = cbind(a_id, base::merge(p_id, s_id))
  param_set_ids = sort(unite(param_set_ids, "x", sep = "")$x)
  
  # Initiate a dataframe of IDs and a likelihhod placeholder
  param_set_df = data.frame(sim_id     = param_set_ids, 
                            sample_id  = sample_num, 
                            param_id   = rep(param_nums, each = o$n_seeds_emulator), 
                            seed       = rep(seed_nums, times = length(param_nums)), 
                            likelihood = NA)
  
  # Bind parameter set values
  param_set_df = cbind(param_set_df, all_samples)
  
  # On initial call we want to constrct param_sets list
  if (is.null(param_sets)) {
    
    # Store samples along with scale and center of parameter standardisation
    param_sets = list(samples = param_set_df, 
                      center  = base::attributes(use_samples)[["scaled:center"]], 
                      scale   = base::attributes(use_samples)[["scaled:scale"]])
    
  } else {  # On last calls we just need to append
    param_sets$samples = rbind(param_sets$samples, param_set_df)
  }
  
  # Store samples in a database - for now just save as an RDS file
  saveRDS(param_sets, paste0(o$pth$samples, "samples", sample_num, ".rds"))
  
  return(param_sets)
}

# ---------------------------------------------------------
# Run the model with these samples and calculate liklihood
# ---------------------------------------------------------
run_samples = function(o, param_sets, sample_num) {
  
  # All parameter sets we are interested on on this sampling iteration
  param_set_df = todo = param_sets$samples %>%
    filter(sample_id == sample_num)
  
  # Simulations that already exist (ie from previous analyses)
  sims_exist = str_remove(list.files(o$pth$samples), ".rds")
  sims_exist = sims_exist[sims_exist %in% param_set_df$sim_id]
  
  # Total number of simulations to run
  n_jobs = nrow(param_set_df)  # Including any already simulated
  n_done = length(sims_exist)  # Number already simulated
  
  message("  > Calculating likelihood for ", thou_sep(n_jobs), " simulations")
  
  # For most applications we'll want to use the cluster
  use_cluster = TRUE
  
  # We will overwrite - explain to user
  if (n_done > 0 && o$overwrite_fit == TRUE)
    message("   - Overwriting ", thou_sep(n_done), " previously completed simulations")
  
  # We won't overwrite but will recalculate - do this only for existing simulations
  if (n_done > 0 && o$overwrite_fit == FALSE && o$recalculate_fit == TRUE) {
    message("   - Recalculating for ", thou_sep(n_done), " existing simulations")
    
    # An application we don't need the cluster for
    use_cluster = FALSE
    
    # Update list of simulations we want to recalculate the likelihood for
    todo   = todo %>% filter(sim_id %in% sims_exist)
    n_jobs = nrow(todo)
    
    # Overwrite currently existing jobs so we can recalculate (this could potentially be just a subset)
    todo_sets = list(samples = todo, 
                     center  = param_sets$center, 
                     scale   = param_sets$scale)
    
    # Do the overwriting
    saveRDS(todo_sets, paste0(o$pth$samples, "samples", sample_num, ".rds"))
  }
  
  # We won't overwrite and won't recalculate - only rerun any missing simulations
  if (n_done > 0 && o$overwrite_fit == FALSE && o$recalculate_fit == FALSE) {
    message("   - Skipping ", thou_sep(n_done), " previously completed simulations")
    
    # Update list of simulations to be run
    todo   = todo %>% filter(!sim_id %in% sims_exist)
    n_jobs = nrow(todo)
    
    # We'll overwrite param_sets with this so remaining jobs can be done on the cluster
    todo_sets = list(samples = todo, 
                     center  = param_sets$center, 
                     scale   = param_sets$scale)
    
    # Do the overwriting
    saveRDS(todo_sets, paste0(o$pth$samples, "samples", sample_num, ".rds"))
  }
  
  # ---- Simulate model for sampled parameter sets ----
  
  # Check that we still need to simulate
  if (n_jobs > 0) {
    
    # Define simulation type using (re)sampling iteration number
    calibration_num = paste("calibration", sample_num, sep = "::")
    
    # Submit these jobs to the cluster (see myRfunctions.R)
    if (use_cluster == TRUE) {
      submit_cluster_jobs(o, n_jobs, "bash_submit.sh", calibration_num)
      
    } else {  # If just recalculating likelihood it's quicker to not use the cluster
      
      # Reset samples path for consistency in do_simulate
      o_pth_samples = o$pth$samples
      o$pth$samples = o$pth$samples_stem
      
      # Initiate a progress bar
      pb = start_progress_bar(n_jobs)
      
      # Loop through all simulations 
      for (i in seq_len(n_jobs)) {
        
        # Call do_simulate directly
        do_simulate(calibration_num, i, o = o)
       
        # Update progress bar
        setTxtProgressBar(pb, i)
      }
      
      # Close progress bar
      close(pb)
      
      # Reapply samples path
      o$pth$samples = o_pth_samples
    }
  }
  
  # ---- Load likelihood values from simulations ----
  
  # Loop through to load and format results - ideally use a DB and simply call from there
  for (param_set_id in param_set_df$sim_id) {
    
    # Construct file name and path
    param_set_file = paste0(o$pth$samples, param_set_id, ".rds")
    
    # Load likelihood from file (one from each cluster job)
    if (file.exists(param_set_file)) {
      result = readRDS(param_set_file)
      
      # Store this in the parameter set dataframe
      param_set_df$likelihood[param_set_df$sim_id == param_set_id] = result$likelihood
    }
  }
  
  # Proportion of simulations that have been successful
  p_success = sum(!is.na(param_set_df$likelihood)) / nrow(param_set_df)
  
  message("  > ", round(100 * p_success, 2), "% of simulations successful")
  
  # Throw an error if too many failures - the emulator needs data!
  if (p_success < (1 - o$failure_tolerance))
    stop("Too many simulations have failed - some investigation needed")
  
  # Overwrite likelihood outcomes in param_sets list
  update_idx = param_sets$samples$sample_id == sample_num
  param_sets$samples[update_idx, ] = param_set_df
  
  # Overwrite these samples in the database - for now just saving as an RDS file
  saveRDS(param_sets, paste0(o$pth$samples, "samples", sample_num, ".rds"))
  
  return(param_sets)
}

# ---------------------------------------------------------
# Simulate model and train emulator given parameter sets
# ---------------------------------------------------------
train_emulator = function(o, param_sets, sample_num) {
  
  # File path for emulator we are about to train
  emulator_pth = paste0(o$pth$emulator, "emulator", sample_num, ".rds")
  
  # If this already exists and we don't want to overwrite, load and return
  if (file.exists(emulator_pth) && o$overwrite_fit == FALSE) {
    emulator = readRDS(emulator_pth)
    
    return(emulator)
  }
  
  # ---- Split into test and train sets ----
  
  # Use only the best max_train_samples parameter sets
  best_samples = param_sets$samples %>%
    group_by(param_id) %>%
    summarise(best = mean(likelihood, na.rm = TRUE)) %>%
    slice_min(best, n = o$max_train_samples)
  
  # Remove failures from data set we'll be working with to train emulator
  success_df = param_sets$samples %>%
    filter(param_id %in% best_samples$param_id,
           !is.na(likelihood))
  
  # Parameter IDs that have had at least one success
  success_ids = unique(success_df$param_id)
  
  # Number of parameter sets we'll use for testing
  n_test = round(length(success_ids) * o$train_test_split)
  
  # Randomly sample test parameter set numbers
  test_sets = sort(sample(success_ids, n_test))
  
  # Shorthand parameter names
  param_names = o$calibration_df$names
  n_params    = length(param_names)
  
  # Seperate out training set
  train_data = success_df %>%
    filter(!(param_id %in% test_sets)) %>%
    select(one_of(param_names), likelihood)
  
  # Dataframe of all parameter sets (excluding seeds) with their IDs
  param_id_df = success_df %>% 
    select(param_id, one_of(param_names)) %>% 
    unique()
  
  # Seperate out testing set and take the mean across seeds
  test_data = success_df %>%
    filter(param_id %in% test_sets) %>%
    group_by(param_id) %>%
    summarise(likelihood = mean(likelihood)) %>%
    left_join(param_id_df, by = "param_id") %>%
    select(one_of(param_names), likelihood) %>%
    as.data.frame()
  
  # Sanity check that these dataframes have the same columns
  if (!identical(names(train_data), names(test_data)))
    stop("Inconsistent training and test sets")
  
  # ---- Train the model ----
  
  message("  > Fitting model emulator")
  
  # Find and remove duplicates (we've already done a transformation so turn these options off)
  prep_data = hetGP::find_reps(X = as.matrix(select(train_data, -likelihood)), 
                               Z = as.matrix(select(train_data, likelihood)), 
                               rescale = FALSE, normalize = FALSE)
  
  # Construct input for mleHetGP model
  X = list(X0 = as.matrix(prep_data$X0), 
           Z0 = as.matrix(prep_data$Z0),
           mult = prep_data$mult)
  
  # Run GP model (mleHetGP from hetGP package)
  model_emulator = hetGP::mleHetGP(X = X, Z = prep_data$Z, 
                                   covtype  = o$gp_kernel,
                                   lower    = rep(0.0001, n_params), 
                                   upper    = rep(10, n_params), 
                                   maxit    = o$gp_max_iter)
  
  # ---- Prepare output ----
  
  # Append model and data used to calculate the likelihood
  emulator = list(model      = model_emulator, 
                  data       = NA,  # Placeholder
                  train_sets = train_data, 
                  test_sets  = test_data, 
                  center     = param_sets$center, 
                  scale      = param_sets$scale)
  
  # Assess the model with both testing and training data
  emulator = test_emulator(emulator, train_data, test_data)
  
  # Save the fitted emulator for this sampling number to file
  saveRDS(emulator, file = emulator_pth)
  
  return(emulator)
}

# ---------------------------------------------------------
# Run GP with test data to determine model error
# ---------------------------------------------------------
test_emulator = function(emulator, train_data, test_data) {
  
  message("  > Assessing emulator accuracy")
  
  # Define indices of parameters and response
  n_params     = ncol(train_data) - 1
  param_idx    = 1 : n_params
  response_idx = n_params + 1
  
  # Apply the model on the test data
  predict_test = stats::predict(x = as.matrix(test_data[, param_idx]), 
                                object = emulator$model)$mean
  
  # The actual response associated with the test data
  actual_test = test_data[, response_idx]
  
  # Repeat this for the training data
  predict_train = stats::predict(x = as.matrix(train_data[, param_idx]), 
                                 object = emulator$model)$mean
  
  # The actual response associated with the test data
  actual_train = train_data[, response_idx]
  
  # Concatenate output into list
  emulator = c(emulator, list(actual_test   = actual_test, 
                              predict_test  = predict_test, 
                              actual_train  = actual_train, 
                              predict_train = predict_train))
  
  # hetGP::scores(emulator$model, Xtest, Ztest, return.rmse = FALSE)
  
  return(emulator)
}

# ---------------------------------------------------------
# Check whether we can prematurely stop adaptive sampling
# ---------------------------------------------------------
stop_criteria = function(o, param_sets, sample_num) {
  
  browser()
  
  stop_flag = FALSE
  
  return(stop_flag)
}

