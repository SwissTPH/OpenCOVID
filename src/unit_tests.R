###########################################################
# UNIT TESTS
#
# A series of self-contained functions to test out different
# aspects of model functionality.
#
###########################################################

# ---------------------------------------------------------
# Test effect of parameters manually
#
# Example usage: 
# > run_model_test(o)  # To test run with default parameters
# > run_model_test(o, test_params = list(seasonality_scaler = 0, 
#                                        import_initial.CH  = 30))
# ---------------------------------------------------------
run_model_test = function(o, test_params = list()) {
  
  # Only continue if specified by do_step
  if (!is.element(NA, o$do_step)) return()
  
  message("* Testing single model simulation")
  
  if (!is.list(test_params))
    stop("Input should be a list of parameter names with associated values")
  
  # Just look at one canton here
  if (length(o$cantons) > 1)
    stop("Ideal use case for this test is one single canton")
  
  # Use prior means for all parameters as default
  calibration_params = setNames(o$calibration_df$values, o$calibration_df$names)
  
  # Overwrite with user-defined test_params
  calibration_params[names(test_params)] = unlist(test_params)
  
  # Load data (see load_data.R)
  d = load_data(o, from_cache = TRUE)
  d = d[[o$cantons]]
  
  # Generate all model parameters (see parameters.R)
  p = get_parameters(o, d, calibration_params, do_plot = FALSE)
  
  # Run model (see model.R)
  m = model(o, p, seed = 1, short_run = TRUE, verbose = "date")
  
  # Collate results - use trivial bounds as we have only a single simulation
  model_output = m %>%
    select(date, metric, grouping, group, 
           mean  = value, median = value,
           lower = value, upper  = value) %>%
    mutate(scenario = "model_test")
  
  # Append to an analysis file along with data
  a = list(fit_data = d, model_input = p, model_output = model_output)
  
  # Save model outcomes as an RDS file
  saveRDS(a, paste0(o$pth$testing, "model_test.rds"))
  
  # Uncomment to load previously simulated file
  # a = readRDS(paste0(o$pth$testing, "model_test.rds"))
  
  # Likelihood value associated with this parameter set
  # likelihood = likelihood_model(o, m = a$model_output)
  
  message(" - Plotting outcomes")
  
  # Plot metrics against the data
  fig_name = c("Test simulation")
  plot_temporal(o, o$cantons, fig_name, plot_file = a, 
                everything = TRUE, likelihood = NA)
  
  # Plot variants over time the data
  fig_name = c("Test simulation - Variants")
  plot_temporal(o, o$cantons, fig_name, plot_file = a, 
                everything = TRUE, plot_by = "variant")
}

# ---------------------------------------------------------
# Examine difference between simulations with high likelihood
# ---------------------------------------------------------
test_high_likelihood = function(likelihood_tol = 1e5, max_sims = 50) {
  
  message("* Testing high likelihood fits")
  
  # Load options (see options.R)
  o = set_options(quiet = TRUE)
  
  # Just look at one canton here
  if (length(o$cantons) > 1)
    stop("Ideal use case for this test is one single canton")
  
  # All parameter sets we've actually simulated
  param_sets = readRDS(paste0(o$pth$emulator, "all_samples.rds"))
  
  # Select the sims with 'high' likelihood value
  high_samples = param_sets$samples %>%
    select(-sample_id, -param_id) %>%
    arrange(likelihood) %>%
    filter(likelihood <= likelihood_tol) %>%
    slice(1 : max_sims)
  
  # Number of simulations identified
  n_high_sims = nrow(high_samples)
  
  # Sanity check that we've got something to analyse
  if (n_high_sims == 0)
    stop("No 'high' likelihood simulations identified - try relaxing tolerance")
  
  # Create a new set of simulation IDs for out testing work
  high_samples$sim_id = paste0("l", str_pad(1 : n_high_sims, 4, pad = "0"))
  
  # Construct a dataframe similar to param_sets so we can run on the cluster
  high_sets = list(samples = high_samples, 
                   center  = param_sets$center, 
                   scale   = param_sets$scale)
  
  # Store this dataframe on file to be loaded when running cluster jobs
  saveRDS(high_sets, paste0(o$pth$testing, "likelihood_samples.rds"))
  
  # ---- Run all simulations on the cluster ----
  
  # Load data directly from source and save in cache
  #
  # NOTE: Workaround as cluster jobs can't seem to directly load data from source
  load_data(o, from_cache = FALSE)  # See load_data.R
  
  message(" - Running model")
  
  # Submit these jobs to the cluster (see myRfunctions.R)
  submit_cluster_jobs(o, n_high_sims, "bash_submit.sh", "test_likelihood")
  
  # ---- Load and format outcomes ----
  
  browser() # Format of model outcomes is now long format...
  
  # Initiate a master dataframe
  all_sims_df = NULL
  
  # Loop through simulations to load and format model outcomes
  for (sim_id in high_samples$sim_id) {
    
    # Load model outcomes from this simulation
    sim_df = readRDS(paste0(o$pth$testing, sim_id, ".rds"))
    
    # Append simulation ID
    sim_df$sim_id = sim_id
    
    # Continuously append to a master dataframe
    all_sims_df = rbind(all_sims_df, sim_df)
  }
  
  # ---- Plot outcomes ----
  
  message(" - Plotting outcomes")
  
  # Melt down metrics
  plot_df = all_sims_df %>%
    pivot_longer(cols = -c(date, sim_id), 
                 names_to = "metric") %>%
    as.data.frame()
  
  # Plot all simulations as lines
  g = ggplot(plot_df, aes(x = date, y = value, colour = sim_id)) + 
    geom_line(size = 1.5) + 
    facet_wrap(~metric, scales = "free_y", 
               labeller = as_labeller(o$metric_dict))
  
  # Set x axes to user-defined limits (see options.R)
  g = g + scale_x_date(date_breaks = o$x_tick_dates, date_labels = "%d %b") + 
    scale_y_continuous(labels = comma, expand = c(0, NA), limits = c(0, NA))
  
  # Save figure (see plotting.R)
  save_fig(o, g, "Test likelihood")
}

# ---------------------------------------------------------
# Examine what extent stochasticity has on uncertainty
# ---------------------------------------------------------
test_stochasticity = function(n_seeds = 10) {
  
  message("* Running stochasticity test")
  
  # Metrics to plot
  plot_metrics = qc(confirmed, deaths, hospital_beds, icu_beds, R_effective)
  
  # Load options (see options.R)
  o = set_options(quiet = TRUE)
  
  # Just look at one canton here
  if (length(o$cantons) > 1)
    stop("Ideal use case for this test is one single canton")
  
  # Generate a load of sim IDs that vary only by seed ID
  sim_ids = paste0("s", str_pad(1 : n_seeds, 4, pad = "0"))
  
  # Initiate a dataframe of sim IDs and seed numbers
  param_sets = data.frame(sim_id = sim_ids, seed = 1 : n_seeds)
  
  # Store this dataframe on file to be loaded when running cluster jobs
  saveRDS(param_sets, paste0(o$pth$testing, "stochasticity_samples.rds"))
  
  # ---- Run all simulations on the cluster ----
  
  # Load data directly from source and save in cache
  #
  # NOTE: Workaround as cluster jobs can't seem to directly load data from source
  load_data(o, from_cache = FALSE)  # See load_data.R
  
  message(" - Running model")
  
  # Submit these jobs to the cluster (see myRfunctions.R)
  submit_cluster_jobs(o, n_seeds, "bash_submit.sh", "test_stochasticity")
  
  # ---- Load and format outcomes ----
  
  browser() # Format of model outcomes is now long format...
  
  # Initiate a master dataframe
  all_sims_df = NULL
  
  # Loop through simulations to load and format model outcomes
  for (sim_id in sim_ids) {
    
    # Load model outcomes from this simulation
    m = readRDS(paste0(o$pth$testing, sim_id, ".rds"))
    
    # Extract metrics and append seed number
    sim_df = m %>% 
      select(date, one_of(plot_metrics)) %>% 
      mutate(seed = str_extract(sim_id, "s\\d+"))
    
    # Continuously append to a master dataframe
    all_sims_df = rbind(all_sims_df, sim_df)
  }
  
  # ---- Plot outcomes ----
  
  message(" - Plotting outcomes")
  
  # Melt down metrics
  plot_df = all_sims_df %>%
    pivot_longer(cols = plot_metrics, 
                 names_to = "metric") %>%
    as.data.frame()
  
  # Plot all simulations as lines
  g = ggplot(plot_df, aes(x = date, y = value, colour = seed)) + 
    geom_line(size = 1.5) + 
    facet_wrap(~metric, scales = "free_y", 
               labeller = as_labeller(o$metric_dict[plot_metrics]))
  
  # Set x axes to user-defined limits (see options.R)
  g = g + expand_limits(y = c(0, 10)) +  # Restrain upper y limit from being really low
    scale_y_continuous(labels = comma, expand = c(0, 0), limits = c(0, NA)) +
    scale_x_date(date_breaks = o$x_tick_dates, date_labels = "%d %b")
  
  # Save figure (see plotting.R)
  save_fig(o, g, "Test stochasticity")
}

