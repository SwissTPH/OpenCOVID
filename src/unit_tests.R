###########################################################
# UNIT TESTS
#
# A series of self-contained functions to test out different
# aspects of model functionality.
#
###########################################################

# ---------------------------------------------------------
# Run a single simulation and produce a few standard outputs
# ---------------------------------------------------------
run_model_test = function(o, scenario = "baseline", rerun = TRUE) {
  
  # Only continue if specified by do_step
  if (!is.element(0, o$do_step)) return()
  
  message("* Testing single model simulation")
  
  # ---- Simulate model ----
  
  # File to save after a simulation
  test_file = paste0(o$pth$testing, "model_test.rds")
  
  # Check whether we want to simulate (or simply plot a previously created file)
  if (rerun == TRUE || !file.exists(test_file)) {
    
    # Load calibration result (if it exists)
    fit_list = load_calibration(o, throw_error = FALSE)$best
    
    # If uncertainty defined, take the average over the distribution(s)
    uncert_list = sample_average(o)  # See uncertainty.R
    
    # Run model for the defined scenario (see model.R)
    result = model(o, scenario,
                   seed    = 1,
                   fit     = fit_list,
                   uncert  = uncert_list,
                   do_plot = FALSE,
                   verbose = "bar")
    
    # Aggregative and summarise raw model output (see postprocess.R)
    result = process_results(result, result$output)
    
    # Save result as an RDS file
    saveRDS(result, test_file)
  }
  
  # Load simulated file - we may have skipped simulating
  result = readRDS(test_file)
  
  # ---- Plot outcomes ----
  
  message(" - Plotting network properties")

  # Network figure 1) Series of network-related properties
  fig_name = c("Test simulation", scenario, "Network properties")
  plot_network_properties(o, fig_name, result$input, result$network)

  # Network figure 2) Age matrix of contact density per age
  fig_name = c("Test simulation", scenario, "Contact matrices")
  plot_contact_matrices(o, fig_name, result$input, result$network)

  message(" - Plotting epidemiological outcomes")

  # Plot all available metrics over time
  fig_name = c("Test simulation", scenario)
  plot_temporal(o, fig_name, plot_file = result, alt_baseline = scenario)

  # Plot all metrics in cumulative form
  fig_name = c("Test simulation", scenario, "Cumulative")
  plot_temporal(o, fig_name, plot_file = result, alt_baseline = scenario,
                cumulative = TRUE)

  # Plot metrics by variant
  fig_name = c("Test simulation", scenario, "Variants")
  plot_temporal(o, fig_name, plot_file = result, alt_baseline = scenario,
                plot_by = "variant")

  # Same plot but using stacked areas
  fig_name = c("Test simulation", scenario, "Variants", "Area")
  plot_temporal(o, fig_name, plot_file = result, alt_baseline = scenario,
                plot_by = "variant", plot_geom = "area")

  # Plot metrics by age group
  fig_name = c("Test simulation", scenario, "Age")
  plot_temporal(o, fig_name, plot_file = result, alt_baseline = scenario,
                plot_by = "age", plot_geom = "area")

  # Plot metrics by priority group
  fig_name = c("Test simulation", scenario, "Priority groups")
  plot_temporal(o, fig_name, plot_file = result, alt_baseline = scenario,
                plot_by = "priority_group", plot_geom = "area")

  # Plot metrics by latest vaccine received
  fig_name = c("Test simulation", scenario, "Vaccine type")
  plot_temporal(o, fig_name, plot_file = result, alt_baseline = scenario,
                plot_by = "vaccine_type", plot_geom = "area")

  # Plot metrics by number of vaccine doses received
  fig_name = c("Test simulation", scenario, "Vaccine doses")
  plot_temporal(o, fig_name, plot_file = result, alt_baseline = scenario,
                plot_by = "vaccine_doses", plot_geom = "area")

  # Histogram of number of infections per person
  fig_name = c("Test simulation", scenario, "Number of infections")
  plot_num_infections(o, fig_name, plot_file = result, alt_baseline = scenario)
}

