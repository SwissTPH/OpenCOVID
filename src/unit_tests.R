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
run_model_test = function(o, scenario = "baseline") {
  
  # TODO: Use calibrated contacts if it exists
  
  # Only continue if specified by do_step
  if (!is.element(0, o$do_step)) return()
  
  message("* Testing single model simulation")
  
  # Run model for the defined scenario (see model.R)
  result = model(o, scenario, seed = 1, do_plot = FALSE, verbose = "bar")
  
  # Aggregative and summarise raw model output (see postprocess.R)
  result = process_results(result, result$output)
  
  # Save result as an RDS file
  saveRDS(result, paste0(o$pth$testing, "model_test.rds"))
  
  # Uncomment to load previously simulated file
  # result = readRDS(paste0(o$pth$testing, "model_test.rds"))
  
  message(" - Plotting network properties")

  # Network figure 1) Series of network-related properties
  fig_name = c("Test simulation", scenario, "Network properties")
  plot_network_properties(o, fig_name, result$input, result$network)
  
  # Network figure 2) Age matrix of contact density per age
  fig_name = c("Test simulation", scenario, "Contact matrices")
  plot_contact_matrices(o, fig_name, result$input, result$network)

  message(" - Plotting epidemiological outcomes")

  # Plot metrics against the data
  fig_name = c("Test simulation", scenario)
  plot_temporal(o, fig_name, plot_file = result, alt_baseline = scenario)

  # Plot variants over time the data
  fig_name = c("Test simulation", scenario, "Variants")
  plot_temporal(o, fig_name, plot_file = result, alt_baseline = scenario,
                plot_by = "variant")

  # Same plot but using stacked areas
  fig_name = c("Test simulation", scenario, "Variants", "Area")
  plot_temporal(o, fig_name, plot_file = result, alt_baseline = scenario,
                plot_by = "variant", plot_geom = "area")

  # Plot variants over time the data
  fig_name = c("Test simulation", scenario, "Age")
  plot_temporal(o, fig_name, plot_file = result, alt_baseline = scenario,
                plot_by = "age", plot_geom = "area")

  # Same plot but using stacked areas
  fig_name = c("Test simulation", scenario, "Priority groups")
  plot_temporal(o, fig_name, plot_file = result, alt_baseline = scenario,
                plot_by = "priority_group", plot_geom = "area")
  
  # Plot metrics differentiated by layer of transmission
  fig_name = c("Test simulation", scenario, "Layers")
  plot_temporal(o, fig_name, plot_file = result, alt_baseline = scenario,
                plot_by = "layer")
  
  # Histogram of number of infections per person
  fig_name = c("Test simulation", scenario, "Number of infections")
  plot_num_infections(o, fig_name, plot_file = result, alt_baseline = scenario)
}

