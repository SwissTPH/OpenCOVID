###########################################################
# RESULTS
#
# Plot standard set of epi results.
#
###########################################################

# ---------------------------------------------------------
# Parent function for plotting standard results set
# ---------------------------------------------------------
run_results = function(o) {
  
  # Only continue if specified by do_step (or forced)
  if (!is.element(4, o$do_step)) return()
  
  message("* Producing results and figures")
  
  # Load baseline file to extract model inputs - useful for certain type of plots
  base_input = try_load(o$pth$scenarios, "baseline")$input
  
  # ---- Baseline ----
  
  # Check plotting flag
  if (o$plot_baseline == TRUE) {
    
    message(" - Baseline")
    
    # First check that baseline file exists
    if (!file.exists(paste0(o$pth$scenarios, "baseline.rds")))
      stop("Cannot find baseline results - check you have run the 'scenarios' step")
    
    # Plot all metrics (per 100,000 people per day by default)
    plot_temporal(o, "Baseline")

    # Plot virus variants
    fig_name = c("Baseline", "Variants")
    plot_temporal(o, fig_name, plot_by = "variant")

    # Plot virus variants - area plot
    fig_name = c("Baseline", "Variants")
    plot_temporal(o, fig_name, plot_by = "variant", plot_geom = "area")

    # Plot age groups
    fig_name = c("Baseline", "Age groups")
    plot_temporal(o, fig_name, plot_by = "age", plot_geom = "area")

    # Plot vaccine and treatment priority groups
    fig_name = c("Baseline", "Priority groups")
    plot_temporal(o, fig_name, plot_by = "priority_group", plot_geom = "area")

    # Plot types of vaccine distributed
    fig_name = c("Baseline", "Vaccine types")
    plot_temporal(o, fig_name, plot_by = "vaccine_type", plot_geom = "area")

    # Plot types of vaccine distributed
    fig_name = c("Baseline", "Vaccine doses")
    plot_temporal(o, fig_name, plot_by = "vaccine_doses", plot_geom = "area")

    # Histogram of number of infections per person
    fig_name = c("Baseline", "Number of infections")
    plot_num_infections(o, fig_name)
    
    # Check whether we want cumulative plots
    if (o$plot_cumulative == TRUE) {
      
      # Plot all metrics in cumulative form
      fig_name = c("Baseline", "Cumulative")
      plot_temporal(o, fig_name, cumulative = TRUE)
    }
  }
  
  # ---- Alternative scenarios ----
  
  # Alternative scenarios associated with this analysis
  all_scenarios = names(parse_yaml(o, "*read*", read_array = TRUE))
  alt_scenarios = all_scenarios[-1]  # Remove reference to baseline
  
  # Which of these are children of array scenarios
  array_parents = get_array_parents(alt_scenarios)
  all_children  = get_array_children(alt_scenarios, array_parents)
  
  # We only want to consider single, non-array scenarios here
  single_scenarios = setdiff(alt_scenarios, all_children)
  
  # Check plotting flag and whether any single alternative scenarios have been defined
  if (o$plot_scenarios == TRUE && length(single_scenarios) > 0) {
    
    message(" - Alternative scenarios")
    
    # Alternative scenarios - all metrics
    plot_temporal(o, "Alternative scenarios", scenarios = single_scenarios) 
    
    # Alternative scenarios - cumulative bars
    fig_name = c("Alternative scenarios", "Impact bars")
    plot_impact(o, fig_name, scenarios = single_scenarios)
    
    # Histogram of number of infections per person
    fig_name = c("Alternative scenarios", "Number of infections")
    plot_num_infections(o, fig_name, scenarios = single_scenarios)
    
    # Alternative scenarios - cumulative
    if (o$plot_cumulative == TRUE) {
      fig_name = c("Alternative scenarios", "Cumulative")
      plot_temporal(o, fig_name, scenarios = single_scenarios, cumulative = TRUE)
    }
  }
  
  # ---- Array scenario bundles ----
  
  # Check plotting flag
  if (o$plot_arrays == TRUE) {
    
    # Loop through all distinct parent array scenarios
    for (array_parent in array_parents) {
      
      # All children scenarios associated with this parent
      array_children = get_array_children(alt_scenarios, array_parent)
      
      # Only produce plot if reasonable number of children
      if (length(array_children) <= o$max_scenarios) {
        
        message(" - Array scenarios: ", array_parent)
        
        # Array scenarios for this parent
        fig_name = c("Array scenarios", array_parent)
        plot_temporal(o, fig_name, scenarios = array_children)
        
        # Array scenarios for this parent - all metrics
        fig_name = c("Array scenarios", "Impact bars", array_parent)
        plot_impact(o, fig_name, scenarios = array_children)
        
        # Array scenarios for this parent
        if (o$plot_cumulative == TRUE) {
          fig_name = c("Array scenarios", "Cumulative", array_parent)
          plot_temporal(o, fig_name, scenarios = array_children, cumulative = TRUE)
        }
      }
    }
  }
  
  # ---- Array LHC endpoints ----
  
  # Check plotting flag
  if (o$plot_endpoints == TRUE) {
    
    # All array parents that are parents of LHC scenarios
    lhc_parents = array_parents[is_parent_lhc(o, array_parents)]
    
    # All modelled endpoints
    endpoints = names(base_input$scenario_lhc_endpoints)
    
    # Loop through all distinct parents
    for (lhc_parent in lhc_parents) {
      
      message(" - Array LHC endpoints: ", lhc_parent)
      
      # Plot diagnostics for all endpoints
      plot_lhc_diagnostics(o, "LHC diagnostics", lhc_parent)
      
      # Array info for this parent
      array_info = try_load(o$pth$array_info, lhc_parent)
      
      # Construct most basic variable list, varying only the first
      vars_list = c("free", rep("mid", nrow(array_info$vars) - 1))
      vars_list = as.list(setNames(vars_list, array_info$vars$variable_id))
      
      # Loop through all endpoints
      for (endpoint in endpoints) {
        
        # Plot the most basic variable configuration for each endpoint
        #
        # NOTE: The user is encouraged to experiment with different variable configurations
        fig_name = c("LHC predictions", lhc_parent, endpoint)
        plot_lhc_endpoints(o, fig_name = fig_name,
                           endpoint    = endpoint,
                           parent      = lhc_parent,
                           vars        = vars_list)
      }
    }
  }
  
  # ---- Multi-dimensional array heat maps ----
  
  # TODO: Also plot heatmaps for LHC scenarios - have to have a think about the conditions
  
  # Check plotting flag
  if (o$plot_heatmaps == TRUE) {
    
    # Identify multi-dimensional array parents
    multidim_arrays = get_array_parents(alt_scenarios, multidim = TRUE)
    
    # Loop through all distinct parents
    for (this_array in multidim_arrays) {
      
      message(" - Multidimensional array heat map: ", this_array)
      
      # Produce heat map of moving parts: Cumulative new infections
      plot_heatmap(o, fig_name  = c("Heat map", this_array, "New infections"),
                   array        = this_array,
                   plot_metrics = "all_new_infections",
                   summarise    = "sum")

      # Produce heat map of moving parts: Cumulative hospital admissions
      plot_heatmap(o, fig_name  = c("Heat map", this_array, "Hospital admissions"),
                   array        = this_array,
                   plot_metrics = "hospital_admissions",
                   summarise    = "sum")

      # Produce heat map of moving parts: Mean Re over time
      plot_heatmap(o, fig_name  = c("Heat map", this_array, "Re"),
                   array        = this_array,
                   plot_metrics = "Re",
                   summarise    = "mean")
    }
  }
  
  # ---- Model structure and assumptions ----
  
  # Plot these assumptions for baseline scenario (although it could be any scenrio)
  scenario_assumptions = "baseline"
  
  # Check plotting flag
  if (o$plot_assumptions == TRUE) {
    
    message(" - Model structure and assumptions")
    
    # Load scenario_assumptions file
    results = try_load(o$pth$scenarios, scenario_assumptions)
    
    # Extract model input and network details
    model_input   = results$input
    model_network = results$network
    
    # Network figure 1) Series of network-related properties
    fig_name = c("Network properties", scenario_assumptions)
    plot_network_properties(o, fig_name, model_input, model_network)

    # Network figure 2) Age matrix of contact density per age (single age bins)
    fig_name = c("Contact matrices", scenario_assumptions)
    plot_contact_matrices(o, fig_name, model_input, model_network)

    # Plot disease state durations
    fig_name = c("Duration distributions", scenario_assumptions)
    plot_durations(o, fig_name, model_input)

    # Plot viral load profile
    fig_name = c("Viral load profile", scenario_assumptions)
    plot_viral_load(o, fig_name, model_input)

    # Plot immunity profiles for natural infections and vaccination
    fig_name = c("Immunity profiles", scenario_assumptions)
    plot_immunity_profiles(o, fig_name, model_input)

    # Plot basic severe disease probabilities by age
    fig_name = c("Severity by age", scenario_assumptions)
    plot_severity_age(o, fig_name, model_input)

    # Plot exposure & dose response effects on disease severity
    fig_name = c("Severity factors", scenario_assumptions)
    plot_severity_factors(o, fig_name, model_input)

    # Plot seasonality profile
    #
    # NOTE: This function is able to handle multiple scenarios
    fig_name = c("Seasonality profile", scenario_assumptions)
    plot_seasonality_profile(o, fig_name, model_input)
    
    # Plot air pollution - susceptibility relationship options
    fig_name = c("Air pollution relationships", scenario_assumptions)
    plot_pollution_relationships(o, fig_name, model_input) 

    # Plot parameter uncertainty distributions (for all scenarios at once)
    plot_uncertainty(o, "Parameter uncertainty distributions")
    
    # Plot simulation time taken (all single scenarios)
    #
    # NOTE: Plots get too big when considering all array children
    plot_simulation_time(o, "Simulation time", scenarios = single_scenarios)
  }
  
  # ---- Calibration performance and diagnostics ----
  
  # Check plotting flag
  if (o$plot_calibration == TRUE) {
    
    message(" - Calibration performance")
    
    # Loop through adaptive sampling rounds
    for (round_idx in paste0("r", 0 : base_input$adaptive_sampling$rounds)) {
      
      # Plot each round, as long as results exist (we may have broken out early)
      round_result = paste0(o$pth$fitting, round_idx, "_result.rds")
      if (file.exists(round_result)) {
        
        # Plot calibration emulator and optimisation performance (outcomes of step 1)
        plot_best_samples(o, "Best simulated samples",   round_idx)
        plot_emulator(o,     "Emulator performance",     round_idx)
        plot_optimisation(o, "Optimisation performance", round_idx)
      }
    }
    
    # Plot calibration weight assumptions for metrics, time, and peaks
    plot_calibration_weights(o, "Calibration weights")
    
    # Plot baseline model output with data laid on top (cumulative for epi_data)
    plot_temporal(o, "Calibration plot", 
                  fit_target = TRUE, 
                  cumulative = base_input$calibration_type == "epi_data")
  }
  
  # ---- Custom figures ----
  
  # Plot air pollution manuscript figures (CROI conference 2023)
  if (o$plot_manuscript == TRUE)
    manuscript_figures(o)
  
  # Also plot custom results defined within my_results.R (if it exists)
  if (o$plot_custom == TRUE && file.exists("my_results.R"))
    my_results(o)
}

