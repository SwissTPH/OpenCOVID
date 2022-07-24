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
  if (!is.element(3, o$do_step)) return()
  
  message("* Producing outputs and figures")
  
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
    fig_name = c("Baseline", "Variants", "Area")
    plot_temporal(o, fig_name, plot_by = "variant", plot_geom = "area")
    
    # Plot age groups
    fig_name = c("Baseline", "Age groups", "Area")
    plot_temporal(o, fig_name, plot_by = "age", plot_geom = "area")
    
    # Plot vaccine and treatment priority groups
    fig_name = c("Baseline", "Priority groups", "Area")
    plot_temporal(o, fig_name, plot_by = "priority_group", plot_geom = "area")
    
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
  alt_scenarios = names(parse_yaml(o, "*read*", read_array = TRUE)[-1])
  
  # Which of these are children of array scenarios
  array_parents  = get_array_parents(alt_scenarios)
  array_children = get_array_children(alt_scenarios, array_parents)
  
  # We only want to consider single, non-array scenarios here
  single_scenarios = setdiff(alt_scenarios, array_children)
  
  # Check plotting flag and whether any single alternative scenarios have been defined
  if (o$plot_scenarios == TRUE && length(single_scenarios) > 0) {
    
    message(" - Alternative scenarios")
    
    # Alternative scenarios
    plot_temporal(o, "Alternative scenarios", scenarios = single_scenarios)
    
    # Alternative scenarios - all metrics
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
  
  # ---- Multi-dimensional array heat maps ----
  
  # Check plotting flag
  if (o$plot_heatmaps == TRUE) {
    
    # Identify multi-dimensional array parents
    multidim_arrays = get_array_parents(alt_scenarios, multidim = TRUE)
    
    # Loop through all distinct parents
    for (this_array in multidim_arrays) {
      
      message(" - Multidimensional array heat map: ", this_array)
      
      # Produce heat map of moving parts: Cumulative new infections
      fig_name = c("Array scenarios", "Heat map", this_array, "New infections")
      plot_heatmap(o, fig_name, this_array, plot_metrics = "all_new_infections", summarise = "sum")
      
      # Produce heat map of moving parts: Cumulative deaths
      fig_name = c("Array scenarios", "Heat map", this_array, "Deaths")
      plot_heatmap(o, fig_name, this_array, plot_metrics = "deaths", summarise = "sum")
      
      # Produce heat map of moving parts: Mean R_eff over time
      fig_name = c("Array scenarios", "Heat map", this_array, "R_eff")
      plot_heatmap(o, fig_name, this_array, plot_metrics = "R_effective", summarise = "mean")
      
      # Produce heat map of moving parts: Cumulative hospital admissions
      fig_name = c("Array scenarios", "Heat map", this_array, "Hospital admissions")
      plot_heatmap(o, fig_name, this_array, plot_metrics = "hospital_admissions", summarise = "sum")
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
    fig_name = c("Duration distributions", scenario_assumptions)
    plot_immunity_profiles(o, fig_name, model_input)
    
    # Plot calibration emulator and optimisation performance (outcomes of step 1)
    plot_emulator(o, "Emulator performance")
    plot_optimisation(o, "Optimisation performance")
  }
  
  # ---- Custom figures ----
  
  # Plot booster manuscript figures (June 2022)
  if (o$plot_manuscript == TRUE)
    manuscript_figures(o)
  
  # Also plot custom results defined within my_results.R (if it exists)
  if (o$plot_custom == TRUE && file.exists("my_results.R"))
    my_results(o)
}

