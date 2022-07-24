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
  
  message("* Producing outputs and figures")
  
  # Set an alternative baseline if desired (set to NULL to plot default baseline)
  alt_baseline = NULL  # eg "tf0_baseline.novax_npikeep"
  
  # Define a few commonly used dates to define when we plot to/from
  data_end = max(o$dates_data)
  
  # ---- Canton baseline figures ----
  
  if (o$plot_baseline == TRUE) {
    
    # Repeat for each canton
    for (canton in o$cantons) {
      
      message(" - Baseline: ", canton)
      
      if (!file.exists(paste0(o$pth$scenario, canton, "_baseline.rds")))
        stop("Cannot find baseline analysis file - check you have run the 'analysis' step")
      
      # Subset of all metrics for a nice full calibration picture
      age_metrics = c("confirmed", "hospital_beds", "icu_beds", "deaths", 
                      "currently_infected", "total_vaccinated")
      
      # Subset of all metrics for a nice full calibration picture
      custom_metrics = c("confirmed", "hospital_beds", "icu_beds", "deaths", 
                         "new_local_infections", "new_importations", "currently_infected", "seasonality", 
                         "seroprevalence", "pop_susceptibility", "npi_effect", "total_vaccinated")
      
      # Set a colour map of interest
      custom_colours = colour_scheme("pals::cols25", n = length(custom_metrics))
      
      # Plot default metrics for each canton (clipped for pretty calibration plot)
      fig_name = c("Calibration", canton)
      plot_temporal(o, canton, fig_name, plot_to = data_end, alt_baseline = alt_baseline)

      # Plot everything in model output for each canton (also clipped)
      fig_name = c("Calibration", canton, "Custom")
      plot_temporal(o, canton, fig_name, plot_to = data_end, plot_metrics = custom_metrics,
                    override_colours = custom_colours, alt_baseline = alt_baseline)

      # Plot everything in model output for each canton (also clipped)
      fig_name = c("Calibration", canton, "Everything")
      plot_temporal(o, canton, fig_name, plot_to = data_end,
                    everything = TRUE, alt_baseline = alt_baseline)
      
      # Plot age group for baseline calibration - area plot
      fig_name = c("Calibration", canton, "Age groups", "Area")
      plot_temporal(o, canton, fig_name, plot_to = data_end, plot_metrics = age_metrics,
                    plot_by = "age", plot_area = TRUE, alt_baseline = alt_baseline)
      
      # Plot virus variants for baseline calibration
      fig_name = c("Calibration", canton, "Variants")
      plot_temporal(o, canton, fig_name, plot_to = data_end, everything = TRUE, 
                    plot_by = "variant", alt_baseline = alt_baseline)
      
      # Plot default metrics for each canton (with future projection)
      fig_name = c("Baseline projection", canton)
      plot_temporal(o, canton, fig_name, alt_baseline = alt_baseline)
      
      # Plot everything in model output for each canton (also clipped)
      fig_name = c("Baseline projection", canton, "Custom")
      plot_temporal(o, canton, fig_name, plot_metrics = custom_metrics, 
                    alt_baseline = alt_baseline)
      
      # Plot all metrics for each canton (with future projection)
      fig_name = c("Baseline projection", canton, "Everything")
      plot_temporal(o, canton, fig_name, everything = TRUE, alt_baseline = alt_baseline)
      
      # Plot virus variants for baseline calibration (with future projection)
      fig_name = c("Baseline projection", canton, "Variants")
      plot_temporal(o, canton, fig_name, everything = TRUE, plot_by = "variant",
                    alt_baseline = alt_baseline)
      
      # Plot virus variants for baseline calibration (with future projection) - area plot
      fig_name = c("Baseline projection", canton, "Variants", "Area")
      plot_temporal(o, canton, fig_name, everything = TRUE, plot_by = "variant",
                    plot_area = TRUE, alt_baseline = alt_baseline)
      
      # Plot age group for baseline calibration (with future projection) - area plot
      fig_name = c("Baseline projection", canton, "Age groups", "Area")
      plot_temporal(o, canton, fig_name, everything = TRUE, plot_by = "age",
                    plot_area = TRUE, alt_baseline = alt_baseline)
      
      # Plot metrics for vaccine priority group for baseline future projection - area plot
      fig_name = c("Baseline projection", canton, "Priority groups", "Area")
      plot_temporal(o, canton, fig_name, everything = TRUE, plot_by = "vaccine_priority",
                    plot_area = TRUE, alt_baseline = alt_baseline)
      
      if (o$plot_cumulative == TRUE) {
        
        # Plot cumulative metrics for each canton
        fig_name = c("Calibration", canton, "Cumulative")
        plot_temporal(o, canton, fig_name, everything = TRUE, cumulative = TRUE, 
                      plot_to = data_end)
        
        # Plot cumulative metrics for each canton (with future projection)
        fig_name = c("Baseline projection", canton, "Cumulative")
        plot_temporal(o, canton, fig_name, everything = TRUE, cumulative = TRUE)
      }
    }
  }
  
  # ---- Plot past evaluation ----
  
  # Check whether any past scenarios have been defined
  if (o$plot_past == TRUE && length(o$scenarios$past) > 0) {
    
    # Repeat for each canton
    for (canton in o$cantons) {
      
      message(" - Past evaluation: ", canton)
      
      # Past evaluation plot
      fig_name = c("Past evaluation", canton)
      plot_temporal(o, canton, fig_name, scenarios = o$scenarios$past, 
                    plot_to = data_end, alt_baseline = alt_baseline)
    }
  }
  
  # ---- Plot future scenario results ----
  
  # Check whether any future scenarios have been defined
  if (o$plot_future == TRUE && length(o$scenarios$future) > 0) {
    
    # Repeat for each canton
    for (canton in o$cantons) {
      
      message(" - Future scenarios: ", canton)
      
      # Future scenario projection
      fig_name = c("Future scenarios", canton)
      plot_temporal(o, canton, fig_name, scenarios = o$scenarios$future, 
                    alt_baseline = alt_baseline)
      
      # Future scenario projection - all metrics
      fig_name = c("Future scenarios", canton, "Everything")
      plot_temporal(o, canton, fig_name, scenarios = o$scenarios$future, 
                    everything = TRUE, alt_baseline = alt_baseline)
      
      if (o$plot_cumulative == TRUE) {
        
        # Future scenario projection - cumulative
        fig_name = c("Future scenarios", canton, "Cumulative")
        plot_temporal(o, canton, fig_name, scenarios = o$scenarios$future, 
                      cumulative = TRUE, alt_baseline = alt_baseline)
      }
    }
  }
  
  # ---- Plot strategies ----
  
  # Shorthand for all strategies modelled
  all_strategies = setdiff(o$strategies, "none")
  
  if (o$plot_strategies == TRUE) {
    
    # Custom metrics for strategy plots
    plot_metrics = c("osi_level", "total_vaccinated", "confirmed", 
                     "hospital_beds", "icu_beds", "deaths")
    
    # Repeat for each canton
    for (canton in o$cantons) {
      
      message(" - All strategies")
      
      # Plot all combinations in this strategy - default metrics relative to baseline
      fig_name = c("All strategies", canton)
      plot_temporal(o, canton, fig_name, strategies = all_strategies,
                    plot_baseline = FALSE, plot_metrics = plot_metrics, legend_rows = 3)
      
      # Plot impact of all combinations in this strategy - default metrics
      fig_name = c("Impact bars", canton)
      plot_impact(o, canton, fig_name, strategies = all_strategies, plot_baseline = FALSE)
      
      # Plot strategies one at a time
      #
      # NOTE: You can plot multiple strategies at once if you like, this is just for default output
      for (strategy in all_strategies) {
        
        message(" - Strategy: ", strategy)
        
        # Plot all combinations in this strategy - default metrics relative to baseline
        fig_name = c("Strategy", canton, strategy)
        plot_temporal(o, canton, fig_name, strategies = strategy, 
                      plot_baseline = FALSE, plot_metrics = plot_metrics)
        
        # Plot all combinations in this strategy - default metrics relative to baseline
        fig_name = c("Strategy", canton, strategy, "Everything")
        plot_temporal(o, canton, fig_name, strategies = strategy, 
                      plot_baseline = FALSE, everything = TRUE)
        
        # Plot age group for baseline calibration - area plot
        fig_name = c("Strategy", canton, strategy, "Age groups", "Area")
        plot_temporal(o, canton, fig_name, strategies = strategy, everything = TRUE,
                      plot_baseline = FALSE, plot_by = "age", plot_area = TRUE)
        
        # Plot virus variants for baseline calibration
        fig_name = c("Strategy", canton, strategy, "Variants")
        plot_temporal(o, canton, fig_name, strategies = strategy, everything = TRUE, 
                      plot_baseline = FALSE, plot_by = "variant")
        
        # Plot impact of all combinations in this strategy - default metrics
        fig_name = c("Impact bars", canton, strategy)
        plot_impact(o, canton, fig_name, strategies = strategy, plot_baseline = FALSE)
        
        # Do we also we to plot full factorial elements
        if (o$plot_elements == TRUE) {
          
          # Load strategy info and extract the different elements that comprise this strategy
          strategy_scenarios = try_load(o$pth$strategy, paste0("strategy.", strategy))$combinations
          elements = names(select(strategy_scenarios, -id, -name))
          
          # Loop through these (currently just 'roll_out' and 'npi_exit')
          for (element in elements) {
            
            # Loop through policies associated with this element
            element_policies = unique(strategy_scenarios[[element]])
            for (element_policy in element_policies) {
              
              # All scenarios associated with this policy
              policy_scenarios = filter(strategy_scenarios, get(element) == element_policy)$id
              
              # Default metrics for all scenarios associated with this policy
              fig_name = c("Strategy", canton, strategy, element_policy)
              plot_temporal(o, canton, fig_name, scenarios = policy_scenarios,
                            plot_baseline = FALSE, everything = TRUE)
              
              # Custom metrics for all scenarios associated with this policy
              fig_name = c("Strategy", canton, strategy, element_policy, "Custom metrics")
              plot_temporal(o, canton, fig_name, scenarios = policy_scenarios, 
                            plot_baseline = FALSE, plot_metrics = plot_metrics)
            }
          }
        }
      }
    }
  }
  
  # ---- Visualise data ----
  
  if (o$plot_data_sources == TRUE) {
    
    message(" - Data sources")
    
    # Reload data in order to produce plot (see load_data.R)
    load_data(o, plot_data_sources = TRUE, quiet = TRUE)
  }
  
  # ---- Visualise emulator performance ----
  
  if (o$plot_emulator == TRUE) {
    
    message(" - Emulator performance")
    
    # Plot emulator performance (train vs test)
    plot_emulator(o)
    
    # Plot emulator performance (resampling iterations)
    plot_emulator(o, by_resampling = TRUE)
  }
  
  # ---- Visualise parameters ----
  
  if (o$plot_parameters == TRUE) {
    
    message(" - Posteriors")
    
    # Plot canton-grouped posteriors vs priors
    plot_posteriors(o)
    
    # Parameter trace plots for all MCMC chains
    plot_posterior_traces(o)
    
    # Plot parameter-likelihood relationship
    plot_param_likelihood(o)  # Only really useful if testing one or two parameters
  }
  
  # ---- Output csv files ----
  
  if (o$output_csv == TRUE) {
    
    message(" - Creating csv files")
    
    # Output files for the media...
    
    # Seperate strategies into sets for ease of interpretation (and size)
    media_strategies = list(a = qc(s0a, s0b, s1a, s1b, s2a, s2b, s3a, s3b), 
                            b = qc(s8a, s8b, s9a, s9b, s10a, s10b))
    
    # Metrics to report (no disaggregation)
    media_metrics = c("osi_level", "total_vaccinated", "n_doses", 
                      "confirmed", "hospital_beds", "icu_beds")
    
    # Loop through results sets
    for (this_set in names(media_strategies)) {
      set_name = paste("Figure ", toupper(this_set))
      
      # Subset of strategies to produce output for
      this_strategies = media_strategies[[this_set]]
      
      # Produce CSV file of default metrics for baseline scenario
      file_name = c("OpenCOVID model output", set_name)
      output_csv(o, "CH", file_name, plot_baseline = FALSE,
                 strategies = this_strategies, plot_metrics = media_metrics)
    }
    
    # Special results for KOF: policy breif 23rd March 2021...
    
    # Seperate strategies into sets for ease of interpretation (and size)
    kof_strategies = list(a = qc(s0a, s0b, s1a, s1b, s2a, s2b, s3a, s3b), 
                          b = qc(s8a, s8b, s9a, s9b, s10a, s10b), 
                          c = qc(s0a, s0b, s1a, s1b, s2a, s2b, s3a, s3b)) # TODO: Reactive measures 
    
    # Metrics to report (no disaggregation)
    kof_metrics = c("osi_level")
    
    # Metrics to report (disaggregated by age)
    kof_metrics_age = c("confirmed", "new_local_infections", 
                        "currently_infected", "currently_symptomatic", 
                        "hospital_admissions", "hospital_beds", 
                        "icu_admissions", "icu_beds", 
                        "total_vaccinated", "n_doses")
    
    # Loop through results sets
    for (this_set in names(kof_strategies)) {
      set_name = paste("Set", toupper(this_set))
      
      # Subset of strategies to produce output for
      this_strategies = kof_strategies[[this_set]]
      
      # Fine granularity for age at death
      o$plot_ages = c(0, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85)
      
      # Produce CSV file of default metrics for baseline scenario
      file_name = c("KOF", "Deaths by age", set_name)
      output_csv(o, "CH", file_name, plot_baseline = FALSE, plot_by = "age", 
                 strategies = this_strategies, plot_metrics = "deaths")
      
      # Broader age grouping for other metrics
      o$plot_ages = c(0, 20, 65)
      
      # Produce CSV file of default metrics for baseline scenario
      file_name = c("KOF", "Other metrics by age", set_name)
      output_csv(o, "CH", file_name, plot_baseline = FALSE, plot_by = "age", 
                 strategies = this_strategies, plot_metrics = kof_metrics_age)
      
      # Produce CSV file of default metrics for baseline scenario
      file_name = c("KOF", "OSI level", set_name)
      output_csv(o, "CH", file_name, plot_baseline = FALSE, 
                 strategies = this_strategies, plot_metrics = kof_metrics)
    }
  }
  
  # ---- Manuscript plots ----
  
  # Produce manuscript figures (see manuscript.R)
  if (o$plot_fig2) fig2(o)
  if (o$plot_fig3) fig3(o)
  if (o$plot_fig4) fig4(o)
  if (o$plot_fig5) fig5(o)
  if (o$plot_fig6) fig6(o)
  if (o$plot_fig7) fig7(o)
  
  # Produce supplementary figures (see manuscript.R)
  if (o$plot_sup1) sup1(o)
}

