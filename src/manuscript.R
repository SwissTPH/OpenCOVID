###########################################################
# MANUSCRIPT FIGURES
#
# Produce figures for air pollution paper.
#
###########################################################

# ---------------------------------------------------------
# Main plot-calling function
# ---------------------------------------------------------
manuscript_figures = function(o) {
  
  message(" - Manuscript figures")
  
  # Toggle plotting flags
  #
  # NOTE: Figure 1 is a methodology figure
  do_fig2 = TRUE  # Cumulative bars
  do_fig3 = TRUE  # Temporal trends: infections
  do_fig4 = TRUE  # Temporal trends: hospitalisations
  do_fig5 = TRUE  # Temporal trends: ICU
  do_fig6 = TRUE  # Temporal trends: deaths
  do_sup1 = TRUE  # Air pollution effect assumptions
  do_sup2 = TRUE  # All baseline metrics
  do_sup3 = TRUE  # Cumulative comparison to baseline: infections
  do_sup4 = TRUE  # Cumulative comparison to baseline: hospitalisations
  do_sup5 = TRUE  # Cumulative comparison to baseline: ICU
  do_sup6 = TRUE  # Cumulative comparison to baseline: deaths
  
  # Toggles for diagnostic plots
  do_time = FALSE
  
  # ---- Plot properties ----
  
  # Initiate list of figure properties
  f = list()
  
  # Burn in period to cut from start of forward simulations
  f$burn_in = 14  # 2 weeks
  
  # Dictionary of metrics to plot
  f$metric_dict = list(
    time  = c(all_new_infections  = "New infections (per 100k people)", 
              hospital_beds       = "Hospital occupancy (per 100k people)", 
              icu_beds            = "ICU occupancy (per 100k people)",
              deaths              = "Deaths (per 100k people)"),
    total = c(all_new_infections  = "New infections (per 100k/year)", 
              hospital_admissions = "Hospital admissions (per 100k/year)", 
              icu_admissions      = "ICU admissions (per 100k/year)",
              deaths              = "Deaths (per 100k/year)"))
  
  # Dictionary for type pollution-susceptibility/severity relationship
  f$relationship_dict = c(log10  = "Log10 increase in effect", 
                          linear = "Linear increase in effect")
  
  # Dictionary for type of air pollution 'effect'
  f$effect_dict = c(susceptibility = "Effect on susceptibility", 
                    severity       = "Effect on severity")
  
  # Set legend titles
  f$legend_relationship = "Pollution-effect relationship"
  f$legend_effect       = "Mechanism of air pollution effect"
  
  f$legend_relationship = "Pollution-effect relationship: "
  f$legend_effect       = "Mechanism of air pollution effect: "
  
  # Colours (one for each 'effect')
  f$colours = list(
    metric = c("red1", "dodgerblue1", "gold1"), 
    levels = "viridis::viridis",
    base   = "grey70")
  
  # Trasparency of shaded area between bounds
  f$transparency = 0.7
  
  # Summary statistic used for 'best estimate' line
  f$stat = "mean"  # OPTIONS: "mean" or median
  
  # ---- Load results ----
  
  message("  > Loading results")
  
  # All scenarios we've simulated
  scenarios    = parse_yaml(o, "*read*")
  scenarios_id = names(scenarios)
  
  # Initiate results list
  results = list(time = list(), total = list())
  
  # Formatting code in a function
  format_result = function(f, result, type1, type2) {
    
    # Filter burn in and plot metrics
    df = result[[type1]] %>%
      filter(metric %in% names(f$metric_dict[[type2]]),
             date > f$burn_in, 
             grouping == "none") %>%
      mutate(date = date - f$burn_in) %>%
      select(scenario, metric, date,
             best = !!f$stat, lower, upper)
  }
  
  # Loop through the scenarios
  for (scenario in scenarios_id) {
    
    # Load the scenario file 
    result = try_load(o$pth$scenarios, scenario)

    # Extract both temporal and cumulative results
    results$total[[scenario]] = format_result(f, result, "cum_output", "total")
    results$time[[scenario]]  = format_result(f, result, "output",     "time")
  }
  
  # Load model inputs (for baseline scenario)
  input = try_load(o$pth$scenarios, "baseline")$input
  
  # Scaler to report metrics per 100k per year
  scale_pop  = 1e5 / input$population_size
  scale_time = 365 / input$n_days
  
  # ---- Construct general plotting dataframe ----
  
  message("  > Formatting results")
  
  # Initiate list to store plotting dataframes
  plot_list = list()
  
  # Iterate through output formats
  for (type in names(results)) {
    
    # Recode metric names as defined above
    results_df = rbindlist(results[[type]]) %>%
      mutate(metric = recode(metric, !!!f$metric_dict[[type]]))
    
    # Baseline scenario repeated
    baseline_df = results_df %>%
      filter(scenario == "baseline") %>%
      expand_grid(relationship = unique(f$relationship_dict), 
                  effect       = unique(f$effect_dict)) %>%
      mutate(pollution = "baseline", 
             x         = 0) %>%
      select(-scenario) %>%
      setDT()
    
    # All alternative scenarios
    scenario_df = results_df %>%
      filter(scenario != "baseline") %>%
      mutate(scenario = recode(scenario, !!!scenarios)) %>%
      separate(col  = scenario, 
               into = qc(pollution, relationship, effect), 
               sep  = " ") %>%
      mutate(relationship = recode(relationship, !!!f$relationship_dict), 
             effect       = recode(effect,       !!!f$effect_dict)) %>%
      mutate(x = str_extract(pollution, "[0-9]{1}"),
             x = ifelse(is.na(x), 0, as.numeric(x)))
    
    # Concatenate these two datatables
    plot_list[[type]] = bind_rows(baseline_df, scenario_df) %>%
      group_by(metric, date, pollution, effect) %>%
      mutate(y_min = min(lower), 
             y_max = max(upper)) %>%  # Bounds stretch across effects
      ungroup() %>%
      mutate(relationship = factor(relationship, f$relationship_dict), 
             effect       = factor(effect,       f$effect_dict)) %>%
      arrange(metric, x, relationship, effect) %>%
      setDT()
  }
  
  # Extract unique pollution levels (to retain ordering)
  pollution_levels = unique(plot_list$total$pollution)
  
  # ---- Figure 2 ----
  
  # Check plotting flag
  if (do_fig2) {
    fig_name = "Figure 2"
    
    message("  > ", fig_name)
    
    # Scale results to per 100k per year
    plot_df = plot_list$total %>%
      filter(date == max(date)) %>%
      mutate(metric = factor(metric, f$metric_dict$total), 
             across(.cols = c(best, lower, upper, y_min, y_max),
                    .fns  = ~ . * scale_pop * scale_time)) %>%
      select(x, best, y_min, y_max, effect, relationship, metric)
    
    # Produce basic plot
    g = ggplot(plot_df) + 
      aes(x = x, 
          y = best, 
          ymin = y_min, 
          ymax = y_max,
          colour   = effect,  
          fill     = effect, 
          linetype = relationship) +
      geom_ribbon(linetype = 0, alpha = 1 - f$transparency) +
      geom_line(size = 2) + 
      facet_wrap(~metric, scales = "free_y")
    
    # Apply colour scheme
    g = g + scale_fill_manual(values = f$colours$metric) +
      scale_colour_manual(values = f$colours$metric)
    
    # Prettify y-axis
    g = g + scale_y_continuous(labels = comma, 
                               limits = c(0, NA), 
                               expand = expansion(mult = c(0, 0.05)))
    
    # Prettify x-axis
    g = g + scale_x_continuous(labels = pollution_levels, 
                               expand = expansion(mult = c(0, 0)))
    
    # Construct legend titles
    legend_colour = guide_legend(title = f$legend_effect)
    legend_line   = guide_legend(title = f$legend_relationship)
    
    # Prettify legend
    g = g + guides(fill     = legend_colour, 
                   colour   = legend_colour, 
                   linetype = legend_line)
    
    # Prettify theme
    g = g + theme_classic() + 
      theme(strip.text    = element_text(size = 20),
            axis.text     = element_text(size = 12),
            axis.text.x   = element_text(hjust = 1, angle = 50), 
            axis.title    = element_blank(),
            axis.line     = element_blank(),
            panel.border  = element_rect(size = 1, colour = "black", fill = NA),
            panel.spacing = unit(1, "lines"),
            legend.text   = element_text(size = 12),
            legend.title  = element_text(size = 14),
            legend.key.height = unit(2, "lines"),
            legend.key.width  = unit(2, "lines"),
            legend.box.background = element_rect(), 
            strip.background      = element_blank())
    
    # Save the plot to file
    fig_save(o, g, fig_name)
  }
  
  # ---- Figures 3-6 ----
  
  # Check plotting flags
  metrics = f$metric_dict$time[c(do_fig3, do_fig4, do_fig5, do_fig6)]
  
  # Loop through metrics to produce plot for
  for (metric in metrics) {
    f$metric = metric
    
    # Construct figure name
    fig_number = which(metric == f$metric_dict$time) + 2
    fig_name   = paste0("Figure ", fig_number)
    
    message("  > ", fig_name)
    
    # Temporal values for this metric - scale to 100k people
    plot_df = plot_list$time %>%
      filter(metric == !!metric) %>%
      mutate(pollution = fct_rev(fct_inorder(pollution)), 
             across(.cols = c(best, lower, upper),
                    .fns  = ~ . * scale_pop)) %>%
      select(date, best, lower, upper, pollution, effect, relationship)
    
    # Call specialist plotting function
    manuscript_temporal(o, f, plot_df, fig_name)
  }
  
  # ---- Supplementary figure 1 ----
  
  # Check plotting flag
  if (do_sup1) {
    fig_name = "Supplementary 1"
    
    message("  > ", fig_name)
    
    # Load a scenario with non-trivial air pollution - susceptibility effect
    sus_scenario = scenarios_id[grep("sus$", scenarios_id)[1]]
    sus_input    = try_load(o$pth$scenarios, sus_scenario)$input
    
    # Load a scenario with non-trivial air pollution - severity effect
    sev_scenario = scenarios_id[grep("sev$", scenarios_id)[1]]
    sev_input    = try_load(o$pth$scenarios, sev_scenario)$input
    
    # Store these two model inputs 
    pollution_factors = list(susceptibility = sus_input$air_pollution$susceptibility, 
                             severity       = sev_input$air_pollution$severity)
    
    # Use this input to plot out air pollution 
    plot_pollution_relationships(o, fig_name, pollution_factors)
  }
  
  # ---- Supplementary figure 2 ----
  
  # Check plotting flag
  if (do_sup2) {
    fig_name = "Supplementary 2"
    
    message("  > ", fig_name)
    
    # Basleine metrics over time: vanilla usage of plot_temporal
    plot_temporal(o, fig_name, plot_from = f$burn_in + 1)
  }
  
  # ---- Supplementary figures 3-6 ----
  
  # Check plotting flags
  metrics = f$metric_dict$total[c(do_sup3, do_sup4, do_sup5, do_sup6)]
  
  # Loop through metrics to produce plot for
  for (metric in metrics) {
    f$metric = metric
    
    # Call specialist plotting function
    fig_number = which(metric == f$metric_dict$total) + 2
    fig_name   = paste0("Supplementary ", fig_number)
    
    message("  > ", fig_name)
    
    # Cumulative difference with baseline - scale to 100k people
    plot_df = plot_list$total %>%
      filter(metric == !!metric) %>%
      group_by(date, relationship, effect) %>%
      mutate(best  = best  - best[pollution  == "baseline"],
             lower = lower - lower[pollution == "baseline"],
             upper = upper - upper[pollution == "baseline"]) %>%
      ungroup() %>%
      filter(pollution != "baseline") %>%
      mutate(across(.cols = c(best, lower, upper),
                    .fns  = ~ . * scale_pop)) %>%
      select(date, best, lower, upper, pollution, effect, relationship) %>%
      setDT()
    
    # Call specialist plotting function
    manuscript_temporal(o, f, plot_df, fig_name)
  }
  
  # ---- Figure: simulation time ----
  
  # Check plotting flag
  if (do_time) {
    fig_name = "Simulation time"
    
    message("  > ", fig_name)
    
    # Use a custom colour scheme (as we have 20+ scenarios)
    colours = colour_scheme("pals::kovesi.rainbow", 
                            n = length(scenarios))
    
    # Plot simulation time for each scenario - simply for diagnostic purposes
    plot_simulation_time(o, fig_name, 
                         scenarios = scenarios_id, 
                         override_colours = colours)
  }
}

# ---------------------------------------------------------
# Temporal plotting function - used multiple times
# ---------------------------------------------------------
manuscript_temporal = function(o, f, plot_df, fig_name) {
  
  # Produce basic plot
  g = ggplot(plot_df) + 
    aes(x = date, 
        y = best, 
        ymin = lower, 
        ymax = upper,
        colour   = pollution,  
        fill     = pollution) +
    geom_ribbon(linetype = 0, alpha = 1 - f$transparency) +
    geom_line(size = 2) + 
    facet_grid(effect~relationship)
  
  # Determine whether baseline is included
  all_levels   = unique(plot_df$pollution)
  inc_baseline = "baseline" %in% all_levels
  
  # Generate colours
  n_colours   = length(all_levels) - sum(inc_baseline)
  set_colours = colour_scheme(f$colours$levels, n = n_colours)
  
  # Append baseline colour if apprpriate
  if (inc_baseline)
    set_colours = c(f$colours$base, set_colours)
  
  # Apply colour scheme
  g = g + scale_fill_manual(values = rev(set_colours)) +
    scale_colour_manual(values = rev(set_colours))
  
  # Prettify axxs
  g = g + 
    scale_y_continuous(labels = comma, 
                       expand = expansion(mult = c(0.05, 0.05))) + 
    scale_x_continuous(expand = expansion(mult = c(0, 0)))
  
  # Set a title
  g = g + labs(title    = f$metric, 
               subtitle = "(cumulative value relative to baseline)", 
               x        = "Day of simulation")
  
  # Prettify legend
  g = g + guides(fill   = guide_legend(reverse = TRUE), 
                 colour = guide_legend(reverse = TRUE))
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(plot.title    = element_text(size = 32, hjust = 0.5),
          plot.subtitle = element_text(size = 22, hjust = 0.5),
          strip.text    = element_text(size = 20),
          axis.text     = element_text(size = 12),
          axis.text.x   = element_text(hjust = 1, angle = 50), 
          axis.title.x  = element_text(size = 24),
          axis.title.y  = element_blank(),
          axis.line     = element_blank(),
          panel.border  = element_rect(size = 1, colour = "black", fill = NA),
          panel.spacing = unit(1, "lines"),
          legend.text   = element_text(size = 12),
          legend.title  = element_blank(), # size = 14
          legend.key.height = unit(2, "lines"),
          legend.key.width  = unit(2, "lines"),
          legend.box.background = element_rect(), 
          strip.background      = element_blank())
  
  # Save the plot to file
  fig_save(o, g, fig_name)
}

