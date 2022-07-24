###########################################################
# MANUSCRIPT
#
# File used to generate figures for March 2021 manuscript.
#
###########################################################

# ---------------------------------------------------------
# Figure 2: Main results figure
# ---------------------------------------------------------
fig2 = function(o, quiet = FALSE, strategies = NULL, override_metrics = NULL, 
                plot_legend = TRUE) {
  
  if (!quiet) message(" - Manuscript figure 2")
  
  # Metrics to plot (set defaults which can be overridden)
  if (!is.null(override_metrics)) metrics = override_metrics
  else metrics = c("osi_level", "total_vaccinated", "n_doses", 
                   "confirmed", "icu_beds", "deaths")
  
  # Temporal plots from
  plot_from = "2020-01-01"
  
  # Strategies to plot (all for each vaccine rollout)
  if (is.null(strategies))
    strategies = c("s3x", "s1x", "s2x", "s0x")
  
  # Wrap strip names longer than n characters
  n_wrap = 28
  
  # ---- Set up dictionaries ----
  
  # Dictionary for vaccination rollouts
  vax_dict = fig_dict("vaccine")
  if (sum(grepl("x", strategies)) == 0)
    vax_dict = vax_dict[unique(str_sub(strategies, 3))]
  
  # Load NPI scenario dictionary and and reorder
  parent_strategies = paste0(str_sub(strategies, 1, 2), "x")
  scen_dict = fig_dict("scenario")[parent_strategies]
  
  # Complete list of scenarios we'll be plotting (strategies x vaccine rollouts)
  vax_strategies = rep(names(vax_dict), each = length(strategies))
  all_strategies = str_replace(strategies, "x", vax_strategies)
  
  # Expand scenario dictionary for vaccine rollouts
  all_scen_dict  = setNames(rep(unname(scen_dict), length(vax_dict)), 
                            str_replace(names(scen_dict), "x", vax_strategies))
  
  # ---- Extract results ----
  
  # Collate and interpret function inputs so we know what to plot
  f = figure_properties(o, strategies = all_strategies, 
                        plot_baseline = FALSE, 
                        plot_metrics = metrics, 
                        plot_from = plot_from, 
                        override_colours = NA)
  
  # Otherwise load baseline file for this canton
  a_baseline = readRDS(paste0(o$pth$scenario, "CH_", f$baseline_name, ".rds"))
  
  # Inherit key settings from analysis file (see store_options() function)
  o[names(a_baseline$options)] = a_baseline$options
  
  # Metric dictionary
  metric_dict  = f$metric_names
  metric_names = unname(f$metric_names)
  
  # Extract and format baseline data
  data_df = format_data(o, f, a_baseline) %>%
    mutate(metric = factor(metric_dict[metric], levels = metric_names))
  
  # Initiate plotting dataframe
  plot_list = list()
  
  # Loop through metrics to plot
  for (scenario in f$scenarios) {
    
    # Load data file for this canton for this scenario
    if (scenario == f$baseline_name) a = a_baseline else {
      a = try_load(o$pth$scenario, paste0("CH_", scenario)) 
    }
    
    # Format model output and store in list to be concatenated
    plot_list[[scenario]] = format_results(o, f, a)
  }
  
  # Parse strategy names for vaccine speed and variables tested 
  plot_df = rbindlist(plot_list) %>%
    mutate(scenario = str_split_fixed(scenario, "\\.", 2)[, 1], 
           vaccine  = unname(vax_dict[str_extract(scenario, "[a,b]+")]), 
           vaccine  = factor(vaccine, levels = unname(vax_dict)), 
           scenario = unname(all_scen_dict[scenario]), 
           scenario = factor(scenario, levels = rev(unname(scen_dict))))
  
  # Uncertainty bounds don't make sense in the context of vacce rollout
  plot_df[metric %in% c("n_doses", "total_vaccinated"), lower := NA]
  plot_df[metric %in% c("n_doses", "total_vaccinated"), upper := NA]
  
  # Plot historical dose data as points rather than jumpy continuos line - much clearer
  dots_df = plot_df %>%
    filter(metric == "n_doses" & date <= max(o$dates_data), value > 0) %>%
    mutate(metric = factor(metric_dict[metric], levels = metric_names))
  
  # Remove vaccination zeros prior to assumed start of vaccine rollout
  plot_df = plot_df %>%
    mutate(value  = ifelse(metric == "n_doses" & date <= max(o$dates_data), NA, value), 
           value  = ifelse(metric == "total_vaccinated" & date < o$vaccine_start, NA, value), 
           metric = factor(metric_dict[metric], levels = metric_names))
  
  # ---- Create figure ----
  
  # Initiate ggplot structure with each metric in it's facet
  g = ggplot(plot_df, aes(x = date, y = value)) + 
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = scenario), 
                linetype = 0, alpha = 0.3, show.legend = plot_legend) + 
    geom_line(aes(colour = scenario, linetype = vaccine), 
              show.legend = plot_legend, size = 1.2) + 
    geom_point(data = dots_df, colour = o$data_colour) + 
    facet_wrap(~metric, scales = "free", 
               labeller = as_labeller(label_wrap_gen(n_wrap)))
  
  # Add data points
  g = g + geom_point(data = data_df, colour = o$data_colour)
  
  # Add a vertical reference line where data end and prediction starts
  g = g + geom_vline(xintercept = max(o$dates_data), linetype = "dashed", colour = o$dash_colour)
  
  # Add dashed line for ICU capacity
  if ("icu_beds" %in% f$metrics) {
    icu_capacity = a_baseline$fit_data$capacity
    
    # Construct a dataframe for capacity + 10% so reference line is clearly visible
    capacity_blank = data.frame(date = as.Date(NA), group = NA, value  = icu_capacity * 1.1, 
                                metric = factor(metric_dict["icu_beds"], levels = metric_names))
    
    # Plot the blank data then a hline for the actual capacity
    g = g + geom_blank(data = capacity_blank) + 
      geom_hline(data = capacity_blank, aes(yintercept = icu_capacity), linetype = "dashed") +
      geom_hline(data = capacity_blank, aes(yintercept = icu_capacity * 0.25), 
                 colour = "darkred", linetype = "dotdash")
  }
  
  # ---- Figure aesthetics ----
  
  l_order  = c("s2x", "s3x", "s1x", "s0x")
  l_order  = intersect(l_order, parent_strategies)
  l_breaks = unname(scen_dict[l_order])
  
  # Apply scenario dict - done here so more obvious if names are not unique
  g = g + scale_fill_manual(values = fig_colours(l_order), breaks = l_breaks) + 
    scale_colour_manual(values = fig_colours(l_order), breaks = l_breaks) + 
    scale_linetype_manual(values = c("solid", "22"))
  
  # Prettify axes
  g = g + expand_limits(y = c(0, 100)) +  # Restrain upper y limit from being really low
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.05))) + 
    scale_x_date(date_breaks = o$x_tick_dates, date_labels = o$x_tick_format) 
  
  # Order legends and reverse NPI policies
  g = g + guides(colour   = guide_legend(byrow = TRUE, order = 1), 
                 fill     = guide_legend(byrow = TRUE, order = 1), 
                 linetype = guide_legend(byrow = TRUE, order = 2))
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(strip.text   = element_text(size = 14),
          axis.text.y  = element_text(size = 12),
          axis.text.x  = element_text(size = 12, hjust = 1, angle = 50),
          axis.title   = element_blank(),
          legend.title = element_blank(),
          legend.text  = element_text(size = 13),
          legend.position = "bottom", 
          legend.key.size = unit(2, "line"), 
          legend.box   = "vertical", 
          axis.line    = element_blank(), 
          panel.border = element_rect(colour = "black", fill = NA, size = 1))
  
  # Save figure
  if (!quiet) save_fig(o, g, "Figure 2")
  
  return(g)
}

# ---------------------------------------------------------
# Figure 3: Scenario impact bars
# ---------------------------------------------------------
fig3 = function(o, quiet = FALSE, override_metrics = NULL) {
  
  if (!quiet) message(" - Manuscript figure 3")
  
  # Strategies to plot (all for each vaccine rollout)
  strategies = c("s2x", "s3x", "s1x", "s0x")
  
  # Metrics to plot (set defaults which can be overridden)
  if (!is.null(override_metrics)) metrics = override_metrics
  else metrics = c("icu_admissions", "deaths", "confirmed", "n_doses")
  
  # Wrap strip names longer than n characters
  n_wrap = 28
  
  # Flag for plotting error bars 
  error_bars = FALSE
  
  # Transparency values for vaccine rollouts
  alpha_vals = c(0.85, 0.4)
  
  # ---- Set up dictionaries ----
  
  # Dictionary for vaccination rollouts
  vax_dict = fig_dict("vaccine")
  
  # Load NPI scenario dictionary and and reorder
  scen_dict = fig_dict("scenario")[strategies]
  
  # Complete list of scenarios we'll be plotting (strategies x vaccine rollouts)
  vax_strategies = rep(names(vax_dict), each = length(strategies))
  all_strategies = str_replace(strategies, "x", vax_strategies)
  
  # Expand scenario dictionary for vaccine rollouts
  all_scen_dict  = setNames(rep(unname(scen_dict), length(vax_dict)), 
                            str_replace(names(scen_dict), "x", vax_strategies))
  
  # ---- Extract model predictions ----
  
  # Collate and interpret function inputs so we know what to plot
  f = figure_properties(o, strategies = all_strategies, 
                        plot_baseline = FALSE, 
                        plot_metrics = metrics, 
                        plot_from = max(o$dates_data) + 1, 
                        override_colours = NA, 
                        cumulative = TRUE)
  
  # Cumulative needed for metric names, turn off for extracting data
  f$cumulative = FALSE
  
  # Load baseline (or alt baseline) file
  a_baseline = readRDS(paste0(o$pth$scenario, "CH_", f$baseline_name, ".rds"))
  
  # Inherit key settings from analysis file (see store_options() function)
  o[names(a_baseline$options)] = a_baseline$options
  
  # Loop through metrics to plot
  plot_list = list()
  for (scenario in f$scenarios) {
    
    # Load data file for this canton for this scenario
    if (scenario == f$baseline_name) a = a_baseline else {
      a = try_load(o$pth$scenario, paste0("CH_", scenario)) 
    }
    
    # Format model output and store in list to be concatenated
    plot_list[[scenario]] = format_results(o, f, a)
  }
  
  # Number of days of projection - used for OSI mean
  n_days = length(setdiff(o$dates_all, o$dates_data))
  
  # Sum each metric over time for each scenario
  plot_df = rbindlist(plot_list) %>%
    group_by(scenario, metric) %>%
    summarise(value = sum(value), 
              lower = sum(lower), 
              upper = sum(upper)) %>%
    mutate(scenario = str_split_fixed(scenario, "\\.", 2)[, 1], 
           vaccine  = unname(vax_dict[str_extract(scenario, "[a,b]+")]), 
           vaccine  = factor(vaccine, levels = unname(vax_dict)), 
           scenario = unname(all_scen_dict[scenario]), 
           scenario = factor(scenario, levels = unname(scen_dict)), 
           value    = ifelse(metric == "osi_level", value / n_days, value),
           metric   = factor(f$metric_names[metric], levels = unname(f$metric_names)))
  
  # Faceted bar chart
  g = ggplot(plot_df, aes(x = scenario, y = value, ymin = lower, ymax = upper, 
                          fill = scenario, alpha = vaccine)) + 
    geom_bar(stat = "identity", position = "dodge2", colour = "black", size = 0.5) + 
    facet_wrap(~metric, scales = "free_y", 
               labeller = as_labeller(label_wrap_gen(n_wrap)))
  
  # Apply error bars if desired
  if (error_bars == TRUE)
    g = g + geom_errorbar(colour = "darkgrey", width = 0.25, size = 0.5)
  
  # Set the colour and transparency scheme
  g = g + scale_alpha_manual(values = alpha_vals) + 
    scale_fill_manual(values = fig_colours(strategies)) +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.05)))
  
  # Apply scenario colour scheme and use descriptive scenario names
  g = g + guides(fill  = guide_legend(order = 1), 
                 alpha = guide_legend(order = 2))
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(strip.text   = element_text(size = 16),
          axis.text.y  = element_text(size = 12),
          axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title   = element_blank(), 
          legend.title = element_blank(),
          legend.text  = element_text(size = 15),
          legend.key.size = unit(2, "line"), 
          axis.line    = element_blank(), 
          panel.grid.major.y = element_line(size = 0.5, colour = "grey80"),
          panel.grid.minor.y = element_line(size = 0.5, colour = "grey80"),
          panel.border = element_rect(colour = "black", fill = NA, size = 1))
  
  # Save figure
  if (!quiet) save_fig(o, g, "Figure 3")
  
  return(g)
}

# ---------------------------------------------------------
# Figure 4: Alternative relax scenarios
# ---------------------------------------------------------
fig4 = function(o) {
  
  message(" - Manuscript figure 4")
  
  # Metrics to plot
  metrics = c("osi_level", "total_vaccinated", "n_doses", 
              "confirmed", "icu_beds", "deaths")
  
  # Temporal plots from
  plot_from = "2020-09-01"
  
  # Strategies to plot (all for each vaccine rollout)
  strategies = c("s1x", "s8x", "s9x", "s10x")
  
  # Wrap strip names longer than n characters
  n_wrap = 28
  
  # ---- Set up dictionaries ----
  
  # Dictionary for vaccination rollouts
  vax_dict = fig_dict("vaccine")
  
  # Load NPI scenario dictionary and and reorder
  scen_dict = fig_dict("scenario")[strategies]
  
  # Complete list of scenarios we'll be plotting (strategies x vaccine rollouts)
  vax_strategies = rep(names(vax_dict), each = length(strategies))
  all_strategies = str_replace(strategies, "x", vax_strategies)
  
  # Expand scenario dictionary for vaccine rollouts
  all_scen_dict  = setNames(rep(unname(scen_dict), length(vax_dict)), 
                            str_replace(names(scen_dict), "x", vax_strategies))
  
  # ---- Extract results ----
  
  # Collate and interpret function inputs so we know what to plot
  f = figure_properties(o, strategies = all_strategies, 
                        plot_baseline = FALSE, 
                        plot_metrics = metrics, 
                        plot_from = plot_from, 
                        override_colours = NA)
  
  # Otherwise load baseline file for this canton
  a_baseline = readRDS(paste0(o$pth$scenario, "CH_", f$baseline_name, ".rds"))
  
  # Inherit key settings from analysis file (see store_options() function)
  o[names(a_baseline$options)] = a_baseline$options
  
  # Metric dictionary
  metric_dict  = f$metric_names
  metric_names = unname(f$metric_names)
  
  # Extract and format baseline data
  data_df = format_data(o, f, a_baseline) %>%
    mutate(metric = factor(metric_dict[metric], levels = metric_names))
  
  # Initiate plotting dataframe
  plot_list = list()
  
  # Loop through metrics to plot
  for (scenario in f$scenarios) {
    
    # Load data file for this canton for this scenario
    if (scenario == f$baseline_name) a = a_baseline else {
      a = try_load(o$pth$scenario, paste0("CH_", scenario)) 
    }
    
    # Format model output and store in list to be concatenated
    plot_list[[scenario]] = format_results(o, f, a)
  }
  
  # Parse strategy names for vaccine speed and variables tested 
  plot_df = rbindlist(plot_list) %>%
    mutate(scenario = str_split_fixed(scenario, "\\.", 2)[, 1], 
           vaccine  = unname(vax_dict[str_extract(scenario, "[a,b]+")]), 
           vaccine  = factor(vaccine, levels = unname(vax_dict)), 
           scenario = unname(all_scen_dict[scenario]), 
           scenario = factor(scenario, levels = rev(unname(scen_dict))))
  
  # Uncertainty bounds don't make sense in the context of vacce rollout
  plot_df[metric %in% c("n_doses", "total_vaccinated"), lower := NA]
  plot_df[metric %in% c("n_doses", "total_vaccinated"), upper := NA]
  
  # Plot historical dose data as points rather than jumpy continuos line - much clearer
  dots_df = plot_df %>%
    filter(metric == "n_doses" & date <= max(o$dates_data), value > 0) %>%
    mutate(metric = factor(metric_dict[metric], levels = metric_names))
  
  # Remove vaccination zeros prior to assumed start of vaccine rollout
  plot_df = plot_df %>%
    mutate(value  = ifelse(metric == "n_doses" & date <= max(o$dates_data), NA, value), 
           value  = ifelse(metric == "total_vaccinated" & date < o$vaccine_start, NA, value), 
           metric = factor(metric_dict[metric], levels = metric_names))
  
  # ---- Create figure ----
  
  # Initiate ggplot structure with each metric in it's facet
  g = ggplot(plot_df, aes(x = date, y = value)) + 
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = scenario), 
                linetype = 0, alpha = 0.3) + 
    geom_line(aes(colour = scenario, linetype = vaccine), size = 1.2) + 
    geom_point(data = dots_df, colour = o$data_colour) + 
    facet_wrap(~metric, scales = "free", 
               labeller = as_labeller(label_wrap_gen(n_wrap)))
  
  # Add data points
  g = g + geom_point(data = data_df, colour = o$data_colour)
  
  # Add a vertical reference line where data end and prediction starts
  g = g + geom_vline(xintercept = max(o$dates_data), linetype = "dashed", colour = o$dash_colour)
  
  # Add dashed line for ICU capacity
  if ("icu_beds" %in% f$metrics) {
    icu_capacity = a_baseline$fit_data$capacity
    
    # Construct a dataframe for capacity + 10% so reference line is clearly visible
    capacity_blank = data.frame(date = as.Date(NA), group = NA, value = 0, # icu_capacity, 
                                metric = factor(metric_dict["icu_beds"], levels = metric_names))
    
    # Plot the blank data then a hline for the actual capacity
    g = g + geom_blank(data = capacity_blank) + 
      # geom_hline(data = capacity_blank, aes(yintercept = icu_capacity), linetype = "dashed") +
      geom_hline(data = capacity_blank, aes(yintercept = icu_capacity * 0.25), 
                 colour = "darkred", linetype = "dotdash")
  }
  
  # ---- Figure aesthetics ----
  
  # Construct colour vector
  colours = c(o$baseline_colour, fig_colours(strategies[-1]))
  
  # Order legend
  legend_breaks = unname(scen_dict[strategies])
  
  # Apply scenario dict - done here so more obvious if names are not unique
  g = g + scale_fill_manual(values = colours, breaks = legend_breaks) + 
    scale_colour_manual(values = colours, breaks = legend_breaks) + 
    scale_linetype_manual(values = c("solid", "22"))
  
  # Prettify axes
  g = g + expand_limits(y = c(0, 100)) +  # Restrain upper y limit from being really low
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.05))) + 
    scale_x_date(date_breaks = o$x_tick_dates, date_labels = o$x_tick_format) 
  
  # Order legends and reverse NPI policies
  g = g + guides(colour   = guide_legend(byrow = TRUE, order = 1), 
                 fill     = guide_legend(byrow = TRUE, order = 1), 
                 linetype = guide_legend(byrow = TRUE, order = 2))
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(strip.text   = element_text(size = 14),
          axis.text.y  = element_text(size = 12),
          axis.text.x  = element_text(size = 12, hjust = 1, angle = 50),
          axis.title   = element_blank(),
          legend.title = element_blank(),
          legend.text  = element_text(size = 13),
          legend.position = "bottom", 
          legend.key.size = unit(2, "line"), 
          legend.box   = "vertical", 
          axis.line    = element_blank(), 
          panel.border = element_rect(colour = "black", fill = NA, size = 1))
  
  # Save figure
  save_fig(o, g, "Figure 4")
}

# ---------------------------------------------------------
# Figure 5: Alternative scenario impact bars
# ---------------------------------------------------------
fig5 = function(o) {
  
  message(" - Manuscript figure 5")
  
  # Strategies to plot (all for each vaccine rollout)
  strategies = c("s1x", "s8x", "s9x", "s10x")
  
  # Metrics to plot
  metrics = c("icu_admissions", "deaths", "confirmed", "n_doses")
  
  # Wrap strip names longer than n characters
  n_wrap = 28
  
  # Flag for plotting error bars 
  error_bars = FALSE
  
  # Transparency values for vaccine rollouts
  alpha_vals = c(0.85, 0.4)
  
  # ---- Set up dictionaries ----
  
  # Dictionary for vaccination rollouts
  vax_dict = fig_dict("vaccine")
  
  # Load NPI scenario dictionary and and reorder
  scen_dict = fig_dict("scenario")[strategies]
  
  # Complete list of scenarios we'll be plotting (strategies x vaccine rollouts)
  vax_strategies = rep(names(vax_dict), each = length(strategies))
  all_strategies = str_replace(strategies, "x", vax_strategies)
  
  # Expand scenario dictionary for vaccine rollouts
  all_scen_dict  = setNames(rep(unname(scen_dict), length(vax_dict)), 
                            str_replace(names(scen_dict), "x", vax_strategies))
  
  # ---- Extract model predictions ----
  
  # Collate and interpret function inputs so we know what to plot
  f = figure_properties(o, strategies = all_strategies, 
                        plot_baseline = FALSE, 
                        plot_metrics = metrics, 
                        plot_from = max(o$dates_data) + 1, 
                        override_colours = NA, 
                        cumulative = TRUE)
  
  # Cumulative needed for metric names, turn off for extracting data
  f$cumulative = FALSE
  
  # Load baseline (or alt baseline) file
  a_baseline = readRDS(paste0(o$pth$scenario, "CH_", f$baseline_name, ".rds"))
  
  # Inherit key settings from analysis file (see store_options() function)
  o[names(a_baseline$options)] = a_baseline$options
  
  # Loop through metrics to plot
  plot_list = list()
  for (scenario in f$scenarios) {
    
    # Load data file for this canton for this scenario
    if (scenario == f$baseline_name) a = a_baseline else {
      a = try_load(o$pth$scenario, paste0("CH_", scenario)) 
    }
    
    # Format model output and store in list to be concatenated
    plot_list[[scenario]] = format_results(o, f, a)
  }
  
  # Sum each metric over time for each scenario
  plot_df = rbindlist(plot_list) %>%
    group_by(scenario, metric) %>%
    summarise(value = sum(value), 
              lower = sum(lower), 
              upper = sum(upper)) %>%
    mutate(scenario = str_split_fixed(scenario, "\\.", 2)[, 1], 
           vaccine  = unname(vax_dict[str_extract(scenario, "[a,b]+")]), 
           vaccine  = factor(vaccine, levels = unname(vax_dict)), 
           scenario = unname(all_scen_dict[scenario]), 
           scenario = factor(scenario, levels = unname(scen_dict)), 
           metric   = factor(f$metric_names[metric], levels = unname(f$metric_names)))
  
  # Faceted bar chart
  g = ggplot(plot_df, aes(x = scenario, y = value, ymin = lower, ymax = upper, 
                          fill = scenario, alpha = vaccine)) + 
    geom_bar(stat = "identity", position = "dodge2", colour = "black", size = 0.5) + 
    facet_wrap(~metric, scales = "free_y", 
               labeller = as_labeller(label_wrap_gen(n_wrap)))
  
  # Apply error bars if desired
  if (error_bars == TRUE)
    g = g + geom_errorbar(colour = "darkgrey", width = 0.25, size = 0.5)
  
  # Construct colour vector
  colours = c(o$baseline_colour, fig_colours(strategies[-1]))
  
  # Set the colour and transparency scheme
  g = g + scale_alpha_manual(values = alpha_vals) + 
    scale_fill_manual(values = colours) +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.05)))
  
  # Apply scenario colour scheme and use descriptive scenario names
  g = g + guides(fill  = guide_legend(order = 1), 
                 alpha = guide_legend(order = 2))
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(strip.text   = element_text(size = 16),
          axis.text.y  = element_text(size = 12),
          axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title   = element_blank(), 
          legend.title = element_blank(),
          legend.text  = element_text(size = 15),
          legend.key.size = unit(2, "line"), 
          axis.line    = element_blank(), 
          panel.grid.major.y = element_line(size = 0.5, colour = "grey80"),
          panel.grid.minor.y = element_line(size = 0.5, colour = "grey80"),
          panel.border = element_rect(colour = "black", fill = NA, size = 1))
  
  # Save figure
  save_fig(o, g, "Figure 5")
}

# ---------------------------------------------------------
# Figure 6: Temporal plots of key uncertainties
# ---------------------------------------------------------
fig6 = function(o) {
  
  message(" - Manuscript figure 6")
  
  # Metrics to plot
  metrics = c("deaths", "confirmed")
  
  # Temporal plots from
  plot_from = "2020-09-01"
  
  # Strategies to plot (baseline first, all for each vaccine rollout)
  strategies = c("s3x", "s5x1", "s5x2", "s6x2", "s6x1", "s4x2", "s4x1", "s7x")
  
  # Wrap strip names longer than n characters
  n_wrap = 25
  
  # ---- Set up dictionaries ----
  
  # Dictionary for vaccination rollouts
  vax_dict = fig_dict("vaccine")
  
  # Load dictionaries and and reorder
  scen_dict = fig_dict("scenario")[strategies]
  var_dict  = fig_dict("variable")[unique(str_extract(strategies, "s[0-9]+"))]
  
  # Complete list of scenarios we'll be plotting (strategies x vaccine rollouts)
  vax_strategies = rep(names(vax_dict), each = length(strategies))
  all_strategies = str_replace(strategies, "x", vax_strategies)
  
  # Expand scenario dictionary for vaccine rollouts
  all_scen_dict  = setNames(rep(unname(scen_dict), length(vax_dict)), 
                            str_replace(names(scen_dict), "x", vax_strategies))
  
  # ---- Extract results ----
  
  # Collate and interpret function inputs so we know what to plot
  f = figure_properties(o, strategies = all_strategies, 
                        plot_baseline = FALSE, 
                        plot_metrics  = metrics, 
                        plot_from = plot_from, 
                        override_colours = NA)
  
  # Otherwise load baseline file for this canton
  a_baseline = readRDS(paste0(o$pth$scenario, "CH_", f$baseline_name, ".rds"))
  
  # Inherit key settings from analysis file (see store_options() function)
  o[names(a_baseline$options)] = a_baseline$options
  
  # Extract and format baseline data
  data_df = format_data(o, f, a_baseline) %>%
    mutate(metric = factor(f$metric_names[metric], levels = unname(f$metric_names)))
  
  # Loop through metrics to plot
  plot_list = list()
  for (scenario in f$scenarios) {
    
    # Load data file for this canton for this scenario
    if (scenario == f$baseline_name) a = a_baseline else {
      a = try_load(o$pth$scenario, paste0("CH_", scenario)) 
    }
    
    # Format model output and store in list to be concatenated
    plot_list[[scenario]] = format_results(o, f, a)
  }
  
  # Parse strategy names for vaccine speed and variables tested 
  plot_df = rbindlist(plot_list) %>%
    mutate(scenario = str_split_fixed(scenario, "\\.", 2)[, 1], 
           variable = unname(var_dict[str_extract(scenario, "s[0-9]+")]), 
           variable = factor(variable, levels = unname(var_dict)), 
           vaccine  = unname(vax_dict[str_extract(scenario, "[a,b]+")]), 
           vaccine  = factor(vaccine, levels = unname(vax_dict)), 
           scenario = unname(all_scen_dict[scenario]), 
           scenario = factor(scenario, levels = unname(scen_dict)), 
           metric   = factor(f$metric_names[metric], levels = unname(f$metric_names)))
  
  # Split into baseline and sensitivity variables dataframes
  var_df  = filter(plot_df, variable != unname(var_dict[1]))
  base_df = filter(plot_df, variable == unname(var_dict[1])) %>%
    select(-scenario, -variable)
  
  # ---- Produce plot ----
  
  # Plot baseline in each facet, then seperate by metric, vaccine, and variables
  g = ggplot(var_df, aes(x = date, y = value)) + 
    geom_ribbon(data = base_df, aes(ymin = lower, ymax = upper),
                colour = o$baseline_colour, linetype = 0, alpha = 0.3, show.legend = FALSE) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = scenario), 
                linetype = 0, alpha = 0.3) + 
    geom_line(data = base_df, colour = o$baseline_colour, size = 1.2, show.legend = FALSE) +
    geom_line(aes(colour = scenario), size = 1.2) + 
    facet_grid(metric + vaccine ~ variable, scales = "free_y", 
               labeller = as_labeller(label_wrap_gen(n_wrap)))
  
  # Add data points
  g = g + geom_point(data = data_df, colour = o$data_colour)
  
  # Set colour scheme
  g = g + scale_colour_manual(values = fig_colours(strategies[-1])) +
    scale_fill_manual(values = fig_colours(strategies[-1]))
  
  # Order legends and reverse NPI policies
  g = g + guides(colour = guide_legend(nrow = 2), 
                 fill   = guide_legend(nrow = 2))
  
  # Prettify axes
  g = g + expand_limits(y = c(0, 100)) +  # Restrain upper y limit from being really low
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.05))) + 
    scale_x_date(date_breaks = o$x_tick_dates, date_labels = o$x_tick_format) 
  
  # Prettify theme
  g = g + theme_classic() +
    theme(strip.text   = element_text(size = 12), 
          axis.text.y  = element_text(size = 11),
          axis.text.x  = element_text(size = 11, hjust = 1, angle = 50),
          axis.title   = element_blank(),
          axis.line    = element_blank(),
          legend.title = element_blank(),
          legend.text  = element_text(size = 12),
          legend.key.size = unit(2, "line"),
          legend.position = "bottom", 
          panel.border = element_rect(colour = "black", fill = NA, size = 1))
  
  # Save figure
  save_fig(o, g, "Figure 6")
}

# ---------------------------------------------------------
# Figure 7: Tornado plots of key uncertainties
# ---------------------------------------------------------
fig7 = function(o) {
  
  message(" - Manuscript figure 7")
  
  # Metrics to plot
  metrics = c("deaths", "confirmed")
  
  # Temporal plots from
  plot_from = "2020-09-01"
  
  # Strategies to plot (baseline first, all for each vaccine rollout)
  strategies = c("s3x", "s4x2", "s4x1", "s7x", "s5x1", "s5x2", "s6x2", "s6x1")
  
  # COntrol order of bars in tornado plot
  bar_order = c("s4", "s7", "s5", "s6", "s3")
  
  # ---- Set up dictionaries ----
  
  # Dictionary for vaccination rollouts
  vax_dict = fig_dict("vaccine")
  
  # Load dictionaries and and reorder
  scen_dict = fig_dict("scenario")[strategies]
  var_dict  = fig_dict("variable")[rev(bar_order)]
  
  # Reset name of baseline scenario
  scen_dict[1] = var_dict[[1]]
  
  # Complete list of scenarios we'll be plotting (strategies x vaccine rollouts)
  vax_strategies = rep(names(vax_dict), each = length(strategies))
  all_strategies = str_replace(strategies, "x", vax_strategies)
  
  # Expand scenario dictionary for vaccine rollouts
  all_scen_dict  = setNames(rep(unname(scen_dict), length(vax_dict)), 
                            str_replace(names(scen_dict), "x", vax_strategies))
  
  # ---- Extract results ----
  
  # Collate and interpret function inputs so we know what to plot
  f = figure_properties(o, strategies = all_strategies, 
                        plot_baseline = FALSE, 
                        plot_metrics  = metrics, 
                        plot_from = plot_from, 
                        override_colours = NA)
  
  # Otherwise load baseline file for this canton
  a_baseline = readRDS(paste0(o$pth$scenario, "CH_", f$baseline_name, ".rds"))
  
  # Inherit key settings from analysis file (see store_options() function)
  o[names(a_baseline$options)] = a_baseline$options
  
  # Loop through metrics to plot
  plot_list = list()
  for (scenario in f$scenarios) {
    
    # Load data file for this canton for this scenario
    if (scenario == f$baseline_name) a = a_baseline else {
      a = try_load(o$pth$scenario, paste0("CH_", scenario)) 
    }
    
    # Format model output and store in list to be concatenated
    plot_list[[scenario]] = format_results(o, f, a)
  }
  
  # Parse strategy names for vaccine speed and variables tested 
  plot_df = rbindlist(plot_list) %>%
    mutate(scenario = str_split_fixed(scenario, "\\.", 2)[, 1], 
           variable = unname(var_dict[str_extract(scenario, "s[0-9]+")]), 
           variable = factor(variable, levels = unname(var_dict)), 
           vaccine  = unname(vax_dict[str_extract(scenario, "[a,b]+")]), 
           vaccine  = factor(vaccine, levels = unname(vax_dict)), 
           scenario = unname(all_scen_dict[scenario]), 
           scenario = factor(scenario, levels = unname(scen_dict)), 
           metric   = factor(f$metric_names[metric], levels = unname(f$metric_names)))
  
  bars_df = plot_df %>%
    group_by(metric, scenario, variable, vaccine) %>%
    summarise(total = sum(value)) %>%
    group_by(metric, vaccine) %>%
    mutate(diff = total - total[variable == var_dict[[1]]]) %>%
    filter(scenario != var_dict[[1]])
  
  blank_df = mutate(bars_df, diff = -diff)
  
  # ---- Produce plot ----
  
  g = ggplot(bars_df, aes(x = variable, y = diff, fill = scenario)) +
    geom_bar(position = "identity", stat = "identity", colour = "black", size = 0.5) +
    geom_blank(data = blank_df) + geom_hline(yintercept = 0) + coord_flip() + 
    facet_grid(vaccine ~ metric, scales = "free")
  
  g = g + scale_fill_manual(values = fig_colours(strategies[-1]))
  
  g = g + scale_y_continuous(name = "Difference relative to best estimate values", labels = comma)
  
  # Prettify theme
  g = g + theme_classic() +
    theme(strip.text   = element_text(size = 16), 
          axis.title = element_blank(),
          axis.text.x  = element_text(size = 12),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          legend.title = element_blank(),
          legend.text  = element_text(size = 14),
          legend.key.size = unit(2, "line"),
          axis.line    = element_blank(),
          panel.grid.major.x = element_line(size = 0.5, colour = "grey80"),
          panel.border = element_rect(colour = "black", fill = NA, size = 1))
  
  # Save figure
  save_fig(o, g, "Figure 7")
}

# ---------------------------------------------------------
# Supplement 1: Population susceptibity over time
# ---------------------------------------------------------
sup1 = function(o) {
  
  message(" - Supplement figure 1")
  
  # Metric to plot
  plot_metric = "pop_susceptibility"
  
  # Temporal plots from
  plot_from = "2020-01-01"
  
  # Strategies to plot (all for each vaccine rollout)
  strategies = c("s3x", "s1x", "s2x", "s0x")
  
  # Flag to plot percentage rather than absolute numbers
  plot_percentage = TRUE
  
  # ---- Set up dictionaries ----
  
  # Dictionary for vaccination rollouts
  vax_dict = fig_dict("vaccine")
  
  # Load NPI scenario dictionary and and reorder
  scen_dict = fig_dict("scenario")[strategies]
  
  # Complete list of scenarios we'll be plotting (strategies x vaccine rollouts)
  vax_strategies = rep(names(vax_dict), each = length(strategies))
  all_strategies = str_replace(strategies, "x", vax_strategies)
  
  # Expand scenario dictionary for vaccine rollouts
  all_scen_dict  = setNames(rep(unname(scen_dict), length(vax_dict)), 
                            str_replace(names(scen_dict), "x", vax_strategies))
  
  # ---- Extract results ----
  
  # Collate and interpret function inputs so we know what to plot
  f = figure_properties(o, strategies = all_strategies, 
                        plot_baseline = FALSE, 
                        plot_metrics = plot_metric, 
                        plot_from = plot_from, 
                        override_colours = NA)
  
  # Otherwise load baseline file for this canton
  a_baseline = readRDS(paste0(o$pth$scenario, "CH_", f$baseline_name, ".rds"))
  
  # Inherit key settings from analysis file (see store_options() function)
  o[names(a_baseline$options)] = a_baseline$options
  
  # Initiate plotting dataframe
  plot_list = list()
  
  # Loop through metrics to plot
  for (scenario in f$scenarios) {
    
    # Load data file for this canton for this scenario
    if (scenario == f$baseline_name) a = a_baseline else {
      a = try_load(o$pth$scenario, paste0("CH_", scenario)) 
    }
    
    # Format model output and store in list to be concatenated
    plot_list[[scenario]] = format_results(o, f, a)
  }
  
  # Parse strategy names for vaccine speed and variables tested 
  plot_df = rbindlist(plot_list) %>%
    mutate(scenario = str_split_fixed(scenario, "\\.", 2)[, 1], 
           vaccine  = unname(vax_dict[str_extract(scenario, "[a,b]+")]), 
           vaccine  = factor(vaccine, levels = unname(vax_dict)), 
           scenario = unname(all_scen_dict[scenario]), 
           scenario = factor(scenario, levels = rev(unname(scen_dict))))
  
  # We can either plot absolute (as is) or convert to pop percentage
  if (plot_percentage == TRUE) {
    
    # Total number of people in the population
    n_pop = sum(a_baseline$fit_data$demog)
    
    # Convert best estimate
    plot_df$value = 100 * plot_df$value / n_pop
    
    # Convert bounds
    plot_df$lower = 100 * plot_df$lower / n_pop
    plot_df$upper = 100 * plot_df$upper / n_pop
  }
  
  # ---- Create figure ----
  
  # Initiate ggplot structure with each metric in it's facet
  g = ggplot(plot_df, aes(x = date, y = value)) + 
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = scenario), 
                linetype = 0, alpha = 0.3) + 
    geom_line(aes(colour = scenario), size = 1.2) + 
    facet_wrap(~vaccine, scales = "free")
  
  # Add a vertical reference line where data end and prediction starts
  g = g + geom_vline(xintercept = max(o$dates_data), linetype = "dashed", colour = o$dash_colour)
  
  # ---- Figure aesthetics ----
  
  l_order  = c("s2x", "s3x", "s1x", "s0x")
  l_breaks = unname(scen_dict[l_order])
  
  # Apply scenario dict - done here so more obvious if names are not unique
  g = g + scale_fill_manual(values = fig_colours(l_order), breaks = l_breaks) + 
    scale_colour_manual(values = fig_colours(l_order), breaks = l_breaks) + 
    scale_linetype_manual(values = c("solid", "22"))
  
  # Construct a figure title
  fig_title = f$metric_names
  if (plot_percentage)
    fig_title = paste(fig_title, "(%)")
  
  # Add the figure title
  g = g + ggtitle(fig_title)
  
  # Prettify axes
  g = g + expand_limits(y = c(0, 100)) +  # Restrain upper y limit from being really low
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0))) + 
    scale_x_date(date_breaks = o$x_tick_dates, date_labels = o$x_tick_format) 
  
  # Order legends and reverse NPI policies
  g = g + guides(colour   = guide_legend(byrow = TRUE, order = 1), 
                 fill     = guide_legend(byrow = TRUE, order = 1), 
                 linetype = guide_legend(byrow = TRUE, order = 2))
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(plot.title   = element_text(size = 32, hjust = 0.5),
          strip.text   = element_text(size = 20),
          axis.text.y  = element_text(size = 13),
          axis.text.x  = element_text(size = 13, hjust = 1, angle = 50),
          axis.title   = element_blank(),
          panel.grid.major.y = element_line(size = 0.5, colour = "grey80"),
          panel.grid.minor.y = element_line(size = 0.5, colour = "grey80"),
          legend.title = element_blank(),
          legend.text  = element_text(size = 14),
          legend.position = "bottom", 
          legend.key.size = unit(2, "line"), 
          legend.box   = "vertical", 
          axis.line    = element_blank(), 
          panel.border = element_rect(colour = "black", fill = NA, size = 1))
  
  # Save figure
  save_fig(o, g, "Population susceptibility")
}

# ---------------------------------------------------------
# Figure dictionaries
# ---------------------------------------------------------
fig_dict = function(dict_type) {
  
  # Dictionary for vaccination rollout speeds
  if (dict_type == "vaccine")
    dict = c(a = "Slower vaccine rollout", 
             b = "Faster vaccine rollout")

  # Scenario dictionary for legend
  if (dict_type == "scenario")
    dict = c(s0x  = "No relax after 1 March", 
             s1x  = "Relax 22 March", 
             s2x  = "Phased relax until 3 May", 
             s3x  = "Phased relax until 5 July", 
             s4x1 = "B.1.1.7 50% more infectious than D614G",
             s4x2 = "B.1.1.7 70% more infectious than D614G",
             s5x1 = "mRNA 60% transmission blocking",
             s5x2 = "mRNA 95% transmission blocking",
             s6x1 = "Vaccine acceptance 90% (P2-P5)", 
             s6x2 = "Vaccine acceptance 60% (P2-P5)", 
             s7x  = "B.1.1.7 50% higher mortality than D614G", 
             s8x  = "Relax 12 April", 
             s9x  = "Relax 3 May", 
             s10x = "Relax 24 May")
  
  # Variable dictionary for sensitivity facets
  if (dict_type == "variable")
    dict = c(s3 = "Baseline", 
             s4 = "VOC infectiousness",
             s5 = "Vaccine transmission blocking",
             s6 = "Vaccine acceptibility", 
             s7 = "VOC mortality")
  
  return(dict)
}

# ---------------------------------------------------------
# Figure colours
# ---------------------------------------------------------
fig_colours = function(scenarios) {
  
  # Colours for all figures
  all_colours = c(s0x  = "springgreen4", 
                  s1x  = "gold1", 
                  s2x  = "firebrick3", 
                  s3x  = "dodgerblue3", 
                  s4x1 = "mediumpurple1", 
                  s4x2 = "mediumpurple4",
                  s5x1 = "darkorange4", 
                  s5x2 = "darkorange1", 
                  s6x1 = "springgreen1", 
                  s6x2 = "springgreen4", 
                  s7x  = "deeppink3", 
                  s8x  = "tan1", 
                  s9x  = "deepskyblue1", 
                  s10x = "deeppink1")
  
  # Subset and convert to character vector
  colours = unname(all_colours[scenarios])
  
  return(colours)
}

