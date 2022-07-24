###########################################################
# PLOTTING
#
# All plotting functions in one place.
#
###########################################################

# ---------------------------------------------------------
# Plot metrics over time for multiple groups, scenarios, or cantons 
# ---------------------------------------------------------
plot_temporal = function(o, canton, fig_name, plot_file = NULL, likelihood = NA, plot_area = FALSE, ...) {
  
  # Collate and interpret function inputs so we know what to plot
  f = figure_properties(o, ...)
  
  # Do not allow relative difference to baseline for temporal plotting
  if (f$relative == TRUE)
    stop("Set 'relative' to FALSE for temporal plots (use impact bars to assess difference to baseline)")
  
  # ---- Load baseline (or alt baseline) file ----
  
  # Handle alternative functionality case - file fed in directly
  if (!is.null(plot_file)) a_baseline = plot_file
  else {
    
    # Otherwise load baseline file for this canton
    a_baseline = readRDS(paste0(o$pth$scenario, canton, "_", f$baseline_name, ".rds"))
    
    # Inherit key settings from analysis file (see store_options() function)
    o[names(a_baseline$options)] = a_baseline$options
  }
  
  # Extract and format baseline data
  data_df = format_data(o, f, a_baseline)
  
  # ---- Extract model predictions ----
  
  # Initiate plotting dataframe
  plot_list = list()
  
  # Loop through metrics to plot
  for (scenario in f$scenarios) {
    
    # Load data file for this canton for this scenario
    if (scenario == f$baseline_name) a = a_baseline else {
      a = try_load(o$pth$scenario, paste0(canton, "_", scenario)) 
    }
    
    # Format model output and store in list to be concatenated
    plot_list[[scenario]] = format_results(o, f, a)
  }
  
  # Concatenate plotting dataframes for all scenarios
  plot_df = rbindlist(plot_list) %>%
    mutate(scenario = factor(scenario, levels = rev(f$scenarios)))
  
  # ---- Create figure ----
  
  # Initiate ggplot structure with each metric in it's facet
  g = ggplot(plot_df, aes_string(x = "date", y = "value", colour = f$plot_type, fill = f$plot_type)) + 
    facet_wrap(~metric, scales = "free", nrow = f$facet_rows, labeller = as_labeller(f$metric_names))
  
  # Either plot by area (useful when plotting by groups)...
  if (plot_area == TRUE) g = g + geom_area() else {
    
    # ... or plot uncertainty bounds and best estimate
    g = g + geom_ribbon(aes(ymin = lower, ymax = upper), linetype = 0, alpha = 0.3) + 
      geom_line(size = 2)
  }
  
  # Check whether there is data to plot
  if (!is.null(data_df)) {
    
    # Add modelled data points to plot
    if (f$plot_type != "group")
      g = g + geom_point(data = data_df, colour = o$data_colour)
    
    # Use group colours (with black edges) if plotting data by group
    if (f$plot_type == "group")
      g = g + geom_point(data = data_df) + geom_point(data = data_df, shape = 1, colour = "black")
  }
  
  # ---- Add reference lines ----
  
  # Add a vertical reference line where data end and prediction starts
  g = g + geom_vline(xintercept = max(o$dates_data), linetype = "dashed", colour = o$dash_colour)
  
  # Add a vertical line at 1 if plotting effective reproduction number
  if ("R_effective" %in% f$metrics)
    g = g + geom_hline(data = subset(plot_df, metric == "R_effective"), aes(yintercept = 1))
  
  # Add dashed line for ICU capacity
  if ("icu_beds" %in% f$metrics) {
    icu_capacity = a_baseline$fit_data$capacity
    
    # Construct a dataframe for capacity + 10% so reference line is clearly visible
    capacity_blank = data.frame(date = as.Date(NA), group = NA, value  = icu_capacity * 1.1, 
                                metric = factor("icu_beds", levels = f$metrics))
    
    # Plot the blank data then a hline for the actual capacity
    g = g + geom_blank(data = capacity_blank) + 
      geom_hline(data = capacity_blank, aes(yintercept = icu_capacity), linetype = "dashed") + 
      geom_hline(data = capacity_blank, aes(yintercept = icu_capacity * 0.25), 
                 colour = "darkred", linetype = "dotdash")
  }
  
  # ---- Figure aesthetics ----
  
  # Set pretty colour scheme (either by metric, group, scenario, or canton)
  if (f$plot_type != "scenario")
    g = g + scale_colour_manual(values = f$colours) + scale_fill_manual(values = f$colours)
  
  # A little more needed when plotting multiple scenarios
  if (f$plot_type == "scenario") {
    
    # Apply scenario dict - done here so more obvious if names are not unique
    g = g + scale_colour_manual(values = rev(f$colours), labels = f$scenario_dict) + 
      scale_fill_manual(values = rev(f$colours), labels = f$scenario_dict)
    
    # Reverse legend order so baseline is at the top
    apply_guide = guide_legend(reverse = TRUE, nrow = f$legend_rows, byrow = TRUE)
    g = g + guides(colour = apply_guide, fill = apply_guide)
  }
  
  # Figure title: canton name (include likelihood if provided)
  f$fig_title = as.character(o$canton_dict[canton])
  if (!is.na(likelihood))
    f$fig_title = paste0(f$fig_title, " (likelihood = ", thou_sep(round(likelihood, 0)), ")")
  
  # Apply this title
  g = g + ggtitle(f$fig_title)
  
  # Set x axes to user-defined limits (see options.R)
  g = g + expand_limits(y = c(0, f$expand_limit)) +  # Restrain upper y limit from being really low
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.05))) + 
    scale_x_date(date_breaks = o$x_tick_dates, date_labels = o$x_tick_format)
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(plot.title   = element_text(size = f$fontsize[1], hjust = 0.5),
          strip.text   = element_text(size = f$fontsize[2]),
          axis.text.y  = element_text(size = f$fontsize[3]),
          axis.text.x  = element_text(size = f$fontsize[3], hjust = 1, angle = 50),
          axis.title   = element_blank(),
          legend.text  = element_text(size = f$fontsize[4]),
          legend.title = element_blank(),
          legend.key   = element_blank(),
          legend.position = "bottom", 
          legend.key.height = unit(2, "lines"),
          legend.box.background = element_rect())
  
  # Remove legend entirely if only plotting multiple metrics (the facets are enough)
  if (f$plot_type == "metric")
    g = g + theme(legend.position = "none")
  
  # Do not save figure if input is trivial  
  if (!is.null(fig_name))
    save_fig(o, g, fig_name)
  
  return(g)
}  

# ---------------------------------------------------------
# Plot metrics over time for multiple groups, scenarios, or cantons 
# ---------------------------------------------------------
plot_impact = function(o, canton, fig_name, error_bars = TRUE, ...) {
  
  # Collate and interpret function inputs so we know what to plot
  f = figure_properties(o, ...)
  
  # ---- Sanity check on inputs provided ----
  
  # These bars are cumulative by default - don't allow summing of cumulative metrics
  if (f$cumulative == TRUE)
    stop("Impact bars are cumulative by definition - set 'cumulative' to FALSE to avoid double counting")
  
  # As standard do not consider the past - cut further up if user requests
  f$plot_from = max(f$plot_from, max(o$dates_data))
  
  # Make sure we have some future dates to work with
  if (f$plot_to < f$plot_from)
    stop("Impact bars are designed to show future impact - set 'plot_to' to a date in the future")
  
  # These bars are cumulative by default - don't allow summing of cumulative metrics
  if (f$plot_type == "group")
    stop("Impact bar plots have not been properly tested for plotting groupings yet")
  
  # ---- Extract model predictions ----
  
  # Load baseline (or alt baseline) file
  a_baseline = readRDS(paste0(o$pth$scenario, canton, "_", f$baseline_name, ".rds"))
  
  # Inherit key settings from analysis file (see store_options() function)
  o[names(a_baseline$options)] = a_baseline$options
  
  # Initiate plotting list
  plot_list = list()
  
  # Loop through metrics to plot
  for (scenario in f$scenarios) {
    
    # Load data file for this canton for this scenario
    if (scenario == f$baseline_name) a = a_baseline else {
      a = try_load(o$pth$scenario, paste0(canton, "_", scenario)) 
    }
    
    # Format model output and store in list to be concatenated
    plot_list[[scenario]] = format_results(o, f, a)
  }
  
  # Sum each metric over time for each scenario
  plot_df = rbindlist(plot_list) %>%
    mutate(scenario = factor(scenario, levels = f$scenarios)) %>%
    group_by(scenario, metric) %>%
    summarise(value = sum(value), 
              lower = sum(lower), 
              upper = sum(upper)) %>%
    as.data.table()
  
  # Check flag for plotting difference relative to baseline
  if (f$relative == TRUE) {
    
    # Difference betwen outcomes of each scenario and baseline
    plot_df = group_by(plot_df, metric) %>%
      mutate(value = value - value[scenario == f$baseline_name]) %>%
      filter(scenario != f$baseline_name) %>%
      as.data.table()
  }
  
  # ---- Produce bar charts of total over time ----
  
  # Multiple scenarios
  if (f$plot_type == "scenario") {
    
    # Faceted bar chart
    g = ggplot(plot_df, aes(x = scenario, y = value, ymin = lower, ymax = upper, fill = scenario)) + 
      geom_bar(stat = "identity", show.legend = FALSE) + 
      facet_wrap(~metric, scales = "free_y", labeller = as_labeller(f$metric_names))
    
    # Apply error bars if desired
    if (error_bars == TRUE)
      g = g + geom_errorbar(colour = "darkgrey", width = 0.25, size = 0.5)
    
    # Apply scenario colour scheme and use descriptive scenario names
    g = g + scale_fill_manual(values = f$colours) + 
      scale_x_discrete(labels = f$scenario_dict) + 
      scale_y_continuous(labels = comma)
  }
  
  # One single scenario (perhaps relative to baseline)
  if (f$plot_type == "metric") {
    
    # TODO: Work out the error bars...
    
    # Simple non-faceted bar chart
    g = ggplot(plot_df, aes(x = metric, y = value, ymin = lower, ymax = upper, fill = metric)) + 
      geom_bar(stat = "identity", position = position_dodge(), show.legend = FALSE)
    # geom_errorbar(position = position_dodge(width = 0.9), 
    #               colour = "darkgrey", width = 0.25, size = 0.5)
    
    # We've gone below zero (ie negative impact), so plot y = 0 reference line
    if (any(plot_df$value < 0))
      g = g + geom_hline(yintercept = 0)
    
    # Apply metric colour scheme and use descriptive scenario names
    g = g + scale_fill_manual(values = f$colours) + 
      scale_x_discrete(labels = f$metric_names) + 
      scale_y_continuous(labels = comma)
  }
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(strip.text  = element_text(size = f$fontsize[2]),
          axis.text.y = element_text(size = f$fontsize[3]),
          axis.text.x = element_text(size = f$fontsize[3], hjust = 1, angle = 50),
          axis.title  = element_blank())
  
  # Save these figures to file
  save_fig(o, g, fig_name)
  
  return(g)
} 

# ---------------------------------------------------------
# Plot normalised disease state of population over time
# ---------------------------------------------------------
plot_disease_state = function(o, p, states, run_time) {
  
  # Convert list to dataframe and scale population
  state_df = as.data.frame(states[o$all_states]) * p$population_scaler
  
  # Remove trivial dates and colunms we're not interested in
  state_df = state_df %>%
    mutate(date = o$dates_all[1 : run_time]) %>%
    filter(susc > 0) %>%
    select(-susc, -none)
  
  # Convert to long format and preserve order with factors
  plot_df = state_df %>%
    pivot_longer(cols      = -date,
                 names_to  = "state", 
                 values_to = "value") %>%
    mutate(state_type = ifelse(state %in% o$care_states, "care", "disease"), 
           state_type = factor(state_type, levels = c("disease", "care")), 
           state      = factor(state,      levels = o$all_states)) %>%
    arrange(state_type, state, date) %>%
    as.data.frame()
  
  # Plot as stacked area plot
  g1 = ggplot(plot_df, aes(x = date, y = value)) +
    geom_area(aes(fill = state)) +
    facet_wrap(~state_type, scales = "free_y") + 
    theme_classic()
  
  # Plot as stacked area plot
  g2 = ggplot(plot_df, aes(x = date, y = value)) +
    geom_line(aes(colour = state)) +
    facet_wrap(~state_type, scales = "free_y") + 
    theme_classic()
  
  # Save ggplot figure to file
  save_fig(o, g1, "Population disease and care states - stacked area")
  save_fig(o, g2, "Population disease and care states - line")
}

# ---------------------------------------------------------
# Plot series of network-related properties
# ---------------------------------------------------------
plot_network_properties = function(o, p, plot_contacts) {
  
  # ---- Figure properties ----
  
  # Bin width for age groups
  age_bin = 10 
  
  # Quantiles to display on violin plots
  violin_quantiles = c(0.5)  # eg c(0.25, 0.5, 0.75) for median and IQR
  
  # Set font sizes (axes labels, axes text, legend text)
  fontsize = c(16, 12, 13)
  
  # ---- Construct plotting dataframe ----
  
  # Number of contacts per person 
  count_contact = table(plot_contacts$from)
  
  # Seperate into IDs and number of contacts
  id = as.numeric(names(unlist(count_contact)))
  n_contacts = as.numeric(unlist(count_contact))
  
  # Breaks to create age bins
  age_breaks = seq(0, max(o$all_ages) + age_bin, by = age_bin)
  
  # Construct age dataframe to be joined to contact count
  age_df = unique(select(plot_contacts, from, age = age1)) %>%
    mutate(age_group = cut(age, age_breaks, include.lowest = TRUE))
  
  # Join number of contacts to age (and age group) of each person
  plot_df = data.table(from = id, n_contacts = n_contacts) %>%
    full_join(age_df, by = "from") %>%
    arrange(from)
  
  # Also construct dataframe of total number in each age group
  demog_df = data.table(age = o$all_ages, n = unname(p$demography)) %>%
    mutate(age_group = cut(age, age_breaks, include.lowest = TRUE)) %>%
    group_by(age_group) %>%
    summarise(n = sum(n)) %>%
    as.data.table()
  
  # ---- Fig a) Number of contacts per person ----
  
  # Construct subplot
  g1 = ggplot(plot_df, aes(x = n_contacts)) + 
    geom_density(aes(y = ..count.., colour = age_group, fill = age_group), 
                 adjust = 2, alpha = 0.1, size = 1.5) + 
    geom_vline(xintercept = p$contacts, size = 1.5, 
               colour = "darkgrey", linetype = "dashed")
  
  # Set axes labels
  g1 = g1 + scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.05))) + 
    xlab("Number of contacts per person") + 
    ylab("Number of people by age group")
  
  # Prettify subplot theme
  g1 = g1 + theme_classic() + 
    theme(axis.title = element_text(size = fontsize[1]),
          axis.text  = element_text(size = fontsize[2]), 
          legend.position = "none")
  
  # ---- Fig b) Number of contacts per person (normalised by age group) ----
  
  # Construct subplot
  g2 = ggplot(plot_df, aes(x = age_group, y = n_contacts, fill = age_group)) + 
    geom_violin(adjust = 2, draw_quantiles = violin_quantiles)
  
  # Set axes labels
  g2 = g2 + guides(fill = guide_legend(nrow = 1)) + 
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
    xlab("Normalised by age group") + 
    ylab("Number of contacts per person")
  
  # Prettify subplot theme (we take legend from this subplot)
  g2 = g2 + theme_classic() + 
    theme(axis.title = element_text(size = fontsize[1]),
          axis.text  = element_text(size = fontsize[2]),
          legend.position = "bottom", 
          legend.title = element_blank(), 
          legend.text  = element_text(size = fontsize[3]),
          legend.key   = element_blank(),
          legend.key.height = unit(1.5, "lines"),
          legend.box.background = element_rect())
  
  # ---- Fig c) Number of contacts per person (normalised by age group) ----
  
  # Construct subplot
  g3 = ggplot(plot_df, aes(x = n_contacts)) + 
    geom_density(aes(y = ..count.., fill = age_group, colour = age_group), 
                 position = "stack", adjust = 2) + 
    geom_density(aes(y = ..count..), adjust = 2, size = 1.2) + 
    geom_vline(xintercept = p$contacts, size = 1.5, 
               colour = "darkgrey", linetype = "dashed")
  
  # Prettify subplot axes
  g3 = g3 + scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.05))) + 
    xlab("Number of contacts per person") + 
    ylab("Total number of people")
  
  # Prettify subplot theme
  g3 = g3 + theme_classic() + 
    theme(axis.title = element_text(size = fontsize[1]),
          axis.text  = element_text(size = fontsize[2]), 
          legend.position = "none")
  
  # ---- Fig d) Number of contacts per person (normalised by age group) ----
  
  # Construct subplot
  g4 = ggplot(demog_df, aes(x = age_group, y = n, fill = age_group)) + 
    geom_bar(stat = "identity", color = "black")
  
  # Prettify subplot axes
  g4 = g4 + xlab("Age group") + ylab("Number of people") + 
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.05)))
  
  # Prettify subplot theme
  g4 = g4 + theme_classic() + 
    theme(axis.title = element_text(size = fontsize[1]),
          axis.text  = element_text(size = fontsize[2]), 
          legend.position = "none")
  
  # ---- Put it all together ----
  
  # Extract legend from subplot b
  g_legend = ggpubr::get_legend(g2)
  
  # Arrange all subplots into a figure with the legend
  g = ggarrange(g1, g2, g3, g4, nrow = 2, ncol = 2, 
                legend.grob = g_legend, legend = "bottom")
  
  # Save figure
  save_fig(o, g, "Network properties")
}

# ---------------------------------------------------------
# Age matrix of contact density per age (single age bins)
# ---------------------------------------------------------
plot_network_matrix = function(o, plot_contacts) {
  
  # ---- Figure properties ----
  
  # Bin width for age groups
  age_bins = c(1, 2, 4)
  
  # Set font sizes (subplot titles and axes text)
  fontsize = c(16, 12)
  
  # ---- Construct figure ----
  
  # Initiate list of ggplots
  plot_list = list()
  
  # Loop through age bins
  for (i in 1 : length(age_bins)) {
    age_bin = age_bins[i]
    
    # Breaks to create age bins
    age_breaks = seq(0, max(o$all_ages) + age_bin, by = age_bin)
    
    # Drop bottom half (as contacts are symmetric) and group by age
    format_df = plot_contacts %>%
      dplyr::slice_head(prop = 0.5) %>%  
      filter(age1 > 0, age2 > 0) %>%
      mutate(age_group1 = cut(age1, age_breaks, include.lowest = TRUE), 
             age_group2 = cut(age2, age_breaks, include.lowest = TRUE), 
             age_idx1 = as.numeric(age_group1), 
             age_idx2 = as.numeric(age_group2), 
             real_age1 = age_breaks[age_idx1] + age_bin / 2, 
             real_age2 = age_breaks[age_idx2] + age_bin / 2)
    
    # Count pairwise contacts by age
    count_df = suppressWarnings(
      widyr::pairwise_count(format_df, age_idx1, age_idx2, diag = TRUE)
    )
    
    # Add some context by calculating by age and in total
    plot_df = count_df %>%
      mutate(age1 = age_breaks[item1] + age_bin / 2, 
             age2 = age_breaks[item2] + age_bin / 2) %>%
      arrange(age1, age2) %>%
      group_by(age1) %>%
      mutate(n_age1 = sum(n)) %>%  # Number of contacts in age group1
      group_by(age2) %>%
      mutate(n_age2 = sum(n)) %>%  # Number of contacts in age group2
      ungroup() %>%
      mutate(p = n / (n_age1 * n_age2)) %>%
      select(-item1, -item2, -n_age1, -n_age2) %>%
      as.data.table()
    
    # Subplot a) Absolute number of contacts by age group
    g1 = ggplot(plot_df, aes(x = age1, y = age2)) + 
      geom_tile(aes(fill = n), show.legend = FALSE) + 
      labs(title = paste0("Number of contacts (bin width = ", age_bin, ")")) + 
      scale_fill_gradient(low = "grey70", high = "dodgerblue", expand = c(0, 0)) + 
      scale_x_continuous(expand = c(0, 0)) + 
      scale_y_continuous(expand = c(0, 0))
    
    # Prettify subplot theme
    g1 = g1 + theme_classic() + 
      theme(plot.title = element_text(size = fontsize[1], hjust = 0.5), 
            axis.title = element_blank(), 
            axis.text  = element_text(size = fontsize[2]))
    
    # Subplot b) Normalised number of contacts by age group
    g2 = ggplot(plot_df, aes(x = age1, y = age2)) + 
      geom_tile(aes(fill = p), show.legend = FALSE) + 
      labs(title = paste0("Normalised by age group (bin width = ", age_bin, ")")) + 
      scale_fill_gradient(low = "grey70", high = "springgreen") + 
      scale_x_continuous(expand = c(0, 0)) + 
      scale_y_continuous(expand = c(0, 0))
    
    # Prettify subplot theme
    g2 = g2 + theme_classic() + 
      theme(plot.title = element_text(size = fontsize[1], hjust = 0.5), 
            axis.title = element_blank(), 
            axis.text  = element_text(size = fontsize[2]))
    
    # Store subplots in plot list
    plot_list[[i]] = g1
    plot_list[[i + length(age_bins)]] = g2
  }
  
  # Arrange all subplots in grid - no need for legend colour bars
  g = plot_grid(plotlist = plot_list, align = "hv", nrow = 2)
  
  # Save figure
  save_fig(o, g, "Network matrix")
}

# ---------------------------------------------------------
# Viral load profile since day of infection
# ---------------------------------------------------------
plot_viral_load = function(o, p) {
  
  # Number of days post infection to plot up until
  plot_days = 30
  
  # Descriptive names for post-infection durations
  period_names = c("Mean latency period", 
                   "Mean presymptomatic period", 
                   "Mean duration of non-severe disease", 
                   "Mean duration of severe disease")
  
  # Colour for viral load curve
  vl_colour = "grey40"
  
  # Number of time to sample latency durations
  n_latency_samples = 100
  
  # Number of decimal places
  n_dp = 1
  
  # ---- Reconstruct viral load curve and sample latency durations ----
  
  # Vector of time points - more granular to cope with decimal places of latency best estimate
  time = seq(0, plot_days, by = 10 ^ -n_dp)
  
  # Construct best estimate, normalised viral load curve over time
  #
  # TODO: Add gamma parameters to model_parameters and remove hardcoding
  vl = dgamma(time, 3, rate = 0.5)
  vl = vl / max(vl)
  
  # Sample latency durations
  latency_samples = p$duration_latency(rep(0, n_latency_samples))
  
  # Number of time points to offset duations by
  offset = list(best  = 10 ^ n_dp * p$latency_days, 
                lower = 10 ^ n_dp * min(latency_samples), 
                upper = 10 ^ n_dp * max(latency_samples))
  
  # Constuct viral load plotting dataframe
  plot_df = data.frame(x  = time, 
                       y  = c(rep(0, offset$best),  head(vl, -offset$best)), 
                       y1 = c(rep(0, offset$lower), head(vl, -offset$lower)), 
                       y2 = c(rep(0, offset$upper), head(vl, -offset$upper))) %>%
    mutate(lower = pmin(y1, y2), upper = pmax(y1, y2))
  
  # Identify peaks to fill between
  vl_peaks = which(1 - plot_df$upper < 1e-6)
  plot_df$upper[min(vl_peaks) : max(vl_peaks)] = 1
  
  # Vector of durations - bound below by zero and above by plot_days
  durations = c(p$latency_days, p$presymptomatic_days, p$infectious_days_mild)
  
  # Construct background duration plotting dataframe
  rect_df = data.frame(period = factor(period_names, levels = period_names), 
                       from = c(0, cumsum(durations)), 
                       to   = c(cumsum(durations), plot_days))
  
  # ---- Produce plot ----
  
  # Plot viral load over rectangles showing durations - use default colours for rectangles
  g = ggplot(plot_df, aes(x = x, y = y)) + 
    geom_rect(data = rect_df, aes(NULL, NULL, xmin = from, xmax = to, fill = period),
              ymin = 0,ymax = 1.05, colour = "white", size = 0.4, alpha = 0.2) +
    geom_ribbon(aes(ymin = lower, ymax = upper), 
                fill = vl_colour, alpha = 0.5, linetype = 0) + 
    geom_line(colour = vl_colour, size = 3)
  
  # Set figure title and axis titles
  g = g + ggtitle("Viral load profile") + 
    ylab("Infectiousness multiplicative factor") + 
    xlab("Days since infection")
  
  # Set x axes to user-defined limits (see options.R)
  g = g + scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
    scale_x_continuous(expand = expansion(mult = c(0, 0)), 
                       breaks = seq(0, plot_days, by = 2))
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(plot.title   = element_text(size = 32, hjust = 0.5),
          axis.title   = element_text(size = 24),
          axis.text.y  = element_text(size = 14),
          axis.text.x  = element_text(size = 14),
          legend.text  = element_text(size = 15),
          legend.title = element_blank(),
          legend.key   = element_blank(),
          legend.position = "bottom", 
          legend.key.height = unit(2, "lines"),
          legend.box.background = element_rect())
  
  # Do not save figure if input is trivial  
  save_fig(o, g, "Viral load profile")
}

# ---------------------------------------------------------
# Plot duration distributions 
# ---------------------------------------------------------
plot_durations = function(o, p) {
  
  # Number of samples from each distribution
  n_samples = 10000
  
  # Colour map and palette to use
  colours = "viridis::viridis" 
  
  # Facet names can be long, wrap after n characters
  n_wrap = 30
  
  # Set a dictionary for facet names and order
  duration_dict = c(latency           = "Duration of latency phase",
                    presymptomatic    = "Duration of presymptomatic phase",
                    non_severe        = "Duration of non-severe disease",
                    severe            = "Duration of severe disease",
                    onset_to_hospital = "Delay between symptom onset and hospitalisation",
                    hospital_stay     = "Days in hospital for non-critical cases", 
                    hospital_to_icu   = "Days in hospital before ICU transfer", 
                    icu_stay          = "Days in intensive care for critical cases", 
                    icu_stay_death    = "Days in intensive care before death", 
                    hospital_transfer = "Days in hospital after ICU discharge",
                    home_death        = "Delay between symptom onset and death outside of hospital")
  
  # ---- Sample distibutions ----
  
  # Preallocate list
  plot_list = list()
  
  # Loop through distributions to sample
  for (this_duration in names(duration_dict)) {
    duration_name = paste0("duration_", this_duration)
    
    # Sample a ton of points 
    duration_samples = p[[duration_name]](rep(0, n_samples))
    
    # Store as a datatable in a list
    plot_list[[this_duration]] = data.frame(distribution = this_duration, 
                                            value = duration_samples)
  }
  
  # Convert to datatable and set factors to preserve ordering
  plot_df = rbindlist(plot_list) %>%
    mutate(distribution = factor(duration_dict[distribution], 
                                 levels = duration_dict))
  
  # ---- Produce plot ----
  
  # Density plots of parameter sample distributions
  g = ggplot(plot_df, aes(x = value, fill = distribution)) + 
    geom_density(adjust = 5, show.legend = FALSE) + 
    facet_wrap(~distribution, scales = "free", 
               labeller = as_labeller(label_wrap_gen(n_wrap)))
  
  # Use pre-defined colour scheme for distribution fills
  g = g + scale_fill_manual(values = rev(colour_scheme(colours, n = length(duration_dict))))
  
  # Set a figure title and axes titles
  g = g + ggtitle("Infection, disease, and hospitalisation duration distributions") + 
    ylab("Probability density") + xlab("Duration, delay, or days until event")
  
  # Set x axes to user-defined limits (see options.R)
  g = g + scale_x_continuous(expand = expansion(mult = c(0, 0))) + 
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(plot.title = element_text(size = 30, hjust = 0.5),
          strip.text = element_text(size = 15),
          axis.title = element_text(size = 24), 
          axis.text  = element_text(size = 12))
  
  # Save to file
  save_fig(o, g, "Duration distributions")
}

# ---------------------------------------------------------
# Plot range of potential seasonality effects
# ---------------------------------------------------------
plot_seasonality = function(o, p, d) {
  
  # Descriptive dictioary for facet strips
  facet_dict = c(temperature = "Temperature (degrees celsius)", 
                 seasonality = "Associated seasonality effect")
  
  # Define temperature and seasonality colours
  colours = c("gold1", "forestgreen")
  
  # Colour for seasoanlity effect bounds
  bound_colour = "grey30"
  
  # Colour for underlying cantonal data informing national level
  canton_colour = "grey70"
  
  # ---- Extract and format data ----
  
  # Scale seasonality effect based on inverse temprerature and seasonality_scaler
  seasonality_fit = 1 + p$seasonality_scaler * (p$temperature_inverse - 1)
  
  # Extract bounds for seasonality factor from calibration parameter table
  seasonality_bounds = filter(o$calibration_df, names == "seasonality_scaler")
  
  # Apply the two bounds to derive the extremes for seasonality effect
  lower = 1 + seasonality_bounds$lower_bound * (p$temperature_inverse - 1)
  upper = 1 + seasonality_bounds$upper_bound * (p$temperature_inverse - 1)
  
  # Construct seasonality dataframe using these vectors - we'll split between lines and area
  seasonality_df = data.table(date = o$dates_all, 
                              lower = lower, upper = upper, 
                              seasonality = seasonality_fit)
  
  # Construct temperature dataframe - including extrapolation
  temp_df = data.table(date = o$dates_all, type = "temperature", 
                       value = p$temperature_extrap)
  
  # Best fit lines for temperature and seasonality
  lines_df = seasonality_df %>%
    select(date, seasonality) %>%
    pivot_longer(cols = seasonality,
                 names_to = "type") %>%
    rbind(temp_df) %>%
    mutate(type = factor(type, levels = names(facet_dict)))
  
  # Construct blank dataframe for colour consistency
  blank_df = data.table(date = o$dates_all[1], 
                        lower = NA, upper = NA, 
                        type = "temperature")
  
  # Constract dataframe for seasonality affect bounds
  area_df = seasonality_df %>%
    select(date, lower, upper) %>%
    mutate(type = "seasonality") %>%
    rbind(blank_df) %>%
    mutate(type = factor(type, levels = names(facet_dict)))
  
  # ---- Produce plot ----
  
  # First plot full range of possible seasonality effects
  g = ggplot(lines_df, aes(x = date)) + 
    geom_ribbon(data = area_df, aes(ymin = lower, ymax = upper, fill = type), 
                alpha = 0.5, show.legend = FALSE) + 
    geom_line(data = area_df, aes(y = lower), colour = bound_colour, size = 0.5) + 
    geom_line(data = area_df, aes(y = upper), colour = bound_colour, size = 0.5) + 
    facet_wrap(~type, scales = "free", labeller = as_labeller(facet_dict), nrow = 2)
  
  # If plotting national level we'll also show how temperature is calculated
  if ("national_weather" %in% names(d)) {
    
    # Format cantonal level temperatures
    canton_df = d$national_weather %>%
      rename(value = temperature) %>%
      mutate(type = "temperature", 
             type = factor(type, levels = names(facet_dict)))
    
    # Plot temperature of each canton we have data for
    g = g + geom_line(data = canton_df, aes(y = value, group = canton), 
                      colour = canton_colour, size = 1, alpha = 0.5)
  }
  
  # Finally, plot overall temperature and best estimate seasonality effect
  g = g + geom_line(data = canton_df, aes(y = value, group = canton), 
                    colour = canton_colour, size = 1, alpha = 0.5) + 
    geom_line(aes(y = value, colour = type), size = 2, show.legend = FALSE) + 
    facet_wrap(~type, scales = "free", labeller = as_labeller(facet_dict), nrow = 2)
  
  # Use pre-defined colour scheme for distribution fills
  g = g + scale_fill_manual(values = colours) + 
    scale_colour_manual(values = colours)
  
  # Add a vertical reference line where data end and prediction starts
  g = g + geom_vline(xintercept = max(o$dates_data), linetype = "dashed", colour = o$dash_colour)
  
  # Set x axes to user-defined limits (see options.R)
  g = g + expand_limits(y = c(0, 1)) + 
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) + 
    scale_x_date(date_breaks = o$x_tick_dates, date_labels = o$x_tick_format)
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(strip.text  = element_text(size = 22),
          axis.title  = element_blank(),
          axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14, hjust = 1, angle = 50))
  
  # Save figure to file
  save_fig(o, g, "Seasonality profile")
}

# ---------------------------------------------------------
# Plot vaccine efficacy profiles 
# ---------------------------------------------------------
plot_vaccine_efficacy = function(o, p) {
  
  # Which vaccines to plot, along with a description
  vaccine_dict = c(pfizer = "Pfizer-BioNTech / Moderna mRNA vaccine", 
                   oxford = "Oxford-AstraZeneca viral vector vaccine")
  
  # Colour scheme
  colours = c("navy", "goldenrod2")
  
  # Number of days to plot (from first dose)
  n_days = 35
  
  # ---- Format data frames ----
  
  # Convert vaccine profiles to percentages
  plot_df = p$vaccine_profile %>% 
    mutate(day = 1 : o$n_dates) %>%
    select(day, one_of(names(vaccine_dict))) %>%
    pivot_longer(cols = -day, names_to = "vaccine") %>%
    arrange(vaccine, day) %>%
    filter(day <= n_days) %>%
    mutate(value = value * 100, 
           vaccine = factor(vaccine_dict[vaccine], 
                            levels = unname(vaccine_dict)))
  
  # Maximum efficacy of each vaccine - we'll plot a reference line
  max_efficacy = unlist(p$vaccine_efficacy[names(vaccine_dict)])
  
  # Construct a plotting dataframe for these reference lines
  max_df = data.frame(vaccine = names(max_efficacy), 
                      max = unname(max_efficacy) * 100) %>%
    mutate(vaccine = factor(vaccine_dict[vaccine], 
                            levels = unname(vaccine_dict)))
  
  # ---- Produce plot ----
  
  # Plot the profiles over time along with maximum efficacies
  g = ggplot(plot_df, aes(x = day, y = value, colour = vaccine)) + 
    geom_hline(data = max_df, aes(yintercept = max, colour = vaccine), 
               linetype = "dotdash", size = 1, alpha = 0.8) + 
    geom_line(size = 3)
  
  # Use pre-defined colour scheme for distribution fills
  g = g + scale_colour_manual(values = colours)
  
  # Set axis labels
  g = g + xlab("Days since first dose") + 
    ylab("Vaccine efficacy (%)")
  
  # Set x axes to user-defined limits (see options.R)
  g = g + expand_limits(y = c(0, 100)) + 
    scale_y_continuous(expand = expansion(mult = c(0, 0)), 
                       breaks = seq(0, 100, by = 10)) + 
    scale_x_continuous(expand = expansion(mult = c(0, 0)), 
                       breaks = seq(0, n_days, by = 5))
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(axis.title   = element_text(size = 26),
          axis.text    = element_text(size = 14),
          legend.text  = element_text(size = 18),
          legend.title = element_blank(),
          legend.key   = element_blank(),
          legend.position = "bottom", 
          legend.key.height = unit(2, "lines"),
          legend.box.background = element_rect())
  
  # Save figure to file
  save_fig(o, g, "Vaccine efficacy profile")
}

# ---------------------------------------------------------
# Plot future testing assumptions
# ---------------------------------------------------------
plot_future_testing = function(o, p) {
  
  # Select past and future colours
  colour_past   = "grey40"
  colour_future = "slateblue4"
  
  # Dates in the past and in the future
  dates_past   = o$dates_all[o$dates_all <= max(o$dates_data)]
  dates_future = o$dates_all[o$dates_all > max(o$dates_data)]
  
  # Diagnosis ratio in the past and in the future
  ratio_past   = p$diagnosis_ratio[o$dates_all %in% dates_past]
  ratio_future = p$future_ratio[o$dates_all %in% dates_future]

  # Construct past and future dataframes
  past_df   = data.frame(date = dates_past,   value = ratio_past)
  future_df = data.frame(date = dates_future, value = ratio_future)
  
  # Plot past as points and future as a line
  g = ggplot(past_df, aes(x = date, y = value)) +
    geom_point(colour = colour_past, size = 3) +
    geom_point(data = future_df, colour = colour_future, size = 4)
  
  # Add a vertical reference line where data end and prediction starts
  g = g + geom_vline(xintercept = max(o$dates_data), linetype = "dashed", colour = o$dash_colour)
  
  # Set axis labels
  g = g + ggtitle("Number of model-calculated diagnoses per infected case")
  
  # Set x axes to user-defined limits (see options.R)
  g = g + expand_limits(y = c(0, NA)) + 
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
    scale_x_date(date_breaks = o$x_tick_dates, date_labels = o$x_tick_format)
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(plot.title   = element_text(size = 32, hjust = 0.5),
          axis.title   = element_blank(),
          axis.text.y  = element_text(size = 14),
          axis.text.x  = element_text(size = 14, hjust = 1, angle = 50))
  
  # Save figure to file
  save_fig(o, g, "Future testing assumptions")
}

# ---------------------------------------------------------
# Plot performance of Gaussian process emulators
# ---------------------------------------------------------
plot_emulator = function(o, sample_num = NULL, by_resampling = FALSE) {
  
  # Zoom into top percentiles of data to see what's going on (values between 1 and 99)
  #
  # NOTE: This works best when you define 3 different levels - 100% is shown as standard
  zoom_in = c(50, 25, 10) # Set to null to turn off
  
  # Remove actual or predicted outliers (3 sd) from plot
  remove_outliers = FALSE
  
  # Should zero by plotted to show distance from perfect fit? 
  show_zero = TRUE
  
  # Define figure title
  fig_title = "Emulator performance"
  
  # ---- Load emulator ----
  
  # Load the full emulator
  if (is.null(sample_num) || by_resampling == TRUE) {
    emulator_name = "model_emulator"
    
  } else {  # ... or the emulator from a given sampling iteration
    emulator_name = paste0("emulator", sample_num)
    
    # Explain as much in the figure title
    fig_title = paste0(fig_title, " (sample ", sample_num, ")")
  }
  
  # Load emulator of choice - throw an error if it doesn't exist explaining step 1 needs to run first
  err_msg = "Attempting to plot emulator performance but cannot find model emulator file"
  emulator = try_load(o$pth$emulator, emulator_name, msg = err_msg) 
  
  # ---- Construct plotting dataframe ----
  
  # Classic use case - construct train and test dataframes then concatenate
  if (by_resampling == FALSE) {
    
    # Construct dataframe of actual likelihood vs predicted: train data
    train_df = data.frame(actual  = emulator$actual_train, 
                          predict = emulator$predict_train, 
                          group   = "Train")
    
    # Construct dataframe of actual likelihood vs predicted: test data
    test_df  = data.frame(actual  = emulator$actual_test,  
                          predict = emulator$predict_test,  
                          group   = "Test")
    
    # Bind together and remove outliers (> 3 std devs from mean)
    samples_df = rbind(train_df, test_df)
    
    # A bit more to do if interested in resampling performance
  } else if (by_resampling == TRUE) {
    
    # Explain as much in the figure title
    fig_title = paste0(fig_title, " (by sampling iteration)")
    
    # Load emulator of choice - throw an error if it doesn't exist explaining step 1 needs to run first
    err_msg = "Attempting to plot emulator performance but cannot find all_samples file"
    param_sets = try_load(o$pth$emulator, "all_samples", msg = err_msg)
    
    # Remove any that failed
    samples_mat = param_sets$samples %>%
      select(-sim_id, -sample_id, -param_id, -seed, -likelihood) %>%
      as.matrix()
    
    # Predict likelihood at all these points using the emulator
    predictions = stats::predict(x = samples_mat, object = emulator$model)$mean
    
    # Construct plotting dataframe which includes (re)sampling numbers
    samples_df = param_sets$samples %>%
      mutate(group   = factor(sample_id, levels = 0 : o$n_init_samples), 
             predict = predictions) %>%
      filter(!is.na(likelihood)) %>%
      select(group, actual = likelihood, predict)
  }
  
  # Zoom in on the best simulations so we can so what's going on
  zoom_frames = c(1, rev(sort(zoom_in)) / 100)
  
  # Threshold values for the defined quantiles
  zoom_vals  = quantile(samples_df$actual, zoom_frames)
  zoom_names = paste0("Top ", names(zoom_vals), " of simulations")
  
  # Loop through all defined quantiles (may just be 100%)
  zoom_list = list()
  for (i in 1 : length(zoom_frames)) {
    
    # Sort the simulations be simulated likelihood
    zoom_list[[i]] = samples_df %>%
      filter(actual <= zoom_vals[i]) %>%
      mutate(zoom_frame = zoom_names[i])
  }
  
  # Convert to factors to retain ordering
  plot_df = rbindlist(zoom_list) %>%
    mutate(zoom_frame = factor(zoom_frame, levels = zoom_names))
  
  # Construct dummy dataframe so we get 
  dummy_df = select(plot_df, actual = predict, predict = actual, group, zoom_frame)
  
  # Remove outliers if desired (> 3 std devs from mean)
  if (remove_outliers == TRUE)
    plot_df = filter(plot_df, 
                     abs(actual  - mean(actual))  < 3 * sd(actual), 
                     abs(predict - mean(predict)) < 3 * sd(predict))
  
  # ---- Produce plot ----
  
  # Plot truth vs predicted (also plot invisible in reverse for square axes)
  g = ggplot(plot_df, aes(x = actual, y = predict, colour = group)) + 
    geom_point() + geom_abline() + geom_blank(data = dummy_df)
  
  if (length(unique(plot_df$zoom_frame)) > 1)
    g = g + facet_wrap(~zoom_frame, scales = "free")
  
  # Overwrite upper limit to zero if desired
  if (show_zero) {
    
    # Also plot dahsed constants through the origin
    g = g + geom_vline(xintercept = 0, linetype = "dashed") + 
      geom_hline(yintercept = 0, linetype = "dashed")
  }
  
  # Prettify axes
  g = g + labs(x = "Actual", y = "Predicted", title = fig_title) + 
    scale_x_continuous(labels = comma) + scale_y_continuous(labels = comma)
  
  # Prettify legend - depends on type of plot
  g = g + scale_colour_discrete(name = ifelse(by_resampling, "Sampling iteration", "Data type"))
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(plot.title   = element_text(size = 24, hjust = 0.5), 
          strip.text.x = element_text(size = 16), 
          strip.text.y = element_text(size = 12), 
          axis.title   = element_text(size = 20), 
          axis.text.x  = element_text(size = 10, hjust = 1, angle = 50), 
          axis.text.y  = element_text(size = 10), 
          legend.text  = element_text(size = 12), 
          legend.title = element_text(size = 14), 
          axis.line    = element_blank(), 
          panel.border = element_rect(colour = "black", fill = NA, size = 1))
  
  # Save figure to file
  save_fig(o, g, fig_title)
}

# ---------------------------------------------------------
# Plot parameter-likelihood relationship
# ---------------------------------------------------------
plot_param_likelihood = function(o) {
  
  # TODO: Prettify these figures...
  
  # Zoom in to focus on top percentage of parameter sets
  zoom_in = 10 # Set to 100 to turn off
  
  # ---- Extract samples ----
  
  # Load emulator of choice - throw an error if it doesn't exist explaining step 1 needs to run first
  err_msg = "Attempting to plot emulator performance but cannot find all_samples file"
  param_sets = try_load(o$pth$emulator, "all_samples", msg = err_msg)
  
  # Re-transform all parameters back to real scale
  param_set_df = param_sets$samples %>% 
    select(-sim_id, -sample_id, -param_id, -seed, -likelihood)
  param_values = t(t(as.matrix(param_set_df)) * param_sets$scale + param_sets$center)
  
  # Extract the top x% best fitting parameter set id
  best_param_id = param_sets$samples %>%
    group_by(param_id) %>%
    summarise(mean_likelihood = mean(likelihood, na.rm = TRUE)) %>%
    slice_min(mean_likelihood, prop = zoom_in / 100)
  
  # Extract parameter set associated with the very best of these
  best_sample = param_sets$samples %>%
    filter(param_id == best_param_id$param_id[1]) %>%
    slice_min(likelihood , n = 1) %>% 
    select(-sim_id, -sample_id, -param_id, -seed, -likelihood)
  
  # Extract parameters and transform back to real scale
  best_simulated = unlist(best_sample) * param_sets$scale + param_sets$center
  
  # ---- Produce plot ----
  
  # Format into plotting dataframe
  plot_df1 = param_sets$samples %>% 
    select(sample = sample_id, param_id, seed, likelihood) %>% 
    mutate(sample = as.factor(sample)) %>% 
    cbind(param_values) %>%
    filter(param_id %in% best_param_id$param_id) %>% 
    pivot_longer(cols = o$calibration_df$names,
                 names_to = "param") %>% 
    rename(x = value) %>% 
    as.data.frame()
  
  # We'll also plot the best fitting parameter set here
  vline_df1 = data.frame(param = names(best_simulated), 
                         x = as.numeric(best_simulated))
  
  # Plot each parameter into a facet
  g1a = ggplot(arrange(plot_df1, likelihood), 
               aes(x = x, y = likelihood)) +
    geom_point(colour = "grey", alpha = 0.2) +
    geom_smooth(formula = y ~ x, method = "loess") + 
    geom_vline(data = vline_df1, aes(xintercept = x), 
               col = "darkblue", size = 2) + 
    facet_wrap(~param, scales = "free_x")
  
  # Plot each parameter into a facet and differentiate by sampling iteration
  g1b = ggplot(arrange(plot_df1, sample), 
               aes(x = x, y = likelihood)) + 
    geom_point(aes(colour = sample), alpha = 0.5) + 
    geom_vline(data = vline_df1, aes(xintercept = x), 
               col = "darkblue", size = 2) + 
    facet_wrap(~param, scales = "free_x")
  
  # Also plot known best parameters if running calibration test
  if (!is.null(o$known_params)) {
    
    # Construct dataframe of known values for each parameter
    vline_df2 = data.frame(param = names(o$known_params), x = unname(o$known_params))
    
    # Add vertical lines to plot
    g1a = g1a + geom_vline(data = vline_df2, aes(xintercept = x), 
                           col = "darkgreen", linetype = "dashed", size = 2)
    g1b = g1b + geom_vline(data = vline_df2, aes(xintercept = x), 
                           col = "darkgreen", linetype = "dashed", size = 2)
  }
  
  # Save figures
  save_fig(o, g1a, paste0("Parameters vs likelihood (top ", zoom_in, "prct)"))
  save_fig(o, g1b, "Parameters vs likelihood (by sampling iteration)")
  
  # ---- Figure B: One-dimension emulated space ----
  
  # Load emulator
  emulator = readRDS(paste0(o$pth$emulator, "model_emulator.rds"))
  
  # Only produce this plot if we have a single parameter
  if (ncol(param_values) == 1) {
    
    # Default plotting settings
    param_idx = 1
    n_samples = 1000
    
    # Take prior means of all parameters, transform and repeat into matrix
    x = (o$calibration_df$values - param_sets$center) / param_sets$scale
    x = t(matrix(x, nrow = nrow(o$calibration_df), ncol = n_samples))
    
    # Uniformly sample parameter of interest between the bounds on the real scale
    x_real = seq(o$calibration_df$lower_bound[param_idx], 
                 o$calibration_df$upper_bound[param_idx], 
                 length.out = n_samples)
    
    # Transform to standardised scale in order to be evaluated and overwrite into matrix
    x[, param_idx] = (x_real - param_sets$center[param_idx]) / param_sets$scale[param_idx]
    
    # Evaluate all points with the model emulator
    z = stats::predict(x = x, object = emulator$model)
    
    # Construct a dataframe including emulator uncertainty bounds on the standardised scale
    plot_df2 = data.frame(x = x_real, likelihood = z$mean, 
                          y_min = qnorm(0.05, z$mean, sqrt(z$sd2 + z$nugs)), 
                          y_max = qnorm(0.95, z$mean, sqrt(z$sd2 + z$nugs)))
    
    # Plot emulator-evaluated points with bounds and overlay model-evaluated training data
    g2 = ggplot(plot_df2, aes(x = x, y = likelihood)) + 
      geom_ribbon(aes(ymin = y_min, ymax = y_max), alpha = 0.4) +
      geom_hline(aes(yintercept = 0), col = "black", linetype = "dashed")
    
    g2a = g2 + geom_point(data = plot_df1, mapping = aes(colour = seed)) + geom_line(size = 2)
    g2b = g2 + geom_point(data = plot_df1, mapping = aes(colour = sample)) + geom_line(size = 2)
    
    if (!is.null(emulator$data$known_params)) {
      g2a = g2a + geom_vline(data = vline_df2, aes(xintercept = x), col = "darkgreen", size = 2)
      g2b = g2b + geom_vline(data = vline_df2, aes(xintercept = x), col = "darkgreen", size = 2)
    }
    
    # Save figure to file
    save_fig(o, g2a, canton, "Emulator vs likelihood (by seed)")
    save_fig(o, g2b, canton, "Emulator vs likelihood (by sampling iteration)")
  }
}

# ---------------------------------------------------------
# Plot grouped calibration posteriors vs priors
# ---------------------------------------------------------
plot_posteriors = function(o, n_bins = 20) {
  
  # Initiate dataframes
  mcmc_df = dens_df = NULL
  
  # Table of calibrated parameters
  calibration_table = o$parameter_table %>% 
    filter(fixed == 0) # %>% select(-description)
  
  # ---- Extract canton chains ----
  
  # Loop through cantons to load chains
  for (canton in o$cantons) {
    
    # Load canton fitting file
    fit = readRDS(paste0(o$pth$fitting, canton, ".rds"))
    
    # Remove reference to canton name
    param_names  = str_replace(names(fit$chain),        paste0(".", canton), "")
    known_params = str_replace(names(fit$known_params), paste0(".", canton), "")
    
    # Re-transform all parameters back to real scale
    param_values = fit$chain %>% select(-sampno, -lnlike)
    
    # Retransform parameters via sample_posteriors function
    canton_df = as.data.frame(param_values)
    names(canton_df) = param_names[param_names %in% calibration_table$names]
    
    # Bind all canton chains together
    canton_df = cbind(canton, canton_df)
    mcmc_df   = bind_rows(mcmc_df, canton_df)
  }
  
  # ---- Sample from distributions ----
  
  # Names of all parameters fitted in this analysis
  param_names = intersect(names(mcmc_df), calibration_table$names)
  
  # Loop through fitted parameters
  for (param_name in param_names) {
    param = calibration_table[calibration_table$names == param_name, ]
    
    # Randomly sample a load of points from prior distribution
    prior_df = data.frame(param_name = param_name, 
                          param_type = "Prior",
                          samples = rnorm(1e6, param$values, param$prior_sd)) %>%
      filter(samples >= param$lower_bound, 
             samples <= param$upper_bound)
    
    # Bind all posterior samples into density dataframe
    dens_df = rbind(dens_df, prior_df)
    
    # ---- Posterior ----
    
    # Each canton may have it's own posterior
    for (canton in o$cantons) {
      
      # Posterior samples come from MCMC chains
      posterior_df = data.frame(param_name = param_name, 
                                param_type = ifelse(param$global == 1, "Global", canton), 
                                samples    = mcmc_df[mcmc_df$canton == canton, param_name])
      
      # Bind all posterior samples into density dataframe
      dens_df = rbind(dens_df, posterior_df)
      
      # TODO: Escape loop here?
      # if (param_type == "Global") break
    }
  }
  
  # Convert types to factors to preserve ordering
  dens_df$param_type = factor(dens_df$param_type, levels = c("Prior", o$cantons, "Global"))
  
  # ---- Produce plot for each canton ----
  
  # Concatenate colour scheme
  colours = c(o$prior_colour, o$posterior_colour, o$global_colour)
  
  # Also plot for each canton
  for (canton in o$cantons) {
    
    # Convert parameter names to factors to preserve ordering
    dens_df$param_name = factor(dens_df$param_name, levels = calibration_table$names)
    
    # Density plots of parameter sample distributions
    g = ggplot(dens_df, aes(x = samples)) + 
      geom_density(aes(y = ..scaled.., fill = param_type), 
                   show.legend = FALSE, alpha = 0.5, colour = NA) +
      facet_wrap(~param_name, scales = "free_x") +
      scale_fill_manual(values = colours)
    
    # Also plot known best parameters if running calibration test
    if (length(known_params) > 0) {
      
      # Construct dataframe of known values for each parameter
      vline_df = data.frame(param_name = known_params, x = unname(fit$known_params))
      
      # Add vertical lines to plot
      g = g + geom_vline(data = vline_df, aes(xintercept = x), col = "darkgreen", size = 2)
    }
    
    # Prettify axes
    g = g + scale_y_continuous(expand = c(0, NA)) + 
      xlab("Parameter value") + ylab("Probability density")
    
    # Prettify theme
    g = g + theme_classic() + 
      theme(strip.text = element_text(size = 12),
            axis.title = element_text(size = 18),
            axis.text.x  = element_text(size = 10),
            axis.text.y  = element_blank(),
            axis.ticks.y = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA, size = 1))
    
    # Save figure
    save_fig(o, g, "Posteriors vs priors", canton)
  }
}

# ---------------------------------------------------------
# Plot multivariate MCMC traces
# ---------------------------------------------------------
plot_posterior_traces = function(o) {
  
  # TODO: Would it make sense to set y limits to parameter bounds?
  
  # Load all MCMC chains
  all_chains = readRDS(paste0(o$pth$fitting, "all_chains.rds"))
  
  # Melt MCMC chains for plotting
  plot_df = all_chains %>% 
    pivot_longer(cols = c(-sampno, -chain_id), 
                 names_to = "param") %>% 
    arrange(param, chain_id, sampno) %>% 
    as.data.frame()
  
  # We'll join the scope of each parameter
  scope_df = o$calibration_df %>%
    select(names, scope) %>%
    rename(param = names) %>%
    as.data.frame() %>%
    rbind(data.frame(param = "lnlike", scope = "lnlike"))
  
  # Do the joining
  plot_df = left_join(plot_df, scope_df, by = "param")
  
  # We want to plot each scope (aside from lnlike)
  all_scopes = setdiff(unique(plot_df$scope), "lnlike")
  
  # Scope is unnecessary if plotting just a single canton
  if (length(o$cantons) == 1)
    plot_df$scope = all_scopes = o$cantons
  
  # Loop through the scopes (inc. global and hyper parameters)
  for (this_scope in all_scopes) {
    plot_scope_df = plot_df[plot_df$scope %in% c(this_scope, "lnlike"), ]
    
    # Rename parameters - doing this for pretty publication plot
    #
    # WARNING: Does this destroy generalisability when plotting multiple cantons?
    plot_scope_df$name = str_remove(plot_scope_df$param, "\\.[A-Z]*")
    
    # Plot all traces from multivariate chain
    g = ggplot(plot_scope_df) +
      geom_line(aes(x = sampno, y = value, col = as.factor(chain_id))) +
      facet_wrap(~name, scales = "free_y", ncol = 4)
    
    # Prettify colour scheme
    g = g + scale_color_viridis_d()
    
    # Prettify axes
    g = g + scale_x_continuous(expand = c(0, 0)) + 
      xlab("MCMC iteration") + ylab("Parameter value")
    
    # Prettify theme
    g = g + theme_classic() + 
      theme(legend.position = "none", 
            strip.text = element_text(size = 12),
            axis.title = element_text(size = 18),
            axis.text  = element_text(size = 10),
            panel.border = element_rect(colour = "black", fill = NA, size = 1))
    
    # Save new figure for each scope
    save_fig(o, g, "Traceplots", first_cap(this_scope))
  }
}

# ---------------------------------------------------------
# Plot epi data to examine different sources
# ---------------------------------------------------------
plot_epi_data = function(o, foph_data, openzh_data) {
  
  # Melt down FOPH data and add source
  plot_foph = foph_data %>% 
    filter(canton %in% o$cantons) %>% 
    pivot_longer(cols = -c("date", "canton"), 
                 names_to = "metric") %>% 
    mutate(source = "FOPH")
  
  # Same for OpenZH data
  plot_openzh = openzh_data %>% 
    filter(canton %in% o$cantons) %>% 
    pivot_longer(cols = -c("date", "canton"), 
                 names_to = "metric") %>% 
    mutate(source = "OpenZH")
  
  # Full list of metrics (full names) we are about to plot - taken from both data sources
  plot_metrics = as.character(o$metric_dict[unique(c(plot_foph$metric, plot_openzh$metric))])
  
  # Bind by rows, which are now consistent
  plot_df = rbind(plot_foph, plot_openzh) %>%
    filter(!is.na(value)) %>%
    mutate(metric = o$metric_dict[metric],  # Only way I can seem to get label_wrap_gen working
           metric = factor(metric, levels = plot_metrics)) %>%
    as.data.frame()
  
  # Plot both data sources in a canton x metric grid
  g = ggplot(plot_df, aes(x = date, y = value, colour = source)) + 
    geom_line(size = 1.5) + 
    geom_vline(xintercept = max(o$dates_data), linetype = "dashed", colour = "black") +
    geom_vline(xintercept = format_date(today()), linetype = "dotdash", colour = "darkgreen") +
    facet_grid(metric ~ canton, scales = "free_y", 
               labeller = as_labeller(label_wrap_gen(15)))
  
  # Set x axes to user-defined limits (see options.R)
  g = g + scale_y_continuous(labels = comma, expand = c(0, 0), limits = c(0, NA)) +
    scale_x_date(limits = c(min(o$dates_data), today() + 7),
                 date_breaks = o$x_tick_dates, date_labels = o$x_tick_format)
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(strip.text = element_text(size = 14),
          axis.title = element_blank(),
          axis.text.y  = element_text(size = 12),
          axis.text.x  = element_text(size = 12, hjust = 1, angle = 50),
          legend.text  = element_text(size = 14),
          legend.title = element_blank(),
          legend.box.background = element_rect(),
          panel.border = element_rect(colour = "black", fill = NA, size = 1))
  
  # Save figure to file
  save_fig(o, g, "Data", "Epi sources")
}

# ---------------------------------------------------------
# Plot testing data and smoothed model inputs
# ---------------------------------------------------------
plot_smoothed_data = function(o, smoothed_df, variable, fig_name) {
  
  # TODO: Could be prettier...
  
  # Format data into seprate dataframe
  data_df = smoothed_df %>%
    select(date, canton, data = one_of(variable)) %>%
    filter(!is.na(data))
  
  # Subset for only extrapolations and melt down
  plot_df = smoothed_df %>% 
    select(date, canton, sma_fitted, sma_centre) %>%
    pivot_longer(cols = -c("date", "canton"), 
                 names_to = "method") %>% 
    filter(!is.na(value), 
           date >= min(o$dates_data), 
           date <= max(o$dates_data)) %>%
    as.data.frame()
  
  # Plot the curves with the data on top
  g = ggplot(plot_df, aes(x = date, y = value, colour = method)) + 
    geom_line(size = 2) +
    geom_point(data = data_df, mapping = aes(y = data), colour = "black") + 
    facet_wrap(~canton, scales = "free_y")
  
  # Set x axes to user-defined limits (see options.R)
  g = g + scale_y_continuous(labels = comma, expand = c(0, 0)) +
    scale_x_date(limits = c(min(o$dates_data), max(o$dates_data)),
                 date_breaks = o$x_tick_dates, date_labels = o$x_tick_format)
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(strip.text = element_text(size = 18),
          axis.title = element_blank(),
          axis.text.y  = element_text(size = 16),
          axis.text.x  = element_text(size = 16, hjust = 1, angle = 50),
          legend.text  = element_text(size = 16), 
          legend.title = element_blank(),
          legend.box.background = element_rect(),
          panel.border = element_rect(colour = "black", fill = NA, size = 1))
  
  # Save figure to file
  save_fig(o, g, "Data smoothing", fig_name)
}

# ---------------------------------------------------------
# Plot testing data & assumptions - mild to severe case ratio
# ---------------------------------------------------------
plot_diagnosis_ratio = function(o, epi_data) {
  
  # Extract number diagnosed and number hospitalised from data
  plot_df = epi_data %>%
    filter(metric %in% c("confirmed", "hospital_admissions")) %>%
    pivot_wider(id_cols = qc(date, canton),  
                names_from  = "metric", 
                values_from = "value") %>%
    right_join(data.frame(date = o$dates_data, canton = o$cantons), 
               by = qc(date, canton)) %>%
    arrange(canton, date) %>%
    rename(diagnosed_hospital = hospital_admissions) %>%
    mutate(diagnosed_other = pmax(confirmed - diagnosed_hospital, 0), 
           case_per_hosp = diagnosed_other / diagnosed_hospital) %>%
    filter(is.finite(case_per_hosp)) %>%
    as.data.table()
  
  # Plot ratio of non-hospitalised cases to each hospitalised case
  g = ggplot(plot_df, aes(x = date, y = case_per_hosp)) + 
    geom_smooth(method = 'lm', formula = y~x) +
    geom_point(aes(colour = confirmed)) + 
    facet_wrap(~canton, ncol = 1, scales = "free_y")
  
  # Set a nice continuous colour scheme for number 
  g = g + scale_colour_gradient(name = "Number of cases", 
                                low = "grey60", high = "firebrick2")
  
  # Set x axes to user-defined limits (see options.R)
  g = g + ggtitle("Non-hospital cases diagnosed per hospitalisation") + 
    scale_y_continuous(labels = comma, expand = c(0, 0)) +
    scale_x_date(limits = c(min(o$dates_data), max(o$dates_data)),
                 date_breaks = o$x_tick_dates, date_labels = o$x_tick_format)
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(plot.title   = element_text(size = 24, hjust = 0.5),
          strip.text   = element_text(size = 18),
          axis.title   = element_blank(),
          axis.text.y  = element_text(size = 16),
          axis.text.x  = element_text(size = 16, hjust = 1, angle = 50),
          legend.text  = element_text(size = 16), 
          legend.title = element_text(size = 18),
          panel.border = element_rect(colour = "black", fill = NA, size = 1))
  
  # Save figure to file
  save_fig(o, g, "Data", "Diagnoses ratio")
}

# ---------------------------------------------------------
# Produce output csv files to accompany figures
# ---------------------------------------------------------
output_csv = function(o, canton, file_name, bounds = TRUE, ...) {
  
  # Collate and interpret function inputs so we know what to output
  f = figure_properties(o, scenario_groupings_ok = TRUE, ...)
  
  # Do not allow relative difference to baseline for csv outputs
  if (f$relative == TRUE)
    stop("Set 'relative' to FALSE for csv outputs - easily calculated from within the file itself")
  
  # Load baseline file for this canton
  a_baseline = readRDS(paste0(o$pth$scenario, canton, "_", f$baseline_name, ".rds"))
  
  # Inherit key settings from analysis file (see store_options() function)
  o[names(a_baseline$options)] = a_baseline$options
  
  # Initiate output datatable
  output_list = list()
  
  # Loop through metrics to plot
  for (scenario in f$scenarios) {
    
    # Load data file for this canton for this scenario
    if (scenario == f$baseline_name) a = a_baseline else {
      a = try_load(o$pth$scenario, paste0(canton, "_", scenario)) 
    }
    
    # Format model output and store in list to be concatenated
    output_list[[scenario]] = format_results(o, f, a)
  }
  
  # Concatenate output datatables for all scenarios
  output_df = rbindlist(output_list) %>%
    mutate(scenario = f$scenario_dict[scenario])
  
  # Remove bounds if requested
  if (bounds == FALSE)
    output_df = select(output_df, -lower, -upper)
  
  # Construct file name and concatenate with file path
  save_name = paste(file_name, collapse = " - ")
  save_pth  = paste0(o$pth$figures, save_name, ".csv")
  
  # Save model outcomes to file
  utils::write.csv(output_df, file = save_pth, row.names = FALSE)
}

# ---------------------------------------------------------
# Collate plotting function inputs so we know what to plot
# ---------------------------------------------------------
figure_properties = function(o, plot_baseline = TRUE, alt_baseline = NULL, scenarios = NULL, strategies = NULL, 
                             plot_metrics = NULL, plot_by = NULL, everything = FALSE, cumulative = FALSE, relative = FALSE, 
                             plot_from = NULL, plot_to = NULL, expand_limit = 1, legend_rows = 2, facet_rows = NULL, 
                             override_colours = NULL, override_fontsize = NULL, scenario_groupings_ok = FALSE) {
  
  # Initiate figure list
  f = list()
  
  # ---- Strategies and scenarios to plot ----
  
  # Initiate vector of all scenarios to run
  all_scenarios = NULL
  
  # Plot baseline (either default baseline or an alternative baseline)
  if (plot_baseline == TRUE) {
    
    # Name of baseline (either default or an alternative)
    if (is.null(alt_baseline)) f$baseline_name = "baseline" else f$baseline_name = alt_baseline
    
    # Append baseline to scenario vector
    all_scenarios = f$baseline_name
  }
  
  # Append all non-baseline scenarios to scenario vector
  all_scenarios = union(all_scenarios, scenarios)
  
  # Initiate scenario dictionary with all scenarios extracted so far
  scenario_dict = scenario_names(o)[all_scenarios[!grepl("\\.", all_scenarios)]]
  
  # Loop through all strategies to plot together (commonly none or one)
  for (strategy in strategies) {
    
    # Load strategy info
    strategy_info = try_load(o$pth$strategy, paste0("strategy.", strategy))
    
    # Scenarios to extract from this strategy and associated names
    these_scenarios = setdiff(strategy_info$combinations$id, all_scenarios)
    these_names = filter(strategy_info$combinations, id %in% these_scenarios)$name
    
    # Append names of all scenarios in strategy
    all_scenarios = c(all_scenarios, these_scenarios)
    
    # Concatenate to dictionary of scenario names and IDs
    scenario_dict = c(scenario_dict, setNames(these_names, these_scenarios))
  }
  
  # Number of scenarios identified
  f$scenarios = all_scenarios
  n_scenarios = length(f$scenarios)
  
  # Do we want to plot difference between scenarios and baseline
  if (relative == TRUE) {
    
    # Sanity check that a baseline is provided
    if (plot_baseline == FALSE)
      stop("You must provide a baseline if 'relative' is set to TRUE")
    
    # Reduce number of scenarios by one as we won't plot baseline explictly
    n_scenarios = n_scenarios - 1
    
    # Sanity check that not only a baseline is provided
    if (n_scenarios == 0)
      stop("To use 'relative' you must provide some alternative scenarios and/or strategies")
  }
  
  # Store flag so outcomes can be relatively calculated if needed
  f$relative = relative
  
  # Throw an error if no scenarios identified
  if (n_scenarios == 0)
    stop("No scenarios identified - plot a baseline and/or scenarios and/or strategies")
  
  # We now need to deal with scenarios that were generated within strategies to extract names
  unnamed_scenarios = setdiff(all_scenarios, names(scenario_dict))
  within_strategies = as.character(na.omit(str_extract(unnamed_scenarios, "^.*(?=(\\.))")))
  
  # Loop through all such strategies 
  for (strategy in within_strategies) {
    
    # Load strategy info
    strategy_info = try_load(o$pth$strategy, paste0("strategy.", strategy))
    
    # Scenarios in this strategy that are thus far unnamed
    these_scenarios = unnamed_scenarios[unnamed_scenarios %in% strategy_info$combinations$id]
    these_names = filter(strategy_info$combinations, id %in% these_scenarios)$name
    
    # Concatenate to dictionary of scenario names and IDs
    scenario_dict = c(scenario_dict, setNames(these_names, these_scenarios))
  }
  
  # If plotting baseline set default name, if not define something from which we can extract data
  if (plot_baseline == TRUE) scenario_dict[f$baseline_name] = o$baseline_name
  else f$baseline_name = all_scenarios[1]
  
  # Reorder dictionary and append
  f$scenario_dict = scenario_dict[all_scenarios]
  
  # Do we have multiple scenarios to plot?
  if (n_scenarios > 1) {
    
    # If so, we can't also have multiple groups on the same plot
    if (!is.null(plot_by)) {
      
      # This constraint can be sidestepped if producing csv files
      if (!scenario_groupings_ok)
        stop("You cannot plot groupings for multiple scenarios at once")
    }
    
    # We'll be plotting by multiple scenarios
    f$plot_type = "scenario"
  }
  
  # ---- Metrics and groupings to plot ----
  
  # Use default if no metric arguments given
  if (is.null(plot_metrics)) f$metrics = o$default_plot_metrics else f$metrics = plot_metrics
  
  # Overwrite with all metrics if everything flag is on
  if (everything == TRUE) f$metrics = names(o$metric_dict)
  
  # Not all metrics can be plotted cumulatively
  if (cumulative == TRUE) {
    cumulative_df = filter(o$all_metrics, cumulative == TRUE)
    
    # Subset metrics to only those suitable for cumulative plotting
    f$metrics = intersect(f$metrics, cumulative_df$metric)
  }
  
  # If no grouping defined we'll be ploting by multiple metrics
  if (is.null(plot_by)) {
    
    # For some metrics it is nonsensical to aggregate
    aggregate_metrics = o$all_metrics$metric[o$all_metrics$aggregate]
    f$metrics = intersect(f$metrics, aggregate_metrics)
    
    # Colours will be by metric
    if (n_scenarios == 1)
      f$plot_type = "metric"
    
  } else {
    
    # We'll be plotting by the grouping defined
    f$plot_type = "group"
    f$plot_by   = plot_by
    
    # Subset metrics that have been requested AND are grouped by plot_by
    grouping_metrics = names(o$metric_groupings[o$metric_groupings == f$plot_by])
    f$metrics = intersect(f$metrics, grouping_metrics)
  }
  
  # Throw an error if there is nothing left
  if (length(f$metrics) == 0)
    stop("No plotting metrics identified")
  
  # ---- Metric descriptions ----
  
  # Metric descriptions can be different if plotting cumulative results
  if (cumulative == TRUE) metric_dict = o$cumulative_dict else metric_dict = o$metric_dict
  
  # Apply and subset the dictionary
  f$metric_names = metric_dict[f$metrics]
  
  # Store flag so outcomes can be cumulative summed if needed
  f$cumulative = cumulative
  
  # ---- Plotting colours ----
  
  # Skip colour generation if override_colours is NA
  if (!is.null(override_colours) && is.na(override_colours)) f$colours = NA
  else {
    
    # Scenario colours: redefine on each function call
    if (f$plot_type == "scenario") {
      
      # Create colour vector from palette, considering baseline if required
      if (plot_baseline == FALSE) f$colours = colour_scheme(o$palette_scenario, n = n_scenarios)
      else f$colours = c(o$baseline_colour, colour_scheme(o$palette_scenario, n = n_scenarios - 1))
    }
    
    # Metric colours: a colour for each possible metric, indexed if not plotting everything
    if (f$plot_type == "metric") {
      
      # Colours for all metrics we may wish to plot
      all_metric_colours = colour_scheme(o$palette_metric, n = length(o$metric_dict))	
      
      # Subset this for metrics we are plotting
      f$colours = all_metric_colours[names(o$metric_dict) %in% f$metrics]
    }
    
    # Colours by group depends on number in grouping
    if (f$plot_type == "group") {
      
      # Number in this grouping
      n_group = length(o[[paste0("count_", f$plot_by)]])
      
      # Grouping by age is a special case
      if (f$plot_by == "age")
        n_group = length(o$plot_ages)
      
      # User defined palette for this grouping
      group_palette = o[[paste0("palette_", f$plot_by)]]
      
      # Generate set of colours
      f$colours = colour_scheme(group_palette, n = n_group)
    }
    
    # Force override colours if desired
    if (!is.null(override_colours)) {
      
      # Throw an error if override does not have correct number of values
      if (length(override_colours) != length(f$colours))
        stop("Inconsistent number of manual colours provided ", 
             "(", length(f$colours), " needed, ", length(override_colours), " provided)")
      
      # Apply the override
      f$colours[!is.na(override_colours)] = na.omit(override_colours)
    }
  }
  
  # ---- Plotting dates ----
  
  # Plotting dates (bounded above and below by o$dates_all)
  f$plot_from = max(min(o$dates_all), format_date(plot_from))
  f$plot_to   = min(max(o$dates_all), format_date(plot_to))
  
  # ---- Plotting aesthetics ----
  
  # Directly apply default or user-defined values
  f$expand_limit = expand_limit
  f$legend_rows  = legend_rows
  f$facet_rows   = facet_rows
  
  # Font size (title, facets, ticks, legend) depends on how many panels we have
  if (length(f$metrics) <= 4)  f$fontsize = c(32, 20, 15, 15) 
  if (length(f$metrics) <= 9)  f$fontsize = c(32, 17, 15, 15) 
  if (length(f$metrics) >= 10) f$fontsize = c(32, 14, 12, 12)
  
  # Force override font sizes if desired
  if (!is.null(override_fontsize)) {
    
    # Throw an error if override does not have correct number of values
    if (length(override_fontsize) != length(f$fontsize))
      stop("Input 'fontsize' must have 4 values (title, facets, ticks, and legend text)")
    
    # Apply the override
    f$fontsize[!is.na(override_fontsize)] = na.omit(override_fontsize)
  }
  
  return(f)
}

# ---------------------------------------------------------
# Save a ggplot figure to file with default settings
# ---------------------------------------------------------
save_fig = function(o, g, ..., path = "figures", width = o$save_width, height = o$save_height) {
  
  # Collapse inputs into vector of strings
  fig_name_parts = unlist(list(...))
  
  # Construct file name and concatenate with file path
  save_name = paste(fig_name_parts, collapse = " - ")
  save_pth  = paste0(o$pth[[path]], save_name, ".png")
  
  # Save figure (size specified in results.R)
  ggsave(save_pth, plot = g, width = width, height = height)
}

