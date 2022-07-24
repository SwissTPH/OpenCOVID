###########################################################
# PLOTTING
#
# All plotting functions in one place.
#
###########################################################

# ---------------------------------------------------------
# Plot metrics over time for multiple metrics, groups, and/or scenarios
# ---------------------------------------------------------
plot_temporal = function(o, fig_name, plot_file = NULL, ...) {
  
  # ---- Load baseline (or alt baseline) file ----
  
  # Scenarios to plot
  f = fig_scenarios(o, list(...))
  
  # Handle alternative functionality case - file fed in directly
  if (!is.null(plot_file)) baseline = plot_file
  else {
    
    # Otherwise load baseline file
    baseline = try_load(o$pth$scenarios, f$baseline_name)
  }
  
  # ---- Figure properties ----
  
  # Collate and interpret function inputs so we know what to plot
  f = fig_properties(o, f, baseline$input, list(...))
  
  # If nothing to plot, return out
  if (f$n_metrics == 0)
    return()
  
  # Do not allow relative difference to baseline for temporal plotting
  if (f$relative == TRUE)
    stop("Set 'relative' to FALSE for temporal plots (use impact bars to assess difference to baseline)")
  
  # ---- Extract model predictions ----
  
  # Initiate plotting dataframe
  plot_list = list()
  
  # Loop through metrics to plot
  for (scenario in f$scenarios) {
    
    # Load data file for this scenario
    if (scenario == f$baseline_name) result = baseline else {
      result = try_load(o$pth$scenarios, scenario) 
    }
    
    # Format model output and store in list to be concatenated
    plot_list[[scenario]] = format_results(o, f, result)  # See postprocess.R
  }
  
  # Normal use case is to reverse scenarios (so baseline is on top)
  if (f$plot_type %in% c("group", "scenario_group")) 
    scenario_levels = f$scenario_names
  else scenario_levels = rev(f$scenario_names)
  
  # Concatenate plotting dataframes for all scenarios
  plot_df = rbindlist(plot_list) %>%
    mutate(metric   = recode(metric,   !!!f$metric_names), 
           scenario = recode(scenario, !!!f$scenario_names), 
           scenario = factor(scenario, levels = scenario_levels))
  
  # Append scenario grouping details if plotting by scenario groups  
  if (f$plot_type == "scenario_group") {
    
    # We need to apply additional factoring to retain correct order
    levels_group = unique(f$scenario_groups$scenario_group)
    levels_type  = unique(f$scenario_groups$scenario_type)
    
    # Append scenario grouping details 
    plot_df = plot_df %>% 
      left_join(f$scenario_groups, by = "scenario") %>%
      mutate(scenario       = factor(scenario,       levels = scenario_levels), 
             scenario_group = factor(scenario_group, levels = levels_group), 
             scenario_type  = factor(scenario_type,  levels = levels_type))
  }
  
  # ---- Apply dates if necessary ----
  
  # Check whether dates are preferred to day numbers
  if (f$plot_dates == TRUE) {
    
    # Vector of all dates simulated
    all_dates = f$date_from + (1 : baseline$input$n_days) - 1 
    
    # Apply these dates
    plot_df$date = all_dates[plot_df$date]
  }
  
  # ---- Create basic figure ----
  
  # In most cases the plot_type difines the colour scheme
  aes_by = list(colour = f$plot_type, fill = f$plot_type)
  
  # We need a slightly different aesthetic configuration if plotting by scenario_group
  if (f$plot_type == "scenario_group") {
    
    # Normal use case is to differentiate types by dashed lines
    aes_by$colour = "scenario"
    aes_by$fill   = "scenario"
    aes_by$line   = "scenario_type"
    
    # However this can also be reversed, so types by colours and groups by dashes
    if (f$aes_reverse == TRUE)
      aes_by = list(colour = "scenario_type", 
                    fill   = "scenario_type", 
                    line   = "scenario_group")
  }
  
  # Turn off linetype is not desired
  if (f$aes_linetype == FALSE)
    aes_by$line = NULL
  
  # Initiate ggplot structure with each metric in it's facet
  g = ggplot(plot_df, aes_string(x = "date", y = "value", 
                                 colour   = aes_by$colour, 
                                 fill     = aes_by$fill, 
                                 linetype = aes_by$line))
  
  # Several options for plotting geom - default is a line graph
  if (f$plot_geom == "line") {
    
    # Plot uncertainty bounds first (if desired)
    if (f$uncertainty == TRUE)
      g = g + geom_ribbon(aes(ymin = lower, ymax = upper), 
                          linetype = 0, alpha = 0.3)
    
    # Plot best estimate line on top
    g = g + geom_line(size = f$line_width)
    
    # Area plot over time - useful when plotting age groups
  } else if (f$plot_geom == "area") {
    g = g + geom_area()
    
  } else {  # Throw an error if any other value is provided
    stop("Value of plot_geom '", f$plot_geom, "' not recognised")
  }
  
  # ---- Generate facets ----
  
  # Possible to define custom faceting function, in the form of evaluable string
  if (!is.null(f$facet_custom))
    eval(parse(text = paste0("g = g + ", f$facet_custom)))
  
  # Most common use case is to use facet_wrap/grid based on what is being varied
  if (is.null(f$facet_custom)) {
    
    # Most common use case is to wrap by metric
    if (f$plot_type %in% c("scenario", "metric"))
      g = g + facet_wrap(~metric, scales = "free", nrow = f$facet_rows, labeller = f$label_wrap)
    
    # Scenario groups use metric x scenario_group grid
    if (f$plot_type == "scenario_group")
      g = g + facet_grid(metric~scenario_group, scales = "free", labeller = f$label_grid)
    
    # If plotting by group, use special faceting (could be busy!)
    if (f$plot_type == "group") {
      
      # Only one metric (and multiple scenarios): wrap by scenarios
      if (f$n_metrics == 1 && f$n_scenarios > 1)
        g = g + facet_wrap(~scenario, nrow = f$facet_rows, labeller = f$label_wrap)
      
      # Only one scenario (n_metrics irrelevant): wrap by metrics
      if (f$n_scenarios == 1)
        g = g + facet_wrap(~metric, scales = "free", nrow = f$facet_rows, labeller = f$label_wrap)
      
      # Multiple metrics and multiple scenarions: grid faceting
      if (f$n_metrics > 1 && f$n_scenarios > 1)
        g = g + facet_grid(metric~scenario, scales = "free", labeller = f$label_grid)
    }
  }
  
  # Tag facets if desired
  if (f$facet_labels == TRUE)
    g = facet_labels(g)
  
  # ---- Effective reproduction number ----
  
  # Add a vertical line at 1 if plotting effective reproduction number
  if ("R_effective" %in% f$metrics) {
    reff_name = f$metric_names[["R_effective"]]
    
    # Plot R_eff = 1, and also the R_eff we've calibrated to
    g = g + 
      geom_hline(data = subset(plot_df, metric == reff_name), aes(yintercept = 1)) + 
      geom_hline(data = subset(plot_df, metric == reff_name), aes(yintercept = baseline$input$r_eff), 
                 linetype = "dashed", colour = "black")
  }
  
  # ---- Figure aesthetics ----
  
  # Set pretty colour scheme (either by metric, group, or scenario)
  if (f$plot_type != "scenario") {
    
    # If plotting by group, extract dictionary (age already considered)
    if (f$plot_type == "group" && f$plot_by != "age") {
      use_labels = baseline$input$dict[[f$plot_by]]
      
    } else {  # Otherwise no additional labels needed
      use_labels = waiver()
    }
    
    # Apply colour scheme (and labels if appropriate)
    g = g + scale_colour_manual(values = f$colours, labels = use_labels) + 
      scale_fill_manual(values = f$colours, labels = use_labels)
  }
  
  # A little more needed when plotting multiple scenarios
  if (f$plot_type == "scenario") {
    
    # Apply scenario dict - done here so more obvious if names are not unique
    g = g + scale_colour_manual(values = rev(f$colours), labels = rev(f$scenario_names)) + 
      scale_fill_manual(values = rev(f$colours), labels = rev(f$scenario_names))
    
    # Reverse legend order so baseline is at the top
    apply_guide = guide_legend(reverse = TRUE, nrow = f$legend_rows, byrow = f$legend_by_row)
    g = g + guides(colour = apply_guide, fill = apply_guide)
  }
  
  # Most legends are irrelevant when plotting by scenario groups
  if (f$plot_type == "scenario_group") {
    apply_guide = guide_legend(nrow = f$legend_rows, byrow = f$legend_by_row)
    
    # Turn off colour legend
    if (!f$aes_reverse)
      g = g + guides(colour = "none", fill = "none", linetype = apply_guide)
    
    # Turn off linetype legend
    if (f$aes_reverse)
      g = g + guides(colour = apply_guide, fill = apply_guide, linetype = "none")  
  }
  
  # Prettify y-axis
  g = g + expand_limits(y = c(0, f$y_min)) +  # Restrain upper y limit from being really low
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, f$y_expand)))
  
  # Prettify x-axis: day numbers
  if (f$plot_dates == FALSE)
    g = g + scale_x_continuous(labels = comma, 
                               limits = c(f$plot_from, f$plot_to), 
                               expand = expansion(mult = c(0, 0)))
  
  # Prettify x-axis: dates
  if (f$plot_dates == TRUE)
    g = g + scale_x_date(date_breaks = f$date_breaks,
                         date_labels = f$date_labels,
                         limits = c(all_dates[f$plot_from], 
                                    all_dates[f$plot_to]), 
                         expand = expansion(mult = c(0, 0)))
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(plot.title    = element_text(size = f$fontsize[1], hjust = 0.5),
          strip.text    = element_text(size = f$fontsize[2], margin = margin(2, 2, 2, 2)),
          axis.text     = element_text(size = f$fontsize[4]),
          axis.title    = element_blank(),
          axis.line     = element_blank(),
          axis.ticks    = element_line(size = 0.25),
          axis.ticks.length = unit(0.1, "lines"),
          panel.border  = element_rect(size = 0.5, colour = "black", fill = NA),
          panel.spacing = unit(0.5, "lines"),
          strip.background = element_blank(), 
          legend.text   = element_text(size = f$fontsize[3]),
          legend.title  = element_blank(),
          legend.key    = element_blank(),
          legend.margin = margin(2, 1, 2, 1),
          legend.position = "bottom", 
          legend.key.height = unit(1, "lines"),
          legend.key.width  = unit(1, "lines"),
          legend.box.background = element_rect(size = 0.5))
  
  # Rotate x tick labels if plotting dates
  if (f$plot_dates == TRUE)
    g = g + theme(axis.text.x = element_text(size = f$fontsize[4], hjust = 1, angle = 50))
  
  # Remove legend entirely if only plotting multiple metrics (the facets are enough)
  if (f$plot_type == "metric")
    g = g + theme(legend.position = "none")
  
  # Do not save figure if input is trivial  
  if (!is.null(fig_name))
    fig_save(o, g, fig_name)
  
  return(g)
}

# ---------------------------------------------------------
# Plot metrics over time for multiple metrics, groups, and/or scenarios
# ---------------------------------------------------------
plot_impact = function(o, fig_name, plot_file = NULL, ...) {
  
  # Collate key word arguments and override geom
  args = list_modify(list(...), plot_geom = "bar")
  
  # ---- Figure properties ----
  
  # Scenarios to plot
  f = fig_scenarios(o, list(...))
  
  # Handle alternative functionality case - file fed in directly
  if (!is.null(plot_file)) baseline = plot_file
  else {
    
    # Otherwise load baseline file
    baseline = try_load(o$pth$scenarios, f$baseline_name)
  }
  
  # Collate and interpret function inputs so we know what to plot
  f = fig_properties(o, f, baseline$input, args)
  
  # If nothing to plot, return out
  if (f$n_metrics == 0)
    return()
  
  # Valid values for summarise argument
  summarise_valid = c("sum", "mean", "min", "max")
  
  # Throw an error if argument unrecognised
  if (!f$summarise %in% summarise_valid)
    stop("Argument 'summarise = ", f$summarise, "' is invalid, use one of: ", 
         paste(summarise_valid, collapse = ", "))
  
  # ---- Extract model predictions ----
  
  # Initiate plotting list
  plot_list = list()
  
  # Loop through metrics to plot
  for (scenario in f$scenarios) {
    
    # Load data file for for this scenario
    if (scenario == f$baseline_name) result = baseline else {
      result = try_load(o$pth$scenarios, scenario)
    }
    
    # Format model output and store in list to be concatenated
    plot_list[[scenario]] = format_results(o, f, result)  # See postprocess.R
  }
  
  # ---- Summarise results ----
  
  # Get summarise function from user-defined string
  #
  # NOTE: As cumulative metrics have already been computed, we can just take
  #       the last (ie the largest) value to avoid (very large) overcounting
  if (f$summarise != "sum") summarise_fn = get(f$summarise)
  if (f$summarise == "sum") summarise_fn = get("max")
  
  # Summarise each metric over time for each scenario
  plot_df = rbindlist(plot_list) %>%
    mutate(metric   = recode(metric,   !!!f$metric_names), 
           scenario = recode(scenario, !!!f$scenario_names), 
           scenario = factor(scenario, levels = f$scenario_names)) %>%
    group_by(scenario, metric, group) %>%
    summarise(value = summarise_fn(value), 
              lower = summarise_fn(lower), 
              upper = summarise_fn(upper)) %>%
    as.data.table()
  
  # Check flag for plotting difference relative to baseline
  if (f$relative == TRUE) {
    
    # Full name of baseline scenario
    baseline_name = f$scenario_names[[f$baseline_name]]
    
    # Difference betwen outcomes of each scenario and baseline
    plot_df = group_by(plot_df, metric, group) %>%
      mutate(value = value - value[scenario == baseline_name], 
             lower = lower - lower[scenario == baseline_name], 
             upper = upper - upper[scenario == baseline_name]) %>%
      filter(scenario != baseline_name) %>%
      as.data.table()
  }
  
  # Append scenario grouping details if plotting by scenario groups  
  if (f$plot_type == "scenario_group") {
    
    # We need to apply additional factoring to retain correct order
    levels_group = unique(f$scenario_groups$scenario_group)
    levels_type  = unique(f$scenario_groups$scenario_type)
    
    # Append scenario grouping details 
    plot_df = plot_df %>% 
      left_join(f$scenario_groups, by = "scenario") %>%
      mutate(scenario       = factor(scenario,       levels = f$scenario_names), 
             scenario_group = factor(scenario_group, levels = levels_group), 
             scenario_type  = factor(scenario_type,  levels = levels_type))
  }
  
  # Do we have any values less than zero (possible if plotting relative to baseline)
  below_zero = ifelse(any(plot_df$value < 0), TRUE, FALSE)
  
  # ---- Create basic figure ----
  
  # In most cases the plot_type difines the colour scheme
  aes_by = list(x = f$plot_type, fill = f$plot_type)
  
  # Slightly different aesthetic configuration if stacking by groups
  if (f$plot_type == "group" && (f$group_stack || f$group_dodge))
    aes_by$x = "scenario"
  
  # Initiate ggplot structure
  g = ggplot(plot_df, aes_string(x = aes_by$x, 
                                 y = "value", 
                                 ymin = "lower", 
                                 ymax = "upper", 
                                 fill = aes_by$fill))
  
  # Bars can be dodged
  if (f$group_dodge == TRUE) {
    bar_position = "dodge"
    err_position = position_dodge(f$bar_width)
  }
  
  # ...  or stacked
  if (f$group_dodge == FALSE) {
    bar_position = "stack"
    err_position = "identity"
  }
  
  # Only need a legend if plotting by groups
  plot_legend = ifelse(f$plot_type == "group" || f$bar_legend, TRUE, FALSE)
  
  # Plot impact bars
  g = g + geom_bar(stat     = "identity", 
                   position = bar_position, 
                   width    = f$bar_width,
                   colour   = "black",
                   size     = 0.25, 
                   show.legend = plot_legend)
  
  # Flip to horizontal bars for tornado plots
  if (f$group_tornado == TRUE)
    g = g + coord_flip()
  
  # Apply error bars if desired
  if (f$uncertainty == TRUE)
    g = g + geom_errorbar(position = err_position,
                          colour   = "darkgrey", 
                          width    = 0.25, 
                          size     = 0.25)
  
  # Plot y = 0 reference line if we've gone below zero (ie negative impact)
  if (below_zero == TRUE)
    g = g + geom_hline(yintercept = 0)
  
  # ---- Generate facets ----
  
  # Possible to define custom faceting function, in the form of evaluable string
  if (!is.null(f$facet_custom))
    eval(parse(text = paste0("g = g + ", f$facet_custom)))
  
  # Most common use case is to use facet_wrap/grid based on what is being varied
  if (is.null(f$facet_custom)) {
    
    # Metric plots (ie only one scenario) do not require facets
    if (f$plot_type != "metric") {
      
      # For scenario grouping use metric x scenario_type grid
      if (f$plot_type == "scenario_group" && !f$group_tornado) {
        g = g + facet_grid(metric~scenario_type, scales = "free", labeller = f$label_grid)
        
        # For grouping use metric x scenario grid
      } else if (f$plot_type == "group") {
        
        # If stacking or dodging by scenario, drop the scenario columns from the grid
        if (f$group_stack || f$group_dodge) {
          g = g + facet_grid(metric~., scales = "free", labeller = f$label_grid)
          
        } else {  # Otherwise apply the metric x scenario grid
          g = g + facet_grid(metric~scenario, scales = "free", labeller = f$label_grid)
        }
        
      } else { # Otherwise wrap by metric with all scenarios in each factet
        g = g + facet_wrap(~metric, scales = "free", labeller = f$label_wrap)
      }
    }
  }
  
  # Tag facets if desired
  if (f$facet_labels == TRUE)
    g = facet_labels(g)
  
  # ---- Figure aesthetics ----
  
  # If plotting by group, extract dictionary (age already considered)
  if (f$plot_type == "group" && f$plot_by != "age") {
    use_labels = baseline$input$dict[[f$plot_by]]
    
  } else {  # Otherwise no additional labels needed
    use_labels = waiver()
  }
  
  # Set colour scheme
  g = g + scale_fill_manual(values = f$colours, labels = use_labels)
  
  # Prettify y axis - expand above by default, but only below if we have negative values
  if (below_zero) y_expand = rep(f$y_expand, 2) else y_expand = c(0, f$y_expand)
  g = g + scale_y_continuous(labels = comma, expand = expansion(mult = y_expand))
  
  # Prettify x axis by wrapping long lines
  g = g + scale_x_discrete(labels = wrap_format(f$x_wrap))
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(plot.title    = element_text(size = f$fontsize[1], hjust = 0.5),
          strip.text    = element_text(size = f$fontsize[2], margin = margin(2, 2, 2, 2)),
          axis.text.y   = element_text(size = f$fontsize[4]),
          axis.text.x   = element_text(size = f$fontsize[4]),
          axis.title    = element_blank(),
          axis.line     = element_blank(),
          axis.ticks    = element_line(size = 0.25),
          axis.ticks.length = unit(0.1, "lines"),
          panel.border  = element_rect(size = 0.5, colour = "black", fill = NA),
          panel.spacing = unit(0.5, "lines"),
          strip.background = element_blank(), 
          legend.text   = element_text(size = f$fontsize[3]),
          legend.title  = element_blank(),
          legend.key    = element_blank(),
          legend.margin = margin(2, 1, 2, 1),
          legend.position = "bottom", 
          legend.key.height = unit(1, "lines"),
          legend.key.width  = unit(1, "lines"),
          legend.box.background = element_rect(size = 0.5))
  
  # Rotate x axis labels if desired
  if (f$x_rotate != 0)
    g = g + theme(axis.text.x = element_text(hjust = 1, angle = f$x_rotate))
  
  # Apply grid line for tornado plots only
  if (f$group_tornado == TRUE)
    g = g + theme(panel.grid.major.x = element_line(), 
                  panel.grid.minor.x = element_line())
  
  # Remove tick labels if necessary (ie facets and legend tell the whole story)
  if (plot_legend == TRUE) {
    
    # We generally remove x axis unless stacking or dodging by group
    keep_xaxis1 = f$plot_type == "group" && (f$group_stack || f$group_dodge)
    keep_xaxis2 = f$group_tornado == TRUE
    
    # Normal use case is to remove the x axis labels...
    if (!keep_xaxis1 && !keep_xaxis2)
      g = g + theme(axis.text.x  = element_blank(), 
                    axis.ticks.x = element_blank())
    
    # ... but for tornado plots it's the y axis labels we want to erase
    if (f$group_tornado == TRUE)
      g = g + theme(axis.text.y  = element_blank(), 
                    axis.ticks.y = element_blank())
    
    # Set appropriate number of legend rows
    g = g + guides(fill = guide_legend(nrow  = f$legend_rows, 
                                       byrow = f$legend_by_row))
  }
  
  # Save these figures to file
  fig_save(o, g, fig_name)
  
  return(g)
} 

# ---------------------------------------------------------
# Plot heat maps for multi-dimensional array scenarios
# ---------------------------------------------------------
plot_heatmap = function(o, fig_name, array, plot_df = NULL, ...) {
  
  # Collate key word arguments and override geom
  args = list_modify(list(...), plot_geom = "tile", override_colours = NA)
  
  # ---- Figure properties ----
  
  # Scenarios to plot
  f = fig_scenarios(o, args)
  
  # Load baseline (or alt baseline) file
  baseline = try_load(o$pth$scenarios, f$baseline_name)
  
  # Collate and interpret function inputs so we know what to plot
  f = fig_properties(o, f, baseline$input, args)
  
  # If nothing to plot, return out
  if (f$n_metrics == 0)
    return()
  
  # Multiple metrics are OK if plotting 2D arrays only
  if (f$n_metrics > 1)
    stop("Multiple metrics are not yet possible")
  
  # Valid values for summarise argument
  summarise_valid = c("sum", "mean", "min", "max")
  
  # Throw an error if argument unrecognised
  if (!f$summarise %in% summarise_valid)
    stop("Argument 'summarise = ", f$summarise, "' is invalid, use one of: ", 
         paste(summarise_valid, collapse = ", "))
  
  # ---- Inspect array details ----
  
  # Determine number of dimensons from array_info
  if (is.null(plot_df)) {
    
    # Load details of this array scenario
    array_info = try_load(o$pth$arrays, array)
    
    # Values of each dimension simulated
    array_df = array_info$values %>% rename(primary = scenario)
    
    # Extract key array scenario details
    array_scenarios = array_info$meta$scenario
    n_dims = nrow(array_info$vars)
    
  } else {  # If plotting dataframe is already provided
    
    # Number of dimensions defined by plot_df columns
    n_dims = length(intersect(names(plot_df), qc(x, y, w, v)))
  }
  
  # Sanity checks: heat maps only support upto 4 dimensions for now
  if (n_dims < 2) stop("Heat maps require array scenarios to be of at least 2 dimensions")
  if (n_dims > 4) stop("Heat maps are currently only able to support up to 4 dimensions")
  
  # ---- Construct plotting dataframe ----
  
  # Skip this entire process if plot_df is already provided
  if (is.null(plot_df)) {
    
    # ---- Inspect details of any 'difference' array ----
    
    # Cleaner to use simple flags for dif plotting
    if (!is.null(f$dif_array)) plot_dif = TRUE else plot_dif = FALSE
    
    # If plotting the difference to another array, some extration needed...
    if (plot_dif) {
      
      # Dif arrays should be defined by <analysis_name>::<array_name>
      dif_split = strsplit(f$dif_array, "::")[[1]]
      
      # Split <analysis_name> and <array_name>
      dif_name = rev(dif_split)[1]
      dif_analysis = dif_split[1]
      
      # If array name is not given, assume it is the current analysis name
      if (length(dif_split) != 2)
        dif_analysis = o$analysis_name
      
      # Paths to dif array info and scenario files
      dif_info_pth = pth_replace(o$pth$arrays,    o$analysis_name, dif_analysis)
      dif_scen_pth = pth_replace(o$pth$scenarios, o$analysis_name, dif_analysis)
      
      # Load dif array info
      dif_array_info = try_load(dif_info_pth, dif_name)
      
      # Extract values - we need to check for consistency
      dif_array_values = dif_array_info$values %>%
        mutate(sub = str_remove(scenario, paste0("^", dif_name, "\\."))) %>%
        rename(difference = scenario, dif_value = value)
      
      # Join these with equivelent from 'primary' array
      dif_df = array_df %>%
        mutate(sub = str_remove(primary, paste0("^", array, "\\."))) %>%
        left_join(dif_array_values, by = c("sub", "variable_id"))
      
      # Check simualted values along each array dimension are identical
      if (any(abs(dif_df$value - dif_df$dif_value) > 1e-6))
        stop("Dif array not consistent with primary array")
      
      # Vector of scenarios from dif array
      dif_scenarios = unique(dif_df$difference)
      
      # Format array df to be consistent if NOT plotting against a dif array
      array_df = dif_df[, .(primary, difference, variable_id, value)]
    }
    
    # ---- Extract model predictions ----
    
    # Initiate plotting list
    result_list = list()
    
    # Initiate progress bar
    pb = start_progress_bar(length(array_scenarios))
    
    # Loop through scenarios to plot
    for (i in 1 : length(array_scenarios)) {
      scenario = array_scenarios[i]
      
      # Attempt to load scenario outcomes
      result = try_load(o$pth$scenarios, scenario, throw_error = FALSE)
      
      # If no result found, we'll just throw a warning
      if (is.null(result)) {
        warning("Scenario '", scenario, "' not found - plotting NA")
        
      } else {  # Otherwise continue...
        
        # Extract what we need (see postprocess.R)
        result = format_results(o, f, result) %>%
          mutate(array_type = "primary")
        
        # Store in list to be concatenated
        result_list[[scenario]] = result
      }
      
      # If we are plotting dif, also load and store this
      if (plot_dif) {
        dif_scenario = dif_scenarios[i]
        
        # Attempt to load dif scenario outcomes
        dif_result = try_load(dif_scen_pth, dif_scenario, throw_error = FALSE)
        
        # If no result found, we'll just throw a warning
        if (is.null(dif_result)) {
          warning("Dif scenario '", dif_scenario, "' not found - plotting NA")
          
        } else {   # Otherwise continue...
          
          # Extract what we need (see postprocess.R)
          dif_result = format_results(o, f, dif_result) %>%
            mutate(array_type = "difference")
          
          # Store in list to be concatenated
          result_list[[paste0("dif_", dif_scenario)]] = dif_result
        }
      }
      
      # Update progress bar
      setTxtProgressBar(pb, i)
    }
    
    # Close progress bar
    close(pb)
    
    # Get summarise function from user-defined string
    summarise_fn = get(f$summarise)
    
    # Summarise results over time (summary function as defined by 'summarise')
    result_df = rbindlist(result_list) %>%
      group_by(array_type, scenario, group) %>%
      summarise(value = summarise_fn(value)) %>%
      as.data.table()
    
    # If plotting by some grouping, we need to select just one
    if (f$plot_type == "group") {
      
      # Select this group (using group_idx)
      plot_group = levels(result_df$group)[f$group_idx]
      result_df = result_df[group == plot_group, ]
    }
    
    # Check flag for plotting difference relative to baseline
    if (f$relative == TRUE) {
      
      # A bit more work needed to test this eventuality
      if (plot_dif == TRUE)
        stop("Plotting difference array AND relative values has not been properly tested")
      
      # Difference betwen outcomes of each scenario and baseline
      result_df[, value := value - result_df[scenario == f$baseline_name, value]]
    }
    
    # ---- Split array variables ----
    
    # Dimension IDs used for ggplot variables and facets
    dim_id = qc(x, y, w, v)[c(1, 2, 2 + seq_len(n_dims - 2))]
    
    # Associated names for the array variables (used for axes labels)
    dim_names = first_cap(array_info$vars$variable_name)
    
    # Format results into plottable dataframe
    plot_df = array_df %>%
      pivot_wider(id_cols    = c(-variable_id, -value),
                  names_from = variable_id) %>%
      rename(setNames(array_info$vars$variable_id, dim_id)) %>%
      pivot_longer(cols      = -all_of(dim_id), 
                   names_to  = "array_type", 
                   values_to = "scenario") %>%
      left_join(result_df, by = c("array_type", "scenario")) %>%
      pivot_wider(id_cols    = c(all_of(dim_id)), 
                  names_from = array_type) %>%
      as.data.table()
    
    # Set dif as trivial if not plotting dif array
    if (!plot_dif)
      plot_df[, difference := 0]
    
    # Plot value can then be taken as the diference
    plot_df[, value := primary - difference]
    
    # Calculate relative to primary outcomes if desired
    if (f$dif_relative)
      plot_df[, value := value / primary]
    
    # Convert 3rd and/or 4th dimensions to percentages if required
    if (n_dims >= 3 && f$w_percent == TRUE) plot_df$w = paste0(plot_df$w * 100, "%")
    if (n_dims >= 4 && f$v_percent == TRUE) plot_df$v = paste0(plot_df$v * 100, "%")
    
    # For higher dimensions, use discrete strings rather than continuous values
    if (n_dims >= 3) plot_df$w = as.factor(paste0(dim_names[3], "=", plot_df$w))
    if (n_dims >= 4) plot_df$v = as.factor(paste0(dim_names[4], "=", plot_df$v))
    
    # ---- Interpolation ----
    
    # Check whether we want to do linear interpolation of data
    if (!is.null(f$n_interpolate)) {
      interp_list = list()
      
      # For higher dimensions, use discrete strings rather than continuous values
      if (n_dims < 3) plot_df$w = 1
      if (n_dims < 4) plot_df$v = 1
      
      # Loop through higher-dimensions - needs to be 2D at a time
      for (this_v in unique(plot_df$v)) {
        for (this_w in unique(plot_df$w)) {
          
          # Data for this facet only
          facet_df = filter(plot_df, v == this_v, w == this_w)
          
          # Our interpolation function cannot handle NAs - there should not be any anyway
          if (all(!is.na(facet_df$value))) {
            
            # Linearly interpolate z values for n_interpolate points within axes limits 
            interp_matrix = interp(x = facet_df$x, y = facet_df$y, z = facet_df$value, 
                                   nx = f$n_interpolate, ny = f$n_interpolate)
            
            # Melt down to dataframe long format
            interp_df = pivot_longer(as.data.table(interp_matrix$z), cols = everything())
            
            # Concatenate x and y with interpolated values and append facet w and v values
            interp_list[[paste0(this_v, "_", this_w)]] = 
              expand_grid(x = interp_matrix$x, 
                          y = interp_matrix$y) %>%
              mutate(w = this_w, v = this_v, 
                     value = interp_df$value)
            
          } else {  # Throw a warning if unable to interpolate
            warning("Unable to interpolate facet due to NA values")
          }
        }
      }
      
      # Throw an error if no facets could be interpolated
      if (length(interp_list) == 0)
        stop("There are too many NAs in your results to perform interpolation")
      
      # Concatenate together and overwrite plotting dataframe
      plot_df = rbindlist(interp_list) %>%
        mutate(w = factor(w, levels = levels(plot_df$w)), 
               v = factor(v, levels = levels(plot_df$v)))
    }
    
    # Check if user has colour limit restrictions
    if (!is.null(f$colour_limits)) {
      
      # Bound values above min colour_limit
      if (!is.na(f$colour_limits[1]))
        plot_df[, value := pmax(value, f$colour_limits[1])]
      
      # Bound values below max colour_limit
      if (!is.na(f$colour_limits[2]))
        plot_df[, value := pmin(value, f$colour_limits[2])]
    }
  }
  
  # ---- Produce plot ----
  
  # Initiate ggplot
  g = ggplot(plot_df, aes(x = x, y = y))
  
  # Set factets if more than 2 dimensions
  if (n_dims == 3) g = g + facet_wrap(~w,  labeller = f$label_wrap)
  if (n_dims == 4) g = g + facet_grid(w~v, labeller = f$label_grid)
  
  # Normal use case: single colourbar
  if (is.null(f$multi_scale)) {
    
    # Plot full dataframe 
    g = g + geom_tile(aes(fill = value))
    
    # Set colour scheme
    g = apply_colourscale(g, f, 1, plot_df)
    
  } else {  # Also able to plot multiple colourbars across one dimension
    
    # Iterate through unique values across chosen dimension
    unique_scales = levels(plot_df[[f$multi_scale]])
    for (i in seq_along(unique_scales)) {
      
      # Subset plotting dataframe
      scale_df = filter(plot_df, !!as.name(f$multi_scale) == unique_scales[i])
      
      # Plot facets for this subset of data only
      g = g + geom_tile(data = scale_df, aes(fill = value))
      
      # Set unique colour scheme (colour_palette a vector, colour_limits an nx2 matrix)
      g = apply_colourscale(g, f, i, scale_df)
      
      # Prepare to add a new scale if more to come
      if (i < length(unique_scales))
        g = g + new_scale_fill()
    }
  }
  
  # Add user-specified contours if desired
  if (!is.na(f$contour_val)) {
    
    # Points to fit to - values within tolerance of target
    contour_df = plot_df %>% 
      filter(value < f$contour_val * (1 + f$contour_tol), 
             value > f$contour_val * (1 - f$contour_tol))
    
    # Contours can have different aesthetics along a dimension if desired (using contour$by)
    if (is.na(f$contour_by)) contours_by = NA
    else contours_by = unique(contour_df[[f$contour_by]])
    
    # Loop through unique contours to plot (trivial if contour$by = NULL)
    for (i in seq_along(contours_by)) {
      
      # Subset contour dataframe (if necessary)
      if (is.na(f$contour_by)) this_contour = contour_df
      else this_contour = filter(contour_df, !!as.name(f$contour_by) == contours_by[i])
      
      # Plot a smoothed contour beyond bounds of subplot
      g = g + geom_smooth(data      = this_contour, 
                          mapping   = aes(x = x, y = y), 
                          method    = "lm", 
                          formula   = f$contour_fn, 
                          fullrange = TRUE, 
                          se        = FALSE, 
                          linetype  = length(contours_by) - i + 1, 
                          size      = 2, 
                          colour    = "black")
    }
    
    # Setting limits within scale_x_continuous doesn't work - this is a workaround
    g = g + coord_cartesian(xlim = c(min(plot_df$x), max(plot_df$x)), 
                            ylim = c(min(plot_df$y), max(plot_df$y)))
  }
  
  # ---- Set title and axes labels ----
  
  # Turn titles off completely if user specifies NA
  if (is.null(f$plot_title))    f$plot_title = waiver()
  if (is.null(f$plot_subtitle)) f$plot_subtitle = waiver()
  
  # Plot title and subtitle
  g = g + labs(title = f$plot_title, subtitle = f$plot_subtitle)
  
  # Plot x label if desired
  if (!is.null(f$x_lab)) g = g + xlab(f$x_lab) else {
    if (exists("dim_names")) g = g + xlab(dim_names[1]) 
    else g = g + xlab(NULL)
  }
  
  # Plot y label if desired
  if (!is.null(f$y_lab)) g = g + ylab(f$y_lab) else {
    if (exists("dim_names")) g = g + ylab(dim_names[2])
    else g = g + ylab(NULL)
  }
  
  # Prettify axes
  g = g + scale_x_continuous(expand = c(0, 0), labels = f$x_scale_fn) +
    scale_y_continuous(expand = c(0, 0), labels = f$y_scale_fn)
  
  # Prettify theme
  g = g + theme_classic() +
    theme(plot.title    = element_text(size = f$fontsize[1], hjust = 0.5),
          plot.subtitle = element_text(size = f$fontsize[2], hjust = 0.5),
          strip.text.x  = element_text(size = f$fontsize[4]),
          strip.text.y  = element_text(size = f$fontsize[5]),
          legend.title  = element_blank(),
          legend.position = "top", 
          legend.key.width = unit(4, "cm"), 
          legend.text   = element_text(size = f$fontsize[6]),
          axis.title    = element_text(size = f$fontsize[3]),
          axis.text     = element_text(size = f$fontsize[7]),
          axis.line     = element_blank(),
          panel.border  = element_rect(size = 1, colour = "black", fill = NA), 
          panel.spacing.x = unit(1.4, "lines"), 
          panel.spacing.y = unit(1.2, "lines"))
  
  # ---- Save figure and plotting details ----
  
  # As this can be expensive, also save plotting info
  plot_info = list(plot_df = plot_df, g = g)
  
  # Save contour points if necessary
  if (!is.na(f$contour_val))
    plot_info = list.append(plot_info, 
                            contour_df  = contour_df, 
                            contour_val = f$contour_val, 
                            contour_tol = f$contour_tol, 
                            contour_fn  = f$contour_fn, 
                            contour_by  = f$contour_by)
  
  # File name to save to (similar, but computer readable)
  plot_file = tolower(paste(fig_name, collapse = "_")) %>%
    str_replace_all(" ", "_")
  
  # Save this to file
  saveRDS(plot_info, paste0(o$pth$figures, plot_file, ".rds"))
  
  # Finally, save this figure to file
  fig_save(o, g, fig_name)
  
  return(plot_info)
} 

# ---------------------------------------------------------
# Plot number of infections for multiple scenarios
# ---------------------------------------------------------
plot_num_infections = function(o, fig_name, plot_file = NULL, ...) {
  
  # Collate key word arguments and override several values
  args = list_modify(list(...), 
                     plot_geom    = "hist", 
                     plot_metrics = "n_infections", 
                     plot_by      = "infections")
  
  # ---- Figure properties ----
  
  # Scenarios to plot
  f = fig_scenarios(o, list(...))
  
  # Handle alternative functionality case - file fed in directly
  if (!is.null(plot_file)) baseline = plot_file
  else {
    
    # Otherwise load baseline file
    baseline = try_load(o$pth$scenarios, f$baseline_name)
  }
  
  # Collate and interpret function inputs so we know what to plot
  f = fig_properties(o, f, baseline$input, args)
  
  # If nothing to plot, return out
  if (f$n_metrics == 0)
    return()
  
  # ---- Extract model predictions ----
  
  # Initiate plotting list
  plot_list = list()
  
  # Loop through metrics to plot
  for (scenario in f$scenarios) {
    
    # Load data file for for this scenario
    if (scenario == f$baseline_name) result = baseline else {
      result = try_load(o$pth$scenarios, scenario)
    }
    
    # Format model output and store in list to be concatenated
    plot_list[[scenario]] = format_results(o, f, result)  # See postprocess.R
  }
  
  # Concatenate dataframes and convert to percentages for all scenarios
  plot_df = rbindlist(plot_list) %>%
    group_by(scenario) %>%
    mutate(lower = lower / sum(value), 
           upper = upper / sum(value), 
           value = value / sum(value)) %>%
    ungroup() %>%
    mutate(group    = as.factor(paste0(group, " infections")),
           metric   = recode(metric,   !!!f$metric_names), 
           scenario = recode(scenario, !!!f$scenario_names), 
           scenario = factor(scenario, levels = f$scenario_names)) %>%
    as.data.table()
  
  # ---- Create basic figure ----
  
  # X axes is the complement of what we're faceting: scenarios or groups
  x_aes = setdiff(c("scenario", "group"), f$facet_by)
  
  # Make sure one or teh other has been defined
  if (length(x_aes) != 1)
    stop("Argument 'facet_by' must be either 'scenario' or 'group'")
  
  # Initiate plot and aesthetics 
  g = ggplot(plot_df, aes_string(x    = x_aes, 
                                 y    = "value", 
                                 ymin = "lower", 
                                 ymax = "upper",
                                 fill = "scenario"))
  
  # Only need a legend if plotting by groups
  plot_legend = ifelse(f$facet_by == "group", TRUE, FALSE)
  
  # Plot histogram, we can use geom_bar here as we've already binned
  g = g + geom_bar(stat = "identity", colour = "black", show.legend = plot_legend)
  
  # Apply error bars if desired
  if (f$uncertainty == TRUE)
    g = g + geom_errorbar(colour = "darkgrey", width = 0.25, size = 0.5)
  
  # Apply faceting function
  g = g + facet_wrap(as.formula(paste("~", f$facet_by)), labeller = f$label_wrap)
  
  # Tag facets if desired
  if (f$facet_labels == TRUE)
    g = facet_labels(g)
  
  # ---- Figure aesthetics ----
  
  # Set colour scheme
  g = g + scale_fill_manual(values = f$colours)
  
  # Prettify y axis - expand above by default
  g = g + scale_y_continuous(labels = percent, expand = expansion(mult = c(0, f$y_expand)))
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(plot.title   = element_text(size = f$fontsize[1], hjust = 0.5),
          strip.text   = element_text(size = f$fontsize[2]),
          axis.text.y  = element_text(size = f$fontsize[4]),
          axis.text.x  = element_text(size = f$fontsize[4], hjust = 1, angle = 50),
          axis.title   = element_blank(),
          axis.line    = element_blank(),
          axis.ticks    = element_line(size = 0.25),
          axis.ticks.length = unit(0.1, "lines"),
          panel.border = element_rect(size = 1, colour = "black", fill = NA),
          panel.spacing = unit(1, "lines"),
          strip.background = element_blank(), 
          legend.text  = element_text(size = f$fontsize[3]),
          legend.title = element_blank(),
          legend.key   = element_blank(),
          legend.position = "bottom",
          legend.key.height = unit(2, "lines"),
          legend.key.width  = unit(2, "lines"),
          legend.box.background = element_rect())
  
  # Remove tick labels if necessary (ie facets and legend tell the whole story)
  if (plot_legend == TRUE) {
    g = g + theme(axis.text.x  = element_blank(), 
                  axis.ticks.x = element_blank())
    
    # Set appropriate number of legend rows
    g = g + guides(fill = guide_legend(nrow  = f$legend_rows, 
                                       byrow = f$legend_by_row))
  }
  
  # Save these figures to file
  fig_save(o, g, fig_name)
  
  return(g)
}

# ---------------------------------------------------------
# Plot normalised disease state of population over time
# ---------------------------------------------------------
plot_disease_state = function(o, fig_name, p, states) {
  
  # Convert list to dataframe and remove what we're not interested in
  state_df = states[p$model_states$all] %>%
    as.data.table() %>%
    mutate(date = 1 : p$n_days) %>%
    filter(susc > 0) %>%
    select(-susc, -none)
  
  # Convert to long format and preserve order with factors
  plot_df = state_df %>%
    pivot_longer(cols = -date,
                 names_to  = "state", 
                 values_to = "value") %>%
    mutate(state_type = ifelse(state %in% p$model_states$care, "care", "disease"), 
           state_type = factor(state_type, levels = c("disease", "care")), 
           state      = factor(state,      levels = p$model_states$all)) %>%
    arrange(state_type, state, date) %>%
    as.data.table()
  
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
  fig_save(o, g1, fig_name, "Stacked area")
  fig_save(o, g2, fig_name, "Line")
}

# ---------------------------------------------------------
# Plot performance of Gaussian process emulators
# ---------------------------------------------------------
plot_emulator = function(o, fig_name) {
  
  # Data types we're working with
  data_types = c("train", "test")
  
  # Corresponding colour scheme
  colours = c("skyblue1", "blue2")
  
  # Flag to remove extreme values
  remove_outliers = FALSE
  
  # ---- Construct plot dataframe ----
  
  # Load emulator - throw an error if it doesn't exist explaining step 1 needs to run first
  err_msg = "Attempting to plot emulator performance but cannot find model emulator file"
  emulator = try_load(o$pth$fitting, "model_emulator", msg = err_msg) 
  
  # Extract plotting dataframe
  plot_df = emulator$performance %>%
    mutate(group = factor(group, levels = data_types))
  
  # Remove outliers if desired (> 3 std devs from mean)
  if (remove_outliers == TRUE)
    plot_df = filter(plot_df, 
                     abs(actual  - mean(actual))  < 3 * sd(actual), 
                     abs(predict - mean(predict)) < 3 * sd(predict))
  
  # Construct dummy dataframe so we get 
  dummy_df = select(plot_df, actual = predict, predict = actual, group)
  
  # ---- Produce plot ----
  
  # Plot truth vs predicted (also plot invisible in reverse for square axes)
  g = ggplot(plot_df, aes(x = actual, y = predict, colour = group)) + 
    geom_blank(data = dummy_df) +
    geom_point(size = 5, alpha = 0.8) + 
    geom_abline()
  
  # Prettify axes
  g = g + labs(x = "Actual", y = "Predicted", title = "Emulator performance") + 
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) + 
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)))
  
  # Prettify legend - depends on type of plot
  g = g + scale_colour_manual(name   = "Data type", 
                              values = colours, 
                              labels = first_cap(data_types))
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(plot.title   = element_text(size = 24, hjust = 0.5), 
          axis.title   = element_text(size = 20), 
          axis.text.x  = element_text(size = 12), 
          axis.text.y  = element_text(size = 12), 
          legend.text  = element_text(size = 16), 
          legend.title = element_text(size = 18), 
          axis.line    = element_blank(), 
          panel.border = element_rect(colour = "black", fill = NA, size = 1))
  
  # Save figure to file
  fig_save(o, g, fig_name)
}

# ---------------------------------------------------------
# Plot optimisation performance and parameter-R_eff relationship
# ---------------------------------------------------------
plot_optimisation = function(o, fig_name) {
  
  # Optimisation colour
  colour_optim = "dodgerblue"
  
  # Colour map for seeds
  colour_map = "viridis::viridis"
  
  # Number of samples for emulator plotting
  n_eval = 100
  
  # Expad axes limits
  ax_expand = 0.01
  
  # ---- Load samples, emulator, and optimisation result ----
  
  # Load emulator - throw an error if it doesn't exist explaining step 1 needs to run first
  err_msg = "Attempting to plot emulator performance but cannot find model emulator file"
  emulator = try_load(o$pth$fitting, "model_emulator", msg = err_msg) 
  
  # Also load simulated samples and optimisation result
  samples_df = try_load(o$pth$fit_samples, "all_samples")
  fit_result = try_load(o$pth$fitting,     "fit_result")
  
  # Shorthand for fitted parameter names and bounds
  fit_params = fit_result$params
  fit_bounds = fit_result$bounds
  
  # ---- Evaluate emulator at a load of points ----
  
  # Bounds of params-dimensional unit cube
  unit_cube = matrix(rep(c(0, 1), each = length(fit_params)), ncol = 2)
  
  # Sample a load of points and evaluate using emulator
  eval_input  = lhs(n_eval, unit_cube) 
  eval_output = predict(emulator, x = eval_input)
  
  # Function for extacting prediction bounds for emulator result
  eval_bounds_fn = function(x, q) qnorm(q, x$mean, sqrt(x$sd2 + x$nugs))
  
  # Extract emulator evaluated best estimate and bounds
  eval_df = as_named_dt(eval_input, fit_params) %>%
    mutate(mean  = eval_output$mean, 
           lower = eval_bounds_fn(eval_output, 0.05), 
           upper = eval_bounds_fn(eval_output, 0.95)) 
  
  # Retransform evaluated points back to real scale
  for (param in fit_params) {
    param_idx = which(fit_params %in% param)
    
    # Use parameter bounds rather than bounds of sampled points
    eval_df[[param]] = normalise_0to1(
      x = eval_df[[param]], 
      x_min = fit_bounds[[param_idx, 1]], 
      x_max = fit_bounds[[param_idx, 2]], 
      direction = "backward")
  }
  
  # ---- Construct plotting dataframes ----
  
  # Melt down simulated samples ready for plotting
  simulated_df = samples_df %>%
    select(-param_id) %>%
    pivot_longer(cols = all_of(fit_params), 
                 names_to = "param") %>%
    mutate(seed = as.factor(seed)) %>%
    as.data.table()
  
  # Melt down emulator evaluations ready for plotting
  emulator_df = eval_df %>%
    pivot_longer(cols = all_of(fit_params), 
                 names_to = "param") %>%
    as.data.table()
  
  # Construct best optimisation results dataframe
  best_df = fit_result$result %>%
    as.data.table() %>% 
    pivot_longer(cols = everything(), 
                 names_to  = "param") %>%
    as.data.table()
  
  # Compile with bounds of oiptimisation process
  optim_df = fit_result$output$x %>%
    as_named_dt(fit_params) %>%
    pivot_longer(cols = everything(), 
                 names_to = "param") %>%
    group_by(param) %>%
    summarise(lower = min(value), 
              upper = max(value)) %>%
    left_join(best_df, by = "param") %>%
    as.data.table()
  
  # ---- Produce plot ---
  
  # Plot simulated points (distinguished by seed)
  g = ggplot(simulated_df, aes(x = value)) +
    geom_point(aes(y = r_eff, colour = seed), 
               size = 2, stroke = 0, alpha = 0.6) +
    facet_wrap(~param)
  
  # Plot emulator best fit and bounds
  g = g + geom_ribbon(data = emulator_df, 
                      mapping = aes(ymin = lower, ymax = upper), alpha = 0.4) +
    geom_line(data = emulator_df, mapping = aes(y = mean), size = 2) 
  
  # Plot target R_eff
  g = g + geom_hline(yintercept = fit_result$input$r_eff, 
                     colour = colour_optim, linetype = "dashed", size = 2)
  
  # Plot optimal params to achieve this R_eff
  g = g + geom_vline(data = optim_df, aes(xintercept = value), 
                     colour = colour_optim, size = 2)
  
  # Finally, plot extremes found from the different optimisation runs
  g = g + geom_rect(data = optim_df, 
                    mapping = aes(xmin = lower, ymin = -Inf, 
                                  xmax = upper, ymax = Inf), 
                    fill = colour_optim, color = NA, alpha = 0.2)
  
  # ---- Prettify plot ----
  
  # Prettify axes
  g = g + labs(x = "Parameter value", 
               y = "Effective reproduction number", 
               title = "Optimisation performance") + 
    scale_x_continuous(expand = expansion(mult = c(ax_expand, ax_expand))) + 
    scale_y_continuous(expand = expansion(mult = c(ax_expand, ax_expand)))
  
  # Create colour scheme for different seeds
  n_seeds = length(levels(simulated_df$seed))
  seed_colours = colour_scheme(colour_map, n = n_seeds)
  
  # Set this colour scheme
  g = g + scale_colour_manual(name = "Seed number", values = seed_colours)
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(plot.title   = element_text(size = 24, hjust = 0.5), 
          strip.text   = element_text(size = 18), 
          axis.title   = element_text(size = 20), 
          axis.text.x  = element_text(size = 12), 
          axis.text.y  = element_text(size = 12), 
          legend.text  = element_text(size = 16), 
          legend.title = element_text(size = 18), 
          axis.line    = element_blank(), 
          panel.border = element_rect(colour = "black", fill = NA, size = 1))
  
  # Save figure to file
  fig_save(o, g, fig_name)
}

# ---------------------------------------------------------
# Plot series of network-related properties
# ---------------------------------------------------------
plot_network_properties = function(o, fig_name, model_input, network) {
  
  # ---- Figure properties ----
  
  # Bin width for age groups
  age_bin = 10 
  
  # Quantiles to display on violin plots
  violin_quantiles = c(0.5)  # eg c(0.25, 0.5, 0.75) for median and IQR
  
  # Set font sizes (axes labels, axes text, legend text)
  fontsize = c(16, 12, 13)
  
  # ---- Construct plotting dataframe ----
  
  # Number of contacts per person 
  count_contact = table(network$from)
  
  # Seperate into IDs and number of contacts
  id = as.numeric(names(unlist(count_contact)))
  n_contacts = as.numeric(unlist(count_contact))
  
  # All possible ages (from model input)  
  ages = model_input$age$all
  
  # Breaks to create age bins
  age_breaks = seq(0, max(ages) + age_bin, by = age_bin)
  
  # Construct age dataframe to be joined to contact count
  age_df = unique(select(network, from, age = from_age)) %>%
    mutate(age_group = cut(age, age_breaks, include.lowest = TRUE))
  
  # Join number of contacts to age (and age group) of each person
  plot_df = data.table(from = id, n_contacts = n_contacts) %>%
    full_join(age_df, by = "from") %>%
    arrange(from)
  
  # Also construct dataframe of total number in each age group
  demog_df = data.table(age = ages, n = model_input$demography) %>%
    mutate(age_group = cut(age, age_breaks, include.lowest = TRUE)) %>%
    group_by(age_group) %>%
    summarise(n = sum(n)) %>%
    as.data.table()
  
  # ---- Fig a) Number of contacts per person ----
  
  # Construct subplot
  g1 = ggplot(plot_df, aes(x = n_contacts)) + 
    geom_density(aes(y = ..count.., colour = age_group, fill = age_group), 
                 adjust = 2, alpha = 0.1, size = 1.5) + 
    geom_vline(xintercept = model_input$contacts, size = 1.5, 
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
    geom_vline(xintercept = model_input$contacts, size = 1.5, 
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
  g = ggarrange(g1, g2, g3, g4, nrow = 2, ncol = 2, align = "hv", 
                legend.grob = g_legend, legend = "bottom")
  
  # Save figure
  fig_save(o, g, fig_name)
}

# ---------------------------------------------------------
# Age-structured contact matrix heat maps
# ---------------------------------------------------------
plot_contact_matrices = function(o, fig_name, model_input, network) {
  
  # ---- Figure properties ----
  
  # Bin width for age groups
  age_bin = 5
  
  # Set font sizes (title, facets, and ticks)
  fontsize = c(30, 18, 12)
  
  # Two scale colour bars
  colour_low  = "grey80"
  colour_mid  = "dodgerblue"
  colour_high = "blue"
  
  # ---- Figure set up ----
  
  # Maximum possible age (from model input)  
  max_age = max(model_input$age$all)
  
  # Breaks to create age bins
  age_breaks = seq(0, max_age + age_bin, by = age_bin)
  
  # Layers in this network
  all_layers = unique(network$layer)
  
  # Append (or overwrite with) all contacts
  if (length(all_layers) == 1) all_layers = "all"
  else all_layers = c(all_layers, "all")
  
  # All layers combined - we'll always plot this
  edgelist_df = network %>% 
    mutate(layer = "all")
  
  # If constructed of multiple layers, plot these too
  if (length(all_layers) > 1) 
    edgelist_df = rbind(edgelist_df, network)
  
  # ---- Construct plotting dataframe ----
  
  # Specify age groups for all contacts in edgelist
  age_bin_df = edgelist_df %>%
    mutate(age1 = cut(from_age, age_breaks, include.lowest = TRUE), 
           age2 = cut(to_age,   age_breaks, include.lowest = TRUE))
  
  # Vector of all age bins
  age_bins = levels(age_bin_df$age1)
  
  # Summarise for count and proportion of contacts between age groups
  grouped_df = age_bin_df %>%
    group_by(layer, age1, age2) %>%
    summarize(count = n()) %>%
    group_by(layer) %>%
    mutate(proportion = count / max(count)) %>%
    ungroup() %>%
    mutate(count = count / max(count)) %>%
    pivot_longer(cols = c(count, proportion), 
                 names_to = "type") %>%
    as.data.table()
  
  # Fill missing with NA for prettier figures
  plot_df = grouped_df %>%
    full_join(
      expand_grid(
        layer = all_layers, 
        age1 = age_bins,
        age2 = age_bins, 
        type = c("count", "proportion")), 
      by = qc(layer, age1, age2, type)) %>%
    mutate(layer = factor(layer, levels = all_layers), 
           age1 = factor(age1, levels = age_bins), 
           age2 = factor(age2, levels = age_bins)) %>%
    arrange(type, layer)
  
  # ---- Create figure ----
  
  # Plot contact matrices for each layer
  g = ggplot(plot_df, aes(x = age1, y = age2)) + 
    geom_tile(aes(fill = value), show.legend = FALSE) + 
    facet_grid(type~layer)
  
  x_breaks = age_bins[seq(1, length(age_breaks), by = 2)]
  
  # Prettyify axes
  g = g + labs(title = paste0("Age-structure contact matrix by network layer")) + 
    scale_fill_gradient2(low = colour_low, mid = colour_mid, high = colour_high, 
                         na.value = colour_low, midpoint = 0.5, expand = c(0, 0)) + 
    scale_x_discrete(expand = c(0, 0), breaks = x_breaks) + 
    scale_y_discrete(expand = c(0, 0))
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(plot.title   = element_text(size = fontsize[1], hjust = 0.5), 
          strip.text   = element_text(size = fontsize[2]), 
          axis.title   = element_blank(), 
          axis.text.x  = element_text(size = fontsize[3], hjust = 1, angle = 50),
          axis.text.y  = element_text(size = fontsize[3]), 
          panel.border = element_rect(size = 1, colour = "black", fill = NA))
  
  
  # Save figure
  fig_save(o, g, fig_name)
}

# ---------------------------------------------------------
# Plot duration distributions 
# ---------------------------------------------------------
plot_durations = function(o, fig_name, model_input) {
  
  # Number of samples from each distribution
  n_samples = 10000
  
  # Colour map and palette to use
  colours = "viridis::viridis" 
  
  # Facet names can be long, wrap after n characters
  n_wrap = 30
  
  # Set a dictionary for facet names and order
  duration_dict = c(latency           = "Duration of latency phase",
                    presymptomatic    = "Duration of presymptomatic phase",
                    infectious_mild   = "Duration of non-severe disease",
                    infectious_severe = "Duration of severe disease",
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
  for (this_duration in names(model_input$durations)) {
    
    # Sample a ton of points 
    duration_samples = model_input$durations[[this_duration]](rep(0, n_samples))
    
    # Store as a datatable in a list
    plot_list[[this_duration]] = data.table(distribution = this_duration, 
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
    facet_wrap(~distribution, scales = "free", labeller = label_wrap_gen(n_wrap))
  
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
  fig_save(o, g, fig_name)
}

# ---------------------------------------------------------
# Viral load profile since day of infection
# ---------------------------------------------------------
plot_viral_load = function(o, fig_name, model_input) {
  
  # Number of days post infection to plot up until
  plot_days = 30
  
  # Descriptive names for post-infection durations
  period_names = c("Mean latency period", 
                   "Mean presymptomatic period", 
                   "Mean duration of non-severe disease", 
                   "Mean duration of severe disease")
  
  # Colour for viral load curve
  vl_colour = "grey40"
  
  # Number of time to sample durations
  n_samples = 100
  
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
  
  # Shorthard for duration sampling function list
  sample_fn = model_input$durations
  
  # Sample latency durations
  latency_samples = sample_fn$latency(rep(0, n_samples))
  
  # Mean latency from these samples
  mean_latency = round(mean(latency_samples), digits = n_dp)
  
  # Number of time points to offset duations by
  offset = list(best  = 10 ^ n_dp * mean_latency,
                lower = 10 ^ n_dp * min(latency_samples), 
                upper = 10 ^ n_dp * max(latency_samples))
  
  # Constuct viral load plotting dataframe
  plot_df = data.table(x  = time, 
                       y  = c(rep(0, offset$best),  head(vl, -offset$best)), 
                       y1 = c(rep(0, offset$lower), head(vl, -offset$lower)), 
                       y2 = c(rep(0, offset$upper), head(vl, -offset$upper))) %>%
    mutate(lower = pmin(y1, y2), upper = pmax(y1, y2))
  
  # Identify peaks to fill between
  vl_peaks = which(1 - plot_df$upper < 1e-6)
  plot_df$upper[min(vl_peaks) : max(vl_peaks)] = 1
  
  # Sample procedding durations
  mean_presymp = mean(sample_fn$presymptomatic(rep(0, n_samples)))
  mean_mild    = mean(sample_fn$infectious_mild(rep(0, n_samples)))
  
  # Vector of durations - bound below by zero and above by plot_days
  durations = c(mean_latency, mean_presymp, mean_mild)
  
  # Construct background duration plotting dataframe
  rect_df = data.table(period = factor(period_names, levels = period_names), 
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
  fig_save(o, g, fig_name)
}

# ---------------------------------------------------------
# Plot vaccine and acquired immunity profiles
# ---------------------------------------------------------
plot_immunity_profiles = function(o, fig_name, model_input) {
  
  # Which vaccines to plot, along with a description
  immunity_dict = c(vaccine  = "Immunity from COVID-19 vaccination", 
                    acquired = "Immunity from natural SARS-CoV-2 infection")
  
  # Colour scheme
  colours = c("navy", "goldenrod2")
  
  # Data point colour
  data_colour = "grey5"
  
  # Number of days to plot (from first dose / recovery from infection)
  n_days = 600
  
  # ---- Vaccine efficacy data ----
  
  # TODO: Add source and explanation...
  
  # Dose 1
  time_dose1 = c(1,  28)
  eff_dose1  = c(40, 20)
  
  # Dose 2
  time_dose2 = c(28, 168)
  eff_dose2  = c(85, 30)
  
  # Booster dose
  time_boost = c(168, 533)
  eff_boost  = c(90,  23)
  
  # ---- Format data frames ----
  
  # Total number of days we can plot for
  n_days_total = model_input$n_days - model_input$n_days_init
  
  # Extract effects for acquired and vaccine (and also booster)
  acquired_effect = model_input$acquired_immunity
  vaccine_effect  = model_input$vaccine$profile
  booster_effect  = model_input$booster_profile
  
  # Time indices for vaccine booster
  boost_idx1 = min(time_boost) : n_days_total
  boost_idx2 = 1 : (n_days_total - min(time_boost) + 1)
  
  # Apply this booster effect
  vaccine_effect[boost_idx1] = booster_effect[boost_idx2]
  
  # Overall efficacy datatable
  line_df = data.table(vaccine  = vaccine_effect, 
                       acquired = acquired_effect, 
                       day = 1 : n_days_total) %>%
    pivot_longer(cols = -day, names_to = "type") %>%
    arrange(type, day) %>%
    filter(day <= n_days) %>%
    mutate(value = value * 100, 
           type = factor(immunity_dict[type], levels = unname(immunity_dict))) %>%
    as.data.table()
  
  # ... of which is transmission blocking
  area_df = line_df %>%
    mutate(value = ifelse(type == immunity_dict[["vaccine"]], 
                          value * model_input$vaccine$transmission_blocking, value))
  
  # Dataframe for vaccine efficacy data points
  data_df = data.table(day   = c(time_dose1, time_dose2, time_boost), 
                       value = c(eff_dose1,  eff_dose2,  eff_boost), 
                       type  = immunity_dict[["vaccine"]])
  
  # ---- Produce plot ----
  
  # Plot the profiles over time along with maximum efficacies
  g = ggplot(area_df, aes(x = day, y = value, fill = type, colour = type)) + 
    geom_area(alpha = 0.5, size = 0) +
    geom_line(data = line_df, size = 2) + 
    geom_point(data = data_df, show.legend = FALSE, 
               colour = data_colour, size = 5) +
    facet_wrap(~type)
  
  # Use pre-defined colour scheme for distribution fills
  g = g + scale_fill_manual(values = colours) + 
    scale_colour_manual(values = colours)
  
  # Set axis labels
  g = g + xlab("Days since recovery from infection / first dose") + 
    ylab("Immunity from severe disease (%)")
  
  # Set x axes to user-defined limits (see options.R)
  g = g + expand_limits(y = c(0, 100)) + 
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) + 
    scale_y_continuous(expand = expansion(mult = c(0.00, 0.00)), 
                       breaks = seq(0, 100, by = 10))
  
  # Set facet labels
  g = facet_labels(g)
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(strip.text = element_text(size = 20),
          axis.title = element_text(size = 20),
          axis.text  = element_text(size = 14),
          axis.line  = element_blank(),
          panel.border = element_rect(size = 1, colour = "black", fill = NA),
          panel.spacing = unit(1, "lines"),
          legend.position = "none",
          strip.background = element_blank())
  
  # Save figure to file
  fig_save(o, g, fig_name)
}

# ---------------------------------------------------------
# Plot seasonality profile best fit and bounds
# ---------------------------------------------------------
plot_seasonality_profile = function(o, fig_name, labels = NULL, colours = NULL, ...) {
  
  # Extract model inputs
  model_inputs = list(...)
  
  # Number of models provided
  n_models = length(model_inputs)
  
  # Preallocate list for seasonality profiles
  season_list = list()
  
  # Loop through models
  for (i in 1 : n_models) {
    this_model = model_inputs[[i]]
    
    # Extract seasonality profile
    season_list[[i]] = data.table(
      date  = 1 : this_model$n_days,
      id    = this_model$.id, 
      name  = this_model$.name, 
      value = this_model$seasonality)
  }
  
  # Bind into single datatable
  season_df = rbindlist(season_list)
  
  # If colours are not provided, use default palette
  if (is.null(colours)) {
    colours = scales::hue_pal()(n_models)
    
  } else { # If colours are provided...
    
    # ... check appropriate number of them
    if (length(colours) != n_models)
      stop("Inconsistent number of colours provided")
  }
  
  # If labels not provided, use scenario names
  if (is.null(labels)) {
    labels = unique(season_df$name)
    
  } else { # If labels are provided...
    
    # ... check appropriate number of them
    if (length(labels) != n_models)
      stop("Inconsistent number of labels provided")
  }
  
  # Construct dictionary using default or custom labels
  label_dict = setNames(labels, unique(season_df$id))
  
  # Recode names if need be (could be a trivial step)
  plot_df = season_df %>%
    mutate(name = recode(id, !!!label_dict), 
           name = factor(name, label_dict), 
           metric = "Seasonal forcing on SARS-CoV-2 infectiousness per contact")
  
  # Minimum value (ie peak summer) for all curves
  threshold_df = plot_df %>%
    group_by(id, name, metric) %>%
    mutate(value = min(value)) %>%
    as.data.table()
  
  # Text to accompany minimum lines
  text_df = plot_df %>%
    group_by(id, name, metric) %>%
    slice_min(value, with_ties = FALSE) %>%
    mutate(text = paste("Peak summer\n effect:", value)) %>%
    as.data.table()
  
  # Plot all seasonality profiles provided
  g = ggplot(plot_df) + 
    aes(x = date, y = value, color = name) + 
    geom_line(size = 3) + 
    geom_hline(yintercept = 1, colour = "black", linetype = "dashed") + 
    facet_wrap(~metric)
  
  # Plot dashed curves to highlight peak summer
  g = g + geom_line(data = threshold_df, linetype = "dashed")
  
  # Text to accompany peak summer threshold lines
  g = g + geom_text(data = text_df, aes(label = text), 
                    hjust = 0.5, vjust = 1, nudge_y = -0.01)
  
  # Set colour scheme (could be trivial step)
  g = g + scale_colour_manual(values = colours)
  
  # Prettify axes
  g = g + expand_limits(y = c(0, 1)) + 
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
    scale_x_continuous(expand = expansion(mult = c(0, 0)))
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(text         = element_text(size = 14), 
          strip.text   = element_text(size = 24),
          axis.title   = element_blank(),
          axis.text    = element_text(size = 14),
          axis.line    = element_blank(),
          panel.border = element_rect(size = 1, colour = "black", fill = NA),
          panel.spacing = unit(1, "lines"),
          strip.background = element_blank(), 
          legend.text  = element_text(size = 16),
          legend.title = element_blank(),
          legend.key   = element_blank(),
          legend.position = "right", 
          legend.key.height = unit(2, "lines"),
          legend.key.width  = unit(2, "lines"),
          legend.box.background = element_rect())
  
  # Do not save figure if input is trivial  
  if (!is.null(fig_name))
    fig_save(o, g, fig_name)
  
  return(g)
}

# ---------------------------------------------------------
# Set heatmap colour scale
# ---------------------------------------------------------
apply_colourscale = function(g, f, idx, plot_df) {
  
  # Colour palette should be one value
  if (length(f$colour_palette) == 1) {
    colour_palette = f$colour_palette
    
  } else {  # ... or can also be n if using multi_scale
    colour_palette = f$colour_palette[idx]
  }
  
  # Colour limits should be 2-element vector
  if (is.vector(f$colour_limits)) 
    colour_limits = f$colour_limits
  
  # ... or can be nx2 matrix if using multi_scale
  if (is.matrix(f$colour_limits))
    colour_limits = f$colour_limits[idx, ]
  
  # If nothing is defined, take limits of data
  if (is.null(f$colour_limits))
    colour_limits = c(min(plot_df$value, na.rm = TRUE), 
                      max(plot_df$value, na.rm = TRUE))
  
  # By default, scale_fill_distiller use negative colour direction
  colour_direction = -1
  
  # However, this can be reversed using "-" as first char of palette name
  if (str_sub(colour_palette, 1, 1) == "-") {
    
    # Reverse colour, then remove the "-" character
    colour_direction = -colour_direction
    colour_palette   = str_sub(colour_palette, 2, -1)
  }
  
  # Default behaviour: continuous scale
  if (is.null(f$colour_bins)) {
    
    # Set continuous colour bar
    g = g + scale_fill_distiller(
      palette   = colour_palette, 
      direction = colour_direction, 
      limits    = colour_limits, 
      labels    = f$z_labels)
    
  } else {  # Alternative: discrete scale
    
    # Unlike scale_fill_distiller, scale_fill_fermenter doesn't like NAs in limits
    if (is.na(colour_limits[1])) colour_limits[1] = min(plot_df$value, na.rm = TRUE)
    if (is.na(colour_limits[2])) colour_limits[2] = max(plot_df$value, na.rm = TRUE)
    
    # Set continuous colour bar
    g = g + scale_fill_fermenter(
      palette   = colour_palette,
      direction = colour_direction,
      limits    = colour_limits,
      labels    = f$z_labels, 
      breaks    = f$colour_bins,
      guide     = guide_coloursteps(even.steps = FALSE)) 
  }
  
  return(g)
}

# ---------------------------------------------------------
# Determine which scenarios to plot based on user inputs
# ---------------------------------------------------------
fig_scenarios = function(o, args) {
  
  # Load default figure properties
  p = fig_defaults()
  
  # Check any inputs provided by the user are valid
  unknown_properties = setdiff(names(args), names(p))
  if (length(unknown_properties) > 0)
    stop("Unrecognised figure property: ", paste(unknown_properties, collapse = ", "))
  
  # Overwrite defaults for user-defined items
  p[names(args)] = args
  
  # Initiate figure list
  f = list(plot_baseline = p$plot_baseline)
  
  # ---- Figure scenarios ----
  
  # Check whether scenarios are defined as a matrix
  if (is.matrix(p$scenarios)) {
    
    # Interpreted as scenario groupings - store these groupings
    p$scenario_matrix = p$scenarios
    
    # Flatten scenarios into a vector
    p$scenarios = as.vector(p$scenarios)
  }
  
  # Check whether we actually want to plot a baseline
  if (f$plot_baseline == TRUE) {
    
    # Name of baseline (either default or an alternative)
    if (is.null(p$alt_baseline)) f$baseline_name = "baseline" else f$baseline_name = p$alt_baseline
    
    # If this baseline is already within the scenarios, do not plot it twice
    f$plot_baseline = !f$baseline_name %in% p$scenarios
    
  } else {
    
    # No baseline - use a pseudo
    f$baseline_name = p$scenarios[1]
    p$scenarios     = p$scenarios[-1]
    
    # Sanity check that a baseline is provided if plotting relative to that baseline
    if (p$relative == TRUE)
      stop("You must provide a baseline if 'relative' is set to TRUE")
  }
  
  # Append all non-baseline scenarios to scenario vector
  f$scenarios   = union(f$baseline_name, p$scenarios)
  f$n_scenarios = length(f$scenarios)
  
  # Full scenario names as defined in yaml file - try cheaper non-array call
  f$scenario_names = parse_yaml(o, "*read*", read_array = FALSE)[f$scenarios]
  
  # If these are array scenarios, we'll need the more expensive read_array version of *read*
  if (any(is.na(f$scenario_names)))
    f$scenario_names = parse_yaml(o, "*read*", read_array = TRUE)[f$scenarios]
  
  # Check if any names are still missing - throw an error if so
  missing_names = f$scenarios[is.na(f$scenario_names)]
  if (length(missing_names) > 0)
    stop("Scenario names not recognised: \n", paste(missing_names, collapse = "\n"))
  
  # Apply custom function to these scenario names if given
  if (!is.null(args$scenario_name_fn))
    f$scenario_names = args$scenario_name_fn(f$scenario_names)
  
  # ---- Relative to a baseline ----
  
  # Do we want to plot difference between scenarios and baseline
  if (p$relative == TRUE) {
    
    # Reduce number of scenarios by one as we won't plot baseline explictly
    f$n_scenarios = f$n_scenarios - 1
    
    # Sanity check that not only a baseline is provided
    if (f$n_scenarios == 0)
      stop("To use 'relative' you must provide some alternative scenarios and/or strategies")
    
    # We won't be plotting the baseline in this case
    f$plot_baseline = FALSE
  }
  
  # Throw an error if no scenarios identified
  if (f$n_scenarios == 0)
    stop("No scenarios identified - plot a baseline and/or scenarios and/or strategies")
  
  # ---- Scenario groupings ----
  
  # Check whether scenarios should be grouped
  if (!is.null(p$scenario_matrix)) {
    
    # Construct function to append baseline details (if needed)
    add_baseline_fn = function(x) {
      
      # Normal use case: a baseline for each scenario_group
      if (!p$aes_reverse)
        y = expand_grid(scenario_group = unique(x$scenario_group),
                        scenario_type  = x$scenario_type[1],
                        scenario_id    = f$baseline_name)
      
      # Alternative use case: baseline in it's own facet
      if (p$aes_reverse)
        y = expand_grid(scenario_group = "Baseline",
                        scenario_type  = "Baseline",
                        scenario_id    = f$baseline_name)
      
      # Bind with existing dataframe
      z = rbind(y, x)
    }
    
    # Generate a dataframe describing the groupings
    f$scenario_groups = as.data.table(p$scenario_matrix) %>%
      mutate(scenario_group = rownames(p$scenario_matrix)) %>%
      pivot_longer(cols = -scenario_group, 
                   names_to  = "scenario_type", 
                   values_to = "scenario_id") %>%
      {if (f$plot_baseline) add_baseline_fn(.) else .} %>%
      left_join(data.table(scenario_id = f$scenarios, 
                           scenario    = f$scenario_names), 
                by = "scenario_id") %>%
      select(scenario, scenario_type, scenario_group) %>%
      as.data.table()
    
    # Store the number of groups for easy referencing
    f$n_scenario_groups = c(nrow(p$scenario_matrix), 
                            ncol(p$scenario_matrix))
  }
  
  # Do we have multiple scenarios to plot?
  if (f$n_scenarios > 1) {
    
    # They can either be grouped...
    if (!is.null(p$scenario_matrix)) {
      f$plot_type = "scenario_group"
      
    } else {  # ... or not
      f$plot_type = "scenario"
    }
  }
  
  return(f)
}

# ---------------------------------------------------------
# Generate figure properties based on user inputs
# ---------------------------------------------------------
fig_properties = function(o, f, yaml, args) {
  
  # Load defaults and overwrite user-defined items
  p = fig_defaults()
  p[names(args)] = args
  
  # ---- Metrics and groupings to plot ----
  
  # Plot all possible metrics if plot_metrics arguments not given
  if (is.null(p$plot_metrics)) p$metrics = yaml$metrics$df$metric
  else p$metrics = p$plot_metrics
  
  # Metrics suitable for cumulative, aggregation, and coverages plots 
  metrics_cum = yaml$metrics$df[cumulative == TRUE, metric]
  metrics_agg = yaml$metrics$df[aggregate == TRUE,  metric]
  metrics_cov = yaml$metrics$df[coverage == TRUE,   metric]
  
  # Two cases where we'll be plotting cumulatively
  # 
  # Case 1: User specifically asks for cumulative metrics
  # Case 2: Bar plots that sum metrics over time
  plot_cumulative = p$cumulative == TRUE ||
    p$plot_geom == "bar" && p$summarise == "sum"
  
  # If plotting cumulatively, filter out unsuitable metrics
  if (plot_cumulative == TRUE) {
    p$metrics = intersect(p$metrics, metrics_cum)
    
    # Also ensure flag is set, important for post processing
    p$cumulative = TRUE
  }
  
  # Some metrics can only be plotted as distributions
  metrics_dist = yaml$metrics$df[temporal == FALSE, metric]
  
  # If plotting histograms, filter out unsuitable metrics
  if (p$plot_geom == "hist")
    p$metrics = intersect(p$metrics, metrics_dist)
  
  # If no grouping defined we'll be ploting by multiple metrics
  if (is.null(p$plot_by)) {
    
    # For some metrics it is nonsensical to aggregate
    p$metrics = intersect(p$metrics, c(metrics_agg, metrics_cov))
    
    # Colours will be by metric
    if (f$n_scenarios == 1)
      f$plot_type = "metric"
    
  } else {
    
    # Throw an error if trying to use scenario groupings too
    if (identical(f$plot_type, "scenario_group"))
      stop("You cannot use 'plot_by' with scenario grouping: ", 
           "provide scenarios as a vector to use 'plot_by'")
    
    # We'll be plotting by the grouping defined
    f$plot_type = "group"
    
    # Subset metrics that have been requested AND are grouped by plot_by
    grouping_idx     = yaml$metrics$groupings == p$plot_by
    grouping_metrics = names(yaml$metrics$groupings[grouping_idx])
    
    # Remove any metrics that do not have this grouping
    p$metrics = intersect(p$metrics, grouping_metrics)
    
    # Do not plot coverages if stacking plots
    if (p$plot_geom != "line")
      p$metrics = setdiff(p$metrics, metrics_cov)
  }
  
  # Warn the user if there are no metrics left to plot
  if (length(p$metrics) == 0) {
    warning("No plotting metrics identified - skipping plot", 
            "\n  - Plot geom: ", p$plot_geom, 
            "\n  - Plot type: ", f$plot_type, 
            "\n  - Plot by: ",   ifelse(!is.null(p$plot_by), p$plot_by, "none"))
    
    # Return out early
    return(list.append(p, n_metrics = 0))
  }
  
  # Number of metrics to plot
  f$n_metrics = length(p$metrics)
  
  # ---- Metric descriptions ----
  
  # Metric descriptions can be different if plotting cumulative results
  if (plot_cumulative) metric_dict = yaml$dict$cumulative 
  else metric_dict = yaml$dict$metric
  
  # Apply and subset the dictionary
  p$metric_names = metric_dict[p$metrics]
  
  # ---- Group stacking ----
  
  # Tornado plots are a special case of impact bars
  if (p$group_tornado == TRUE) {
    
    # Throw an error if not plotting relative to some baseline
    if (p$relative == FALSE)
      stop("Tornado plots should be used for sensitivity analyses: set 'relative = TRUE'")
    
    # Ideally a matrix of scenarios is defined, although this shouldn't be an error
    if (f$plot_type != "scenario_group")
      warning("Tornado plots are most effective when used with a matrix of scenarios")
    
    # Ensure groups are being stacked
    p$group_stack = TRUE
  }
  
  # ---- Colours ----
  
  # Skip colour generation if override_colours is NA
  if (!is.null(p$override_colours) && is.na(p$override_colours)) p$colours = NA
  else {
    
    # Scenario colours: one per scenario, plus baseline if required
    if (f$plot_type == "scenario" || p$plot_geom == "hist") {
      
      # One colour per scenario
      n_colours = f$n_scenarios 
      
      # Create colour vector from palette
      p$colours = colour_scheme(o$palette_scenario, n = n_colours)
      
      # Append baseline colour if required
      if (f$plot_baseline == TRUE) 
        p$colours = c(o$baseline_colour, p$colours[-n_colours])
    }
    
    # Scenario group colours: repeat colours for each scenario group
    if (f$plot_type == "scenario_group") {
      
      # Default or reverse aesthetic ordering (colours-linetype or linetype-colours)
      if (!p$aes_reverse) n_colours = f$n_scenario_groups[1]  # One colour per group
      if (p$aes_reverse)  n_colours = f$n_scenario_groups[2]  # One colour per type
      
      # Create colour vector from palette
      p$colours = colour_scheme(o$palette_scenario, n = n_colours)
      
      # Append baseline colour if required
      if (f$plot_baseline == TRUE) {
        p$colours = c(o$baseline_colour, p$colours)
        
        # Increment number of colours required
        n_colours = n_colours + 1
      }
    }
    
    # Metric colours: a colour for each possible metric, indexed if not plotting everything
    if (f$plot_type == "metric") {
      n_colours = length(p$metrics)
      
      # Colours for all metrics we may wish to plot
      all_metric_colours = colour_scheme(o$palette_metric, n = length(yaml$dict$metric))	
      
      # Subset this for metrics we are plotting
      p$colours = all_metric_colours[names(yaml$dict$metric) %in% p$metrics]
    }
    
    # Colours by group depends on number in grouping
    if (f$plot_type == "group" && p$plot_geom != "hist") {
      
      # Number in this grouping
      n_colours = length(yaml$count[[p$plot_by]])
      
      # Grouping by age is a special case
      if (p$plot_by == "age")
        n_colours = length(o$plot_ages)
      
      # User defined palette for this grouping
      group_palette = o[[paste0("palette_", p$plot_by)]]
      
      # Generate set of colours
      p$colours = colour_scheme(group_palette, n = n_colours)
    }
    
    # Force override colours if desired
    if (!is.null(p$override_colours)) {
      
      # Throw an error if override does not have correct number of values
      if (length(p$override_colours) != n_colours)
        stop("Inconsistent number of manual colours provided ", 
             "(", n_colours, " needed, ", length(p$override_colours), " provided)")
      
      # Apply the override
      p$colours[!is.na(p$override_colours)] = na.omit(p$override_colours)
    }
  }
  
  # Special case for fill colours: repeat when grouping scenarios
  if (f$plot_type == "scenario_group") {
    n_rep = f$n_scenario_groups[2]
    
    # If no baseline, simply repeat n_rep times
    if (f$plot_baseline == FALSE)
      p$colours = rep(p$colours, times = n_rep)
    
    # If a baseline, all except for this should be repeated
    if (f$plot_baseline == TRUE)
      p$colours = c(p$colours[1], rep(p$colours[-1], times = n_rep))
  }
  
  # ---- Facets ----
  
  # If only one wrap value given, x and y strips use the same value
  if (length(p$n_wrap) == 1)
    p$n_wrap = c(p$n_wrap, p$n_wrap)
  
  # Shorthand for facet_wrap strip wrapping
  p$label_wrap = label_wrap_gen(p$n_wrap[1])
  
  # Shorthand for facet_grid strip wrapping
  #
  # NOTE: This can handle different values for x and y strips
  p$label_grid = labeller(.cols = label_wrap_gen(p$n_wrap[1]), 
                          .rows = label_wrap_gen(p$n_wrap[2]))
  
  # If facet_custom is defined, ensure it is a string
  if (!is.null(f$facet_custom) && !is.character(f$facet_custom))
    stop("Argument 'facet_custom' must be an evaluable string")
  
  # ---- Plotting days or dates ----
  
  # Plotting dates
  p$plot_from = max(as.numeric(p$plot_from), 1)
  p$plot_to   = min(as.numeric(p$plot_to), yaml$n_days)
  
  # Plot day numbers not dates by default
  p$plot_dates = FALSE
  
  # If date_from is given, turn flag on
  if (!is.null(p$date_from)) {
    p$plot_dates = TRUE
    
    # Ensure yyyy-mm-dd format
    p$date_from = format_date(p$date_from) - p$plot_from
  }
  
  # ---- Fontsizes ----
  
  # Heatmaps have title-subtitle and axes labels
  if (p$plot_geom == "tile") {
    p$fontsize = c(40, 24, 24, 16, 16, 16, 14)
    
  } else {  # All other plots don't...
    
    # Font size (title, facets, legend, ticks) depends on how many panels we have
    if (length(p$metrics) >= 10) p$fontsize = c(32, 12, 11, 11)
    if (length(p$metrics) <= 9)  p$fontsize = c(32, 17, 14, 14) 
    if (length(p$metrics) <= 4)  p$fontsize = c(32, 20, 15, 15) 
  }
  
  # Force override font sizes if desired
  if (!is.null(p$override_fontsize)) {
    
    # Heatmaps have 6 text size values
    if (p$plot_geom == "tile") {
      err_msg = "title, subtitle, axes, x-strip, y-strip, legend, ticks"
      
    } else {  # All other plots have 4
      err_msg = "title, facets, legend, ticks"
    }
    
    # Throw an error if override does not have correct number of values
    if (length(p$override_fontsize) != length(p$fontsize))
      stop("Input 'override_fontsize' must have ", length(p$fontsize), " values (", err_msg, ")")
    
    # Apply the override
    p$fontsize[!is.na(p$override_fontsize)] = na.omit(p$override_fontsize)
  }
  
  # ---- Append all details ----
  
  # Combine lists (any scenario info already in f takes precedence)
  f = utils::modifyList(p, f)
  
  # Order list items alphabetically
  f = f[sort(names(f))]
  
  return(f)
}

# ---------------------------------------------------------
# Default values for figure properties
# ---------------------------------------------------------
fig_defaults = function() {
  
  # Define default values
  defaults = list(
    plot_geom     = "line",      # GG plot geom
    summarise     = "sum",       # Summarise data function for bar and tile plots
    plot_baseline = TRUE,        # Flag for plotting a baseline
    alt_baseline  = NULL,        # Define some alternative 'baseline'
    scenarios     = NULL,        # Alternative scenario names (a vector or matrix)
    relative      = FALSE,       # Outcomes plotted relative to baseline
    aes_reverse   = FALSE,       # Reverse default grouped scenarios aesthetics (colour-linetype) 
    aes_linetype  = TRUE,        # Whether to differentiate grouped scenarios by dashed lines
    plot_metrics  = NULL,        # Vector of model metrics to plot (all possible metrics by default)
    plot_by       = NULL,        # Grouping to plot by (eg age or variant)
    group_stack   = TRUE,        # Impact bars stacked by group (alternative is side by side)
    group_dodge   = FALSE,       # Impact bars unstacked, but grouped by scenario on single x axis (overrules group_stack)
    group_tornado = FALSE,       # Impact bars as a tornado plot (use with a matrix of scenarios and 'relative')
    group_idx     = 1,           # For plots that can only show one group, which group
    person_days   = 1e5,         # Metrics per `person_days` per day
    cumulative    = FALSE,       # Plot cumulative outcomes
    uncertainty   = TRUE,        # Plot uncertainty bounds (error bars for geom_bar)
    plot_from     = 1,           # First time point for plotting
    date_from     = NULL,        # Convert days to dates starting from date_from
    date_labels   = "%b %y",     # Date label format (month-year: "%b %y", day-month: "%d %b")
    date_breaks   = "1 month",   # Distance between date labels on x axis for temporal plots
    plot_to       = Inf,         # Cut plotting after so many days
    y_min         = 1,           # Set a minimum for y_max for all metrics
    y_expand      = 0.05,        # Expand y axis as a multiple of largest data point
    bar_width     = 0.9,         # Bar width of impact bars 
    bar_legend    = FALSE,       # Use a legend (rather than tick labels) to define scenarios in bar plots
    legend_rows   = 2,           # Number of legend rows
    legend_by_row = TRUE,        # Order legend entries by row (rather than by column)
    facet_rows    = NULL,        # Number of facet rows for facet_wrap plots
    facet_custom  = NULL,        # Define custom facets as evaluable string, eg using ggh4x::facet_grid2
    facet_by      = "scenario",  # Facet by either 'scenario' or 'group' for number of infection histograms
    facet_labels  = TRUE,        # Flag for whether to tag/label each facet with a capital letter
    n_wrap        = 20,          # Number of chars per line per strip (grid_plots can take 2nd value)
    x_wrap        = 30,          # Number of chars per line per x axis label in relevant bar plots
    x_rotate      = 50,          # Degrees to rotate x axis labels in relevant bar plots
    line_width    = 2,           # Line thickness for geom_line plots
    override_colours  = NULL,    # Vector of custom colours
    override_fontsize = NULL,    # Vector of custom font sizes 
    scenario_name_fn  = NULL)    # Custom scenario naming function
  
  # Define default values specfic for heatmaps
  defaults_heat = list(
    dif_array      = NULL,        # Optional input for calculating difference between two arrays
    dif_relative   = FALSE,       # Should the difference be caclulated relative to primary
    n_interpolate  = NULL,        # Interpolate to nxn points in each facet
    multi_scale    = NULL,        # Produce a new colour scale for each value along one dimension ('v' or 'w')
    colour_palette = "Spectral",  # Colour scale (use "-xxx" to reverse direction)*
    colour_limits  = NULL,        # Data and colour bounds (2 value vector)**
    colour_bins    = NULL,        # Discretise colour bar into blocks (can be uneven)
    contour_val    = NA,          # Plot a contour at a given value
    contour_tol    = 0.05,        # Tolerance around contour_val for point selection
    contour_fn     = "y ~ x",     # Fit contour through points within contour_tol of contour_val
    contour_by     = NA,          # Different linetype for all contours across one dimension ('v' or 'w')
    x_lab          = NULL,        # X-axis label
    y_lab          = NULL,        # Y-axis label
    x_scale_fn     = waiver(),    # X-axis tick labelling/scaling function
    y_scale_fn     = waiver(),    # Y-axis tick labelling/scaling function
    z_labels       = waiver(),    # Colourbar tick labelling/scaling function
    w_percent      = FALSE,       # Convert 3rd dim values to percentage
    v_percent      = FALSE,       # Convert 4th dim values to percentage
    plot_title     = NULL,        # Override plot title
    plot_subtitle  = NULL)        # Override plot subtitle
  
  # NOTES:
  #   * Use 1 or n values for multi_scale plots
  #  ** Use a 2-element vector or nx2 matrix for multi_scale plots
  
  # Concatenate all defaults
  all_defaults = c(defaults, defaults_heat)
  
  return(all_defaults)
}

# ---------------------------------------------------------
# Save a ggplot figure to file with default settings
# ---------------------------------------------------------
fig_save = function(o, g, ..., path = "figures", width = o$save_width, height = o$save_height) {
  
  # Collapse inputs into vector of strings
  fig_name_parts = unlist(list(...))
  
  # Construct file name to concatenate with file path
  save_name = paste(fig_name_parts, collapse = " - ")
  
  # Repeat the saving process for each image format in figure_format
  for (fig_format in o$figure_format) {
    save_pth  = paste0(o$pth[[path]], save_name, ".", fig_format)
    
    # Save figure (size specified in options.R)
    ggsave(save_pth, 
           plot   = g, 
           device = fig_format, 
           dpi    = o$save_resolution, 
           width  = width, 
           height = height, 
           units  = o$save_units)
  }
}

