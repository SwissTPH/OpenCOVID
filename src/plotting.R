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
  if (length(f$metrics) == 0)
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
  
  # Normal use case is to reverse scenarios (so they are on top)
  if (f$plot_type == "group") scenario_levels = f$scenario_names
  else scenario_levels = rev(f$scenario_names)
  
  # Concatenate plotting dataframes for all scenarios
  plot_df = rbindlist(plot_list) %>%
    mutate(metric   = recode(metric,   !!!f$metric_names), 
           scenario = recode(scenario, !!!f$scenario_names), 
           scenario = factor(scenario, levels = scenario_levels))
  
  # ---- Apply dates if necessary ----
  
  # Check whether dates are preferred to day numbers
  if (f$plot_dates == TRUE) {
    
    # Vector of all dates simulated
    all_dates = f$date_from + (1 : baseline$input$n_days) - 1 
    
    # Apply these dates
    plot_df$date = all_dates[plot_df$date]
  }
  
  # ---- Create figure ----
  
  # Initiate ggplot structure with each metric in it's facet
  g = ggplot(plot_df, aes_string(x = "date", y = "value", colour = f$plot_type, fill = f$plot_type))
  
  # Several options for plotting geom - default is a line graph
  if (f$plot_geom == "line") {
    
    # Plot best estimate line with uncertainty bounds
    g = g + geom_ribbon(aes(ymin = lower, ymax = upper), linetype = 0, alpha = 0.3) + 
      geom_line(size = 2)
    
    # Area plot over time - useful when plotting age groups
  } else if (f$plot_geom == "area") {
    g = g + geom_area()
    
  } else {  # Throw an error if any other value is provided
    stop("Value of plot_geom '", f$plot_geom, "' not recognised")
  }
  
  # If plotting by group, use facet grid (could be busy!)
  if (f$plot_type == "group") {
    g = g + facet_grid(metric~scenario, scales = "free", 
                       labeller = label_wrap_gen(f$n_wrap))
    
  } else {  # Otherwise just wrap by metric
    g = g + facet_wrap(~metric, scales = "free", nrow = f$facet_rows, 
                       labeller = label_wrap_gen(f$n_wrap))
  }
  
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
    apply_guide = guide_legend(reverse = TRUE, nrow = f$legend_rows, byrow = TRUE)
    g = g + guides(colour = apply_guide, fill = apply_guide)
  }
  
  # Prettify y-axis
  g = g + expand_limits(y = c(0, f$expand_limit)) +  # Restrain upper y limit from being really low
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.05)))
  
  # Prettify x-axis: day numbers
  if (f$plot_dates == FALSE)
    g = g + scale_x_continuous(labels = comma, expand = expansion(mult = c(0, 0)))
  
  # Prettify x-axis: dates
  if (f$plot_dates == TRUE)
    g = g + scale_x_date(date_breaks = "1 month",  # "2 weeks"
                         date_labels = "%b %y",  # Month-year: "%b %y"; Day-month: "%d %b"
                         expand = expansion(mult = c(0, 0)))
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(plot.title   = element_text(size = f$fontsize[1], hjust = 0.5),
          strip.text   = element_text(size = f$fontsize[2]),
          axis.text    = element_text(size = f$fontsize[4]),
          axis.title   = element_blank(),
          legend.text  = element_text(size = f$fontsize[3]),
          legend.title = element_blank(),
          legend.key   = element_blank(),
          legend.position = "bottom", 
          legend.key.height = unit(2, "lines"),
          legend.box.background = element_rect())
  
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
plot_impact = function(o, fig_name, error_bars = TRUE, ...) {
  
  # Scenarios to plot
  f = fig_scenarios(o, list(...))
  
  # Load baseline (or alt baseline) file
  baseline = try_load(o$pth$scenarios, f$baseline_name)
  
  # ---- Figure properties ----
  
  # Collate key word arguments and set cumulative to true by default
  args = list_modify(list(...), plot_geom = "bar", cumulative = TRUE)
  
  # Collate and interpret function inputs so we know what to plot
  f = fig_properties(o, f, baseline$input, args)
  
  # If nothing to plot, return out
  if (length(f$metrics) == 0)
    return()
  
  # A bit more work would be needed for this
  if (f$plot_type == "group")
    stop("Impact bar plots have not been properly tested for plotting groupings yet")
  
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
      scale_x_discrete(labels = f$scenarios) + 
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
          axis.text.y = element_text(size = f$fontsize[4]),
          axis.text.x = element_text(size = f$fontsize[4], hjust = 1, angle = 50),
          axis.title  = element_blank())
  
  # Save these figures to file
  fig_save(o, g, fig_name)
  
  return(g)
} 

# ---------------------------------------------------------
# Plot heat maps for multi-dimensional array scenarios
# ---------------------------------------------------------
plot_heatmap = function(o, fig_name, array, plot_df = NULL, ...) {
  
  # ---- Figure properties ----
  
  # Collate key word arguments and set scenarios
  args = list_modify(list(...),
                     plot_geom = "tile",
                     override_colours = NA)
  
  # Scenarios to plot
  f = fig_scenarios(o, args)
  
  # Load baseline (or alt baseline) file
  baseline = try_load(o$pth$scenarios, f$baseline_name)
  
  # Collate and interpret function inputs so we know what to plot
  f = fig_properties(o, f, baseline$input, args)
  
  # If nothing to plot, return out
  if (length(f$metrics) == 0)
    return()
  
  # Multiple metrics are OK if plotting 2D arrays only
  if (length(f$metrics) > 1)
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
  if (n_dims == 3) g = g + facet_wrap(~w,  labeller = label_wrap_gen(f$n_wrap))
  if (n_dims == 4) g = g + facet_grid(w~v, labeller = label_wrap_gen(f$n_wrap))
  
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
    if (exists("dim_names")) g = g + ylab(dim_names[1])
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
# Plot normalised disease state of population over time
# ---------------------------------------------------------
plot_disease_state = function(o, p, states) {
  
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
  fig_save(o, g1, "Population disease and care states - stacked area")
  fig_save(o, g2, "Population disease and care states - line")
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
plot_durations = function(o, model_input) {
  
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
  fig_save(o, g, "Duration distributions")
}

# ---------------------------------------------------------
# Viral load profile since day of infection
# ---------------------------------------------------------
plot_viral_load = function(o, model_input) {
  
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
  fig_save(o, g, "Viral load profile")
}

# ---------------------------------------------------------
# Plot vaccine and acquired immunity profiles
# ---------------------------------------------------------
plot_immunity_profiles = function(o, model_input) {
  
  # Which vaccines to plot, along with a description
  immunity_dict = c(vaccine  = "Immunity from vaccination", 
                    acquired = "Immunity from natural infection")
  
  # Colour scheme
  colours = c("navy", "goldenrod2")
  
  # Number of days to plot (from first dose / recovery from infection)
  n_days = 600
  
  # ---- Format data frames ----
  
  # Total number of days we can plot for
  n_days_total = model_input$n_days - model_input$n_days_init
  
  # Overall efficacy datatable
  line_df = data.table(vaccine  = model_input$vaccine$profile, 
                       acquired = model_input$acquired_immunity, 
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
  
  # ---- Produce plot ----
  
  # Plot the profiles over time along with maximum efficacies
  g = ggplot(area_df, aes(x = day, y = value, fill = type, colour = type)) + 
    geom_area(alpha = 0.5, size = 0) +
    geom_line(data = line_df, size = 3) + 
    facet_wrap(~type)
  
  # Use pre-defined colour scheme for distribution fills
  g = g + scale_fill_manual(values = colours) + 
    scale_colour_manual(values = colours)
  
  # Set axis labels
  g = g + xlab("Days since first dose / recovery from infection") + 
    ylab("Immunity from severe disease (%)")
  
  # Set x axes to user-defined limits (see options.R)
  g = g + expand_limits(y = c(0, 100)) + 
    scale_y_continuous(expand = expansion(mult = c(0, 0)), breaks = seq(0, 100, by = 10)) + 
    scale_x_continuous(expand = expansion(mult = c(0, 0))) #, breaks = seq(0, n_days, by = 5))
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(strip.text   = element_text(size = 20),
          axis.title   = element_text(size = 26),
          axis.text    = element_text(size = 14),
          legend.text  = element_text(size = 18),
          legend.title = element_blank(),
          legend.key   = element_blank(),
          legend.position = "bottom", 
          legend.key.height = unit(2, "lines"),
          legend.box.background = element_rect())
  
  # Save figure to file
  fig_save(o, g, "Immunity profiles")
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
  
  # ---- Figure scenarios ----
  
  # Initiate figure list
  f = list(plot_baseline = p$plot_baseline)
  
  # Check whether we actually want to plot a baseline
  if (f$plot_baseline == TRUE) {
    
    # Name of baseline (either default or an alternative)
    if (is.null(p$alt_baseline)) f$baseline_name = "baseline" else f$baseline_name = p$alt_baseline
    
  } else {
    
    # No baseline - use a pseudo
    f$baseline_name = p$scenarios[1]
    p$scenarios     = p$scenarios[-1]
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
  
  # Do we want to plot difference between scenarios and baseline
  if (p$relative == TRUE) {
    
    # Sanity check that a baseline is provided
    if (f$plot_baseline == FALSE)
      stop("You must provide a baseline if 'relative' is set to TRUE")
    
    # Reduce number of scenarios by one as we won't plot baseline explictly
    f$n_scenarios = f$n_scenarios - 1
    
    # Sanity check that not only a baseline is provided
    if (f$n_scenarios == 0)
      stop("To use 'relative' you must provide some alternative scenarios and/or strategies")
  }
  
  # Store flag so outcomes can be relatively calculated if needed
  f$relative = p$relative
  
  # Throw an error if no scenarios identified
  if (f$n_scenarios == 0)
    stop("No scenarios identified - plot a baseline and/or scenarios and/or strategies")
  
  # Do we have multiple scenarios to plot?
  if (f$n_scenarios > 1) {
    
    # If so, we can't also have multiple groups on the same plot
    # if (!is.null(p$plot_by))
    #   stop("You cannot plot groupings for multiple scenarios at once")
    
    # We'll be plotting by multiple scenarios
    f$plot_type = "scenario"
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
  
  # Not all metrics can be plotted cumulatively
  if (p$cumulative == TRUE) {
    
    # Subset metrics to only those suitable for cumulative plotting
    cumulative_df = filter(yaml$metrics$df, cumulative == TRUE)
    p$metrics     = intersect(p$metrics, cumulative_df$metric)
  }
  
  # If no grouping defined we'll be ploting by multiple metrics
  if (is.null(p$plot_by)) {
    
    # For some metrics it is nonsensical to aggregate
    aggregate_df = filter(yaml$metrics$df, aggregate == TRUE)
    p$metrics    = intersect(p$metrics, aggregate_df$metric)
    
    # Colours will be by metric
    if (f$n_scenarios == 1)
      f$plot_type = "metric"
    
  } else {
    
    # We'll be plotting by the grouping defined
    f$plot_type = "group"
    
    # Subset metrics that have been requested AND are grouped by plot_by
    grouping_idx     = yaml$metrics$groupings == p$plot_by
    grouping_metrics = names(yaml$metrics$groupings[grouping_idx])
    
    # Remove any metrics that do not have this grouping
    p$metrics = intersect(p$metrics, grouping_metrics)
  }
  
  # Warn the user if there are no metrics left to plot
  if (length(p$metrics) == 0) {
    warning("No plotting metrics identified - skipping plot", 
            "\n  - Plot geom: ", p$plot_geom, 
            "\n  - Plot type: ", f$plot_type, 
            "\n  - Plot by: ",   ifelse(!is.null(p$plot_by), p$plot_by, "none"))
    
    # Return out early
    return(p)
  }
  
  # ---- Metric descriptions ----
  
  # Metric descriptions can be different if plotting cumulative results
  if (p$cumulative == TRUE) metric_dict = yaml$dict$cumulative 
  else metric_dict = yaml$dict$metric
  
  # Apply and subset the dictionary
  p$metric_names = metric_dict[p$metrics]
  
  # ---- Colours ----
  
  # Skip colour generation if override_colours is NA
  if (!is.null(p$override_colours) && is.na(p$override_colours)) p$colours = NA
  else {
    
    # Scenario colours: redefine on each function call
    if (f$plot_type == "scenario") {
      
      # Create colour vector from palette, considering baseline if required
      if (p$plot_baseline == FALSE) p$colours = colour_scheme(o$palette_scenario, n = f$n_scenarios)
      else p$colours = c(o$baseline_colour, colour_scheme(o$palette_scenario, n = f$n_scenarios - 1))
    }
    
    # Metric colours: a colour for each possible metric, indexed if not plotting everything
    if (f$plot_type == "metric") {
      
      # Colours for all metrics we may wish to plot
      all_metric_colours = colour_scheme(o$palette_metric, n = length(yaml$dict$metric))	
      
      # Subset this for metrics we are plotting
      p$colours = all_metric_colours[names(yaml$dict$metric) %in% p$metrics]
    }
    
    # Colours by group depends on number in grouping
    if (f$plot_type == "group") {
      
      # Number in this grouping
      n_group = length(yaml$count[[p$plot_by]])
      
      # Grouping by age is a special case
      if (p$plot_by == "age")
        n_group = length(o$plot_ages)
      
      # User defined palette for this grouping
      group_palette = o[[paste0("palette_", p$plot_by)]]
      
      # Generate set of colours
      p$colours = colour_scheme(group_palette, n = n_group)
    }
    
    # Force override colours if desired
    if (!is.null(p$override_colours)) {
      
      # Throw an error if override does not have correct number of values
      if (length(p$override_colours) != length(p$colours))
        stop("Inconsistent number of manual colours provided ", 
             "(", length(p$colours), " needed, ", length(p$override_colours), " provided)")
      
      # Apply the override
      p$colours[!is.na(p$override_colours)] = na.omit(p$override_colours)
    }
  }
  
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
    plot_geom     = "line",    # GG plot geom
    plot_baseline = TRUE,      # Flag for plotting a baseline
    alt_baseline  = NULL,      # Define some alternative 'baseline'
    scenarios     = NULL,      # Vector of alternative scenario names
    relative      = FALSE,     # Outcomes plotted relative to baseline
    plot_metrics  = NULL,      # Vector of model metrics to plot (all possible metrics by default)
    plot_by       = NULL,      # Grouping to plot by (eg age or variant)
    group_idx     = 1,         # For plots that can only show one group, which group
    person_days   = 1e5,       # Metrics per `person_days` per day
    cumulative    = FALSE,     # Plot cumulative outcomes
    plot_from     = 1,         # First time point for plotting
    date_from     = NULL,      # Convert days to dates starting from date_from
    plot_to       = Inf,       # Cut plotting after so many days
    expand_limit  = 1,         # Set a minimum for y_max for all metrics
    legend_rows   = 2,         # Number of legend rows
    facet_rows    = NULL,      # Number of facet rows for facet_wrap plots
    n_wrap        = 20,        # Number of chars per strip line in grid plots
    override_colours  = NULL,  # Vector of custom colours
    override_fontsize = NULL,  # Vector of custom font sizes 
    scenario_name_fn  = NULL)
  
  # Define default values specfic for heatmaps
  defaults_heat = list(
    summarise      = "sum",       # Summarise function to squash temporal to constant
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
    
    # Save figure (size specified in results.R)
    ggsave(save_pth, plot = g, width = width, height = height)
  }
}

