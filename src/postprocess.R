###########################################################
# POST PROCESS
#
# Process model output ready for plotting.
#
###########################################################

# ---------------------------------------------------------
# Calculate effective reproduction number
# ---------------------------------------------------------
calculate_Re = function(incidence, si, Re_window) {
  
  # See: https://cran.r-project.org/web/packages/EpiEstim/vignettes/demo.html
  
  # Number of daily observations
  n_days = length(incidence)
  
  # Use Re_window here for non-default time window
  Re_win  = Re_window - 1
  Re_days = seq(2, n_days - Re_win) # Starting at 2 as conditional on the past observations
  
  # Check if serial interval has already been summarised
  if (is.list(si)) si_dist = si[c("mean", "std")]
  
  # If still a vector, calculate normal distribution parameters
  if (is.numeric(si)) si_dist = list(mean = mean(si), std = sd(si))
  
  # We'll use the mean and standard deviation of this data
  si_list = list(mean_si = si_dist$mean, 
                 std_si  = si_dist$std, 
                 t_start = Re_days,
                 t_end   = Re_days + Re_win)
  
  # Use estimate_R function from EpiEstim package
  Re_list = estimate_R(incid  = incidence, 
                       method = "parametric_si",
                       config = make_config(si_list))
  
  # Extract Re estimate along with 95% CI bounds
  Re_df = Re_list$R %>%
    select(date  = t_start,
           value = "Mean(R)", 
           lower = "Quantile.0.025(R)", 
           upper = "Quantile.0.975(R)") %>%
    filter(date >= Re_window) %>%
    right_join(data.frame(date = 1 : n_days), 
               by = "date") %>%
    arrange(date) %>%
    mutate(incidence = incidence)
  
  # Append serial distribution 
  Re_info = list(Re_df = Re_df, si = si_dist)
  
  return(Re_info)
}

# ---------------------------------------------------------
# Aggregate groupings for appropriate metrics
# ---------------------------------------------------------
aggregate_results = function(input, raw_output) {
  
  # NOTE: This is the first step when processing raw results.
  #       It is in a seperate function because we may want to
  #       only aggregate, and not summarise raw results.
  
  # Metrics that can be aggregated
  agg_metrics = input$metrics$df[aggregate == TRUE, metric]
  
  # Need to aggregate grouped metrics: get first grouping per metric to avoid double counting
  first_grouping = raw_output %>% 
    filter(metric %in% agg_metrics) %>% 
    select(metric, grouping) %>%
    unique() %>%
    filter(!grouping %in% c("none", "na")) %>%
    group_by(metric) %>%
    slice(1) %>%
    ungroup() %>%
    setDT()
  
  # Aggregate values of first grouping
  agg_df = first_grouping %>%
    left_join(raw_output,
              by = c("metric", "grouping")) %>%
    group_by(metric, date, scenario, seed) %>%
    summarise(value = sum(value)) %>%
    ungroup() %>%
    mutate(grouping = "none",
           group    = NA) %>%
    bind_rows(raw_output) %>%
    setDT()
  
  return(agg_df)
}

# ---------------------------------------------------------
# Post process raw model output and append to 'result' list
# ---------------------------------------------------------
process_results = function(result, raw_output) {
  
  # NOTE: Datatable syntax used here for speed - can be an expensive process
  
  # Aggregate groupings for appropriate metrics
  agg_df = aggregate_results(result$input, raw_output)
  
  # Create key summary statistics across seeds
  result$output =
    agg_df[, .(mean   = mean(value),
               median = quantile(value, 0.5, na.rm = TRUE),
               lower  = quantile(value, o$quantiles[1], na.rm = TRUE),
               upper  = quantile(value, o$quantiles[2], na.rm = TRUE)),
           keyby = .(date, metric, grouping, group, scenario)
    ][order(metric, grouping, group, date)]
  
  # Metrics that can be cumulatively summed
  cum_metrics = result$input$metrics$df[cumulative == TRUE, metric]
  
  # Cumulatively sum all values
  cum_df = agg_df[metric %in% cum_metrics, ]
  cum_df[, cum_value := cumsum(value), by = .(metric, grouping, group, scenario, seed)]
  
  # ...then compute quantiles across seeds
  result$cum_output = 
    cum_df[, .(mean   = mean(cum_value), 
               median = quantile(cum_value, 0.5), 
               lower  = quantile(cum_value, o$quantiles[1]), 
               upper  = quantile(cum_value, o$quantiles[2])), 
           keyby = .(date, metric, grouping, group, scenario)
    ][order(metric, grouping, group, date)]
  
  return(result)
}

# ---------------------------------------------------------
# Format model outcomes ready for plotting
# ---------------------------------------------------------
format_results = function(o, f, results) {
  
  # Which results df to load depends on cumulative flag
  which_output = ifelse(f$cumulative, "cum_output", "output")
  
  # Reduce model output down to what we're interested in
  output_df = results[[which_output]] %>%
    select(date, metric, grouping, group, scenario, 
           value = o$best_estimate_simulation, lower, upper) %>%
    filter(is.na(date) | date >= f$plot_from, 
           is.na(date) | date <= f$plot_to, 
           metric %in% f$metrics, 
           !is.na(value))
  
  # ---- Do scaling and convert proportion to percentage ----
  
  # Subset of metric details: coverages and scaled metrics
  metric_df = results$input$metrics$df %>%
    select(metric, scale, coverage)
  
  # Scaler required (based on number simulated)
  scaler = f$person_days / results$input$population_size
  
  # For appropriate metrics, apply scaler and multiple by 100 
  output_df = output_df %>%
    left_join(metric_df, by = "metric") %>%
    mutate(value = ifelse(scale, value * scaler, value), 
           lower = ifelse(scale, lower * scaler, lower), 
           upper = ifelse(scale, upper * scaler, upper)) %>% 
    mutate(value = ifelse(coverage, value * 100, value), 
           lower = ifelse(coverage, lower * 100, lower), 
           upper = ifelse(coverage, upper * 100, upper)) %>%
    select(-scale, -coverage)
  
  # ---- Select appropriate grouping ----
  
  # Plotting by some disaggregation
  if (f$plot_type == "group") {
    
    # Filter for only this grouping
    output_df = output_df[grouping == f$plot_by, ]
    
    # Levels of factors
    group_levels = results$input$count[[f$plot_by]]
    output_df[, group := factor(group, group_levels)]
    
    # Plotting by age is a special case - group by age classification
    if (f$plot_by == "age")
      output_df = group_ages(o, output_df)
  }
  
  # No grouping - we want to remove all disaggregations
  if (f$plot_type != "group")
    output_df = output_df[is.na(group), ]
  
  # ---- Format output ----
  
  # Final formatting touches
  output_df = output_df %>%
    select(date, metric, group, value, lower, upper, scenario) %>%
    arrange(scenario, metric, group, date) %>% 
    mutate(metric = factor(metric, levels = f$metrics))
  
  return(output_df)
}

# ---------------------------------------------------------
# Group by user-defined age classifications
# ---------------------------------------------------------
group_ages = function(o, df, summarised = TRUE) {
  
  # Preallocate and bind age group variable to datatable
  age_df = df %>% 
    mutate(group = as.numeric(group), 
           age_group = "none")
  
  # Loop through age groups to classify into (from youngest to oldest)
  n_ages = length(o$plot_ages)
  for (i in 1 : n_ages) {
    
    # Upper and lower bound of this age group
    age_lower = sort(o$plot_ages)[i]
    age_upper = sort(o$plot_ages)[i + 1]
    
    # Name of this age classification
    if (i == n_ages) age_name = paste0(age_lower, "+")
    else age_name = paste0(age_lower, "-", age_upper)
    
    # Apply this grouping for all qualified individuals
    age_df[group >= age_lower, age_group := age_name]
  }
  
  # Recode with age group rather than age value
  age_df = select(age_df, -group) %>%
    rename(group = age_group) %>%
    mutate(group = paste("Age", group))
  
  # Dealing with already summarised output
  if (summarised == TRUE) {
    
    # Summarise by summing across age groups
    age_df = age_df %>%
      group_by(date, metric, scenario, group) %>%
      summarise(value = sum(value),
                lower = sum(lower),
                upper = sum(upper)) %>%
      ungroup() %>%
      setDT()
  }
  
  # Dealing with raw, non-summarised output
  if (summarised == FALSE) {
    
    # Summarise by summing across age groups but not simulations
    age_df = age_df %>%
      group_by(date, metric, seed, scenario, group) %>%
      summarise(value = sum(value)) %>%
      ungroup() %>%
      setDT()
  }
  
  return(age_df)
}

