###########################################################
# POST PROCESS
#
# Process model output and data ready for plotting.
#
###########################################################

# ---------------------------------------------------------
# Format model outcomes ready for plotting
# ---------------------------------------------------------
format_results = function(o, f, a) {
  
  # Remove trivial model output before outbreak start
  plot_from = max(f$plot_from, o$dates_all[a$fit_data$outbreak_start])
  
  # Reduce model output down to what we're interested in
  model_df = a$model_output %>%
    rename(value = o$best_estimate_simulation) %>%
    filter(date >= plot_from, 
           date <= f$plot_to, 
           metric %in% f$metrics) %>% 
    mutate(metric = factor(metric, levels = f$metrics))
  
  # Plotting by some disaggregation - filter for only such results
  if (f$plot_type == "group") {
    model_df = filter(model_df, grouping == f$plot_by) %>% select(-grouping)
    
    # Plotting by age is a special case - group by age classification
    if (f$plot_by == "age")
      model_df = group_ages(o, model_df)
  }
  
  # No grouping - we want to remove all disaggregations
  if (f$plot_type != "group") {
    
    # Summarise first grouping per metric to avoid double counting
    first_grouping = model_df %>% 
      select(metric, grouping) %>%
      unique() %>%
      group_by(metric) %>%
      slice(1)
    
    # Sum across the groupings for each metric
    model_df = left_join(first_grouping, model_df, 
                         by = c("metric", "grouping")) %>%
      group_by(date, metric, scenario) %>%
      summarise(value = sum(value), 
                lower = sum(lower),
                upper = sum(upper)) %>%
      mutate(group = NA)
  }
  
  # Check cumulative flag
  if (f$cumulative == TRUE) {
    
    # Cumulatively sum over time if requested
    model_df = group_by(model_df, metric, scenario, group) %>%
      mutate(value = cumsum(value), 
             lower = cumsum(lower), 
             upper = cumsum(upper))
  }
  
  # Final formatting touches
  model_df = as.data.table(model_df) %>%
    select(date, metric, group, value, lower, upper, scenario) %>%
    arrange(scenario, metric, group, date)
  
  return(model_df)
}

# ---------------------------------------------------------
# Format data ready for plotting
# ---------------------------------------------------------
format_data = function(o, f, a) {
  
  # Check that data is appended
  if (is.null(a$fit_data))
    stop("Data not inherited from fitting process")
  
  # Grouping of data requested
  if (is.null(f$plot_by)) this_grouping = c("none", "na") else this_grouping = f$plot_by
  
  # Filter to leave only what we're interested in
  data_df = a$fit_data$epi %>%
    filter(grouping %in% this_grouping, 
           metric %in% f$metrics, 
           date >= f$plot_from,
           date <= f$plot_to, 
           !is.na(value)) %>%
    select(-grouping) %>%
    mutate(metric = factor(metric, levels = f$metrics))
  
  # If grouping is not trivial, set factors to preserve ordering
  if (!is.null(f$plot_by)) {
    
    # Levels of factors
    group_levels = o[[paste0("count_", f$plot_by)]]
    
    # Apply these levels
    data_df = mutate(data_df, group = factor(group, levels = group_levels))
  }
  
  # Check cumulative flag
  if (f$cumulative == TRUE) {
    
    # Cumulatively sum over time if requested
    data_df = group_by(data_df, metric, group) %>%
      mutate(value = cumsum(value)) %>%
      as.data.table()
  }
  
  # If no data identified, return NULL value
  if (nrow(data_df) == 0)
    data_df = NULL
  
  return(data_df)
}

# ---------------------------------------------------------
# Group by user-defined age classifications
# ---------------------------------------------------------
group_ages = function(o, df) {
  
  # Preallocate and bind age group variable to datatable
  age_df = mutate(df, group = as.numeric(group), age_group = "none")
  
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
  
  # Summarise by summing across age groups
  age_df = select(age_df, -group) %>%
    rename(group = age_group) %>%
    mutate(group = paste("Age", group)) %>%
    group_by(date, metric, scenario, group) %>%
    summarise(value = sum(value),
              lower = sum(lower),
              upper = sum(upper)) %>%
    as.data.frame()
  
  return(age_df)
}

