###########################################################
# POST PROCESS
#
# Process model output ready for plotting.
#
###########################################################

# ---------------------------------------------------------
# Format model outcomes ready for plotting
# ---------------------------------------------------------
calculate_Reff = function(incidence, si, r_eff_window) {
  
  # See: https://cran.r-project.org/web/packages/EpiEstim/vignettes/demo.html
  
  # Number of daily observations
  n_days = length(incidence$local)

  # Use r_eff_window here for non-default time window
  r_win  = r_eff_window - 1
  r_days = seq(2, n_days - r_win) # Starting at 2 as conditional on the past observations

  # Check if serial interval has already been summarised
  if (is.list(si)) si_dist = si[c("mean", "std")]
  
  # If still a vector, calculate normal distribution parameters
  if (is.numeric(si)) si_dist = list(mean = mean(si), std = sd(si))
  
  # We'll use the mean and standard deviation of this data
  si_list = list(mean_si = si_dist$mean, 
                 std_si  = si_dist$std, 
                 t_start = r_days,
                 t_end   = r_days + r_win)
  
  # Use estimate_R function from EpiEstim package
  R_eff_df = estimate_R(incid  = incidence$local, 
                        method = "parametric_si",
                        config = make_config(si_list))
  
  # Extract R_eff estimate along with 95% CI bounds
  R_eff = R_eff_df$R %>%
    select(date  = t_start,
           value = "Mean(R)", 
           lower = "Quantile.0.025(R)", 
           upper = "Quantile.0.975(R)") %>%
    filter(date >= r_eff_window) %>%
    right_join(data.frame(date = 1 : n_days), 
               by = "date") %>%
    arrange(date) %>%
    mutate(incidence = incidence$local)
  
  # Append serial distribution 
  R_info = list(R_eff = R_eff, si = si_dist)
  
  return(R_info)
}

# ---------------------------------------------------------
# Format model outcomes ready for plotting
# ---------------------------------------------------------
format_results = function(o, f, results) {
  
  # Check if we're plotting a single simulation
  if ("value" %in% names(results$output)) {
    
    # If so, append trivial bounds
    output_df = results$output %>%
      mutate(lower = value, 
             upper = value, .after = value)
    
  } else {  # Otherwise select summary statistic of choice
    output_df = results$output %>%
      rename(value = o$best_estimate_simulation)
  }
  
  # Reduce model output down to what we're interested in
  model_df = output_df %>%
    filter(date >= f$plot_from, 
           date <= f$plot_to, 
           metric %in% f$metrics, 
           !is.na(value)) %>% 
    mutate(metric = factor(metric, levels = f$metrics))
  
  # ---- Scale metrics per n person days ----
  
  # Scaler required (based on number simulated)
  scaler = f$person_days / results$input$population_size
  
  # Vector of metrics to be scaled
  scale_metrics = results$input$metrics$df %>%
    filter(scaled == TRUE) %>%
    pull(metric)
  
  # Row indices of model_df to be scaled
  scale_idx = model_df$metric %in% scale_metrics
  
  # Apply this scaler to the relevant metrics
  model_df[scale_idx, ] = model_df[scale_idx, ] %>%
    mutate(value = value * scaler, value, 
           lower = lower * scaler, lower, 
           upper = upper * scaler, upper)
  
  # ---- Aggregate if desired ----
  
  # Plotting by some disaggregation
  if (f$plot_type == "group") {
    
    # Levels of factors
    group_levels = results$input$count[[f$plot_by]]
    
    # Filter for only such results and convert to factors
    model_df = model_df %>%
      filter(grouping == f$plot_by) %>% 
      mutate(group = factor(group, levels = group_levels)) %>% 
      select(-grouping)
    
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
  
  # ---- Format output ----
  
  # Check cumulative flag
  if (f$cumulative == TRUE) {
    
    # Cumulatively sum over time if requested
    model_df = model_df %>% 
      group_by(metric, scenario, group) %>%
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
    as.data.table()
  
  return(age_df)
}

