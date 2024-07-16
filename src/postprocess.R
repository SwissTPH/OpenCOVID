###########################################################
# POST PROCESS
#
# Process model output ready for plotting.
#
###########################################################

# ---------------------------------------------------------
# Calculate effective reproduction number
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
# Aggregate groupings for appropriate metrics
# ---------------------------------------------------------
aggregate_results = function(input, raw_output) {
  
  # NOTE: This is the first step when processing raw results.
  #       It is in a seperate function because we may want to
  #       only aggregate, and not summarise raw results.
  
  # Metrics that can be aggregated
  agg_metrics = input$metrics$df[aggregate == TRUE,  metric]
  
  # Need to aggregate grouped metrics: get first grouping per metric to avoid double counting
  first_grouping = raw_output %>% 
    filter(metric %in% agg_metrics) %>% 
    select(metric, grouping) %>%
    unique() %>%
    filter(!grouping %in% c("none", "na")) %>%
    group_by(metric) %>%
    slice(1)
  
  # Aggregate values of first grouping
  agg_df = first_grouping %>%
    left_join(raw_output,
              by = c("metric", "grouping")) %>%
    group_by(metric, date, scenario, seed) %>%
    summarise(value = sum(value)) %>%
    mutate(grouping = "none",
           group    = NA) %>%
    bind_rows(raw_output) %>%
    as.data.table()
  
  return(agg_df)
}

# ---------------------------------------------------------
# Post process raw model output and append to 'result' list
# ---------------------------------------------------------
process_results = function(result, raw_output) {
  
  # Aggregate groupings for appropriate metrics
  agg_df = aggregate_results(result$input, raw_output)
  
  # Create key summary statistics across seeds
  result$output = agg_df %>% 
    group_by(date, metric, grouping, group, scenario) %>% 
    summarise(mean   = mean(value),
              median = quantile(value, 0.5, na.rm = TRUE),
              lower  = quantile(value, o$quantiles[1], na.rm = TRUE),
              upper  = quantile(value, o$quantiles[2], na.rm = TRUE)) %>% 
    arrange(metric, grouping, group, date) %>% 
    as.data.table()
  
  # Metrics that can be cumulatively summed
  cum_metrics = result$input$metrics$df[cumulative == TRUE, metric]
  
  # Cumulatively sum all values then compute quantiles across seeds
  result$cum_output = agg_df %>%
    filter(metric %in% cum_metrics) %>%
    group_by(metric, grouping, group, scenario, seed) %>%
    mutate(cum_value = cumsum(value)) %>%  
    group_by(date, metric, grouping, group, scenario) %>% 
    summarise(mean   = mean(cum_value),
              median = quantile(cum_value, 0.5),
              lower  = quantile(cum_value, o$quantiles[1]),
              upper  = quantile(cum_value, o$quantiles[2])) %>% 
    arrange(metric, grouping, group, date) %>% 
    as.data.table()
  
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
    select(metric, scaled, coverage)
  
  # Scaler required (based on number simulated)
  scaler = f$person_days / results$input$population_size
  
  # Variables that need to be mutated
  vars = c("value", "lower", "upper")
  
  # For appropriate metrics, apply scaler and multiple by 100 
  output_df = output_df %>%
    left_join(metric_df, by = "metric") %>%
    mutate(across(all_of(vars), function(x) ifelse(scaled, x * scaler, x)), 
           across(all_of(vars), function(x) ifelse(coverage, x * 100, x))) %>%
    select(-scaled, -coverage)
  
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

# ---------------------------------------------------------
# Calculate network summary statistics
# ---------------------------------------------------------
calculate_network_statistics <- function(result, raw_df, metrics = c("all_new_infections", "hospital_admissions")) {
  
  skip_days <- 100
  
  # ---- Calculate peak metrics ----
  
  raw_df <- raw_df[metric %in% metrics & grouping == "none"]
  
  raw_df_trim <- copy(raw_df[date > skip_days, ])
  
  # Select highest value by group
  descriptives <- raw_df_trim[raw_df_trim[order(date), .I[value == max(value)], by=.(metric, scenario, seed)]$V1][order(date), .SD, by = .(metric, scenario, seed)]
  
  # Select always the first entry - To discuss!
  descriptives <- descriptives[descriptives[, .I[1], by = .(metric, scenario, seed)]$V1]
  
  descriptives <- descriptives[, .(metric, scenario, seed, peak_date = date, peak_value = value)]
  
  cum_df <- copy(raw_df)
  
  cum_df[order(date), cum_value := cumsum(value), by = .(metric, scenario, seed)]
  
  cum_df <- cum_df[date == max(date), .SD, by = .(metric, scenario, seed)]
  
  results_df <- cum_df[descriptives, on = .(metric, scenario, seed)][, .(metric, scenario, seed, peak_date, peak_value, cum_value)]
  
  results_df <- pivot_longer(results_df, c("peak_date", "peak_value", "cum_value"), names_to = "variable") %>%
    as.data.table()
  
  stats_df <- results_df[, list(mean = mean(value), 
                                  median = median(value), 
                                  lower = quantile(value, probs = o$quantiles[1]), 
                                  upper = quantile(value, probs = o$quantiles[2])), 
                           by = list(metric, scenario, variable)]
  
  result$clustering$stats <- stats_df
  
  # ---- Calculate network metrics (e.g., clustering) ----
  
  network <- result$network %>% select(from, to, layer)
  
  layer_weights <- result$input$layer_properties$beta_scaler %>%
    as.data.frame() %>%
    pivot_longer(everything(), names_to = "layer", values_to = "weight") %>%
    as.data.table()
  
  # Join layer weight from layer_weights table
  network <- layer_weights[network[, .(from, to, layer)], on = .(layer)][, .(i = from, j = to, w = weight)]
  
  # Remove duplicate entries
  network <- network[, .(w, .N), by = .(i, j)][N == 1, .(i, j, w)]
  
  # Calculate "simple" clustering coefficient (igraph package)
  result$clustering$cc <- network[, .(from = i, to = j)] %>%
    graph_from_data_frame(directed = FALSE) %>%
    transitivity()
  
  # Calculate generalized clustering coefficient (tnet package)
  cc_generalized <- network %>%
    as.tnet(type = "weighted one-mode tnet") %>%
    clustering_w()
  
  # Use arithmetic mean only
  result$clustering$cc_generalized <- cc_generalized[["am"]]
  
  return(result)
  
}

