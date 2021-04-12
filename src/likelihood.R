###########################################################
# LIKELIHOOD
#
# Define P(model | data) likelihood function. Use primarily
# when training model emulator.
#
###########################################################

# ---------------------------------------------------------
# Likleihood function for actual model ouput - used to fit emulator
# ---------------------------------------------------------
likelihood_model = function(o, param_set = NULL, seed_num = NULL, m = NULL, 
                            data_from_cache = TRUE, quiet = FALSE) {
  
  # Load recently cached data (see load_data.R)
  #
  # NOTE: If running calibration test, syntehtic data is loaded here
  d = load_data(o, from_cache = data_from_cache, quiet = quiet)
  
  # Only need to run model for as long as we have data for
  if (!is.null(o$known_params)) short_run = FALSE else short_run = TRUE
  
  # Construct weighting of data points over time (see recent_data_weight option)
  time_weight    = max(min(1 - o$recent_data_weight, 1), 0)
  time_weight_df = data.table(date = o$dates_data, 
                              time_weight = seq(from = time_weight, to = 1, 
                                                length.out = length(o$dates_data)))
  
  # Construct weighting of different metrics (can be canton specific)
  metric_weight_df = o$calibration_multiplier %>%
    pivot_longer(cols = -canton, 
                 names_to = "metric", 
                 values_to = "metric_weight") %>%
    as.data.table()
  
  # We can really reduce metrics outputed by the model down to the essentials
  o$all_metrics = o$all_metrics %>%
    filter(metric %in% names(o$calibration_multiplier[-1])) %>%
    select(-grouping) %>%
    rename(grouping = calibrate_grouping)
  
  # Accordingly simplify metric grouping so update_output function applies correctly
  o$metric_groupings = setNames(o$all_metrics$grouping, o$all_metrics$metric)
  
  # Initiate likelihood value
  likelihood = 0
  
  # TODO: Reassess this loop if/when we start fitting to multiple cantons simultaneously...
  
  # Loop through cantons
  for (this_canton in o$cantons) {
    
    # Simulate model for this canton (unless it is provided)
    if (is.null(m)) {
      
      # Extract parameter values for this canton
      canton_idx    = o$calibration_df$scope %in% c("global", this_canton)
      canton_params = param_set[canton_idx]
      
      # Generate all model parameters (see parameters.R)
      p = get_parameters(o, d[[this_canton]], canton_params)
      
      # Run model (see model.R)
      m = model(o, p, seed = seed_num, short_run = short_run)
    }
    
    # Dates when we want to assess fit to data (ignore some points after epidemic outbreak)
    count_from = o$dates_all[d[[this_canton]]$outbreak_start + o$ignore_data_days]
    
    # Join model to data and calculate likelihood
    #
    # NOTE: Now using gamma distribution to have continuous support (rate parameter defaults to 1)
    likelihood_df = d[[this_canton]]$epi %>%
      filter(date >= count_from) %>%
      left_join(m, by = c("date", "metric", "grouping", "group")) %>%
      mutate(canton = this_canton, 
             data   = pmax(value.x, 1e-6), 
             model  = pmax(value.y, 1e-6)) %>%
      select(-value.x, -value.y) %>%
      mutate(likelihood = dgamma(data, model, log = TRUE))
    
    # Apply weightings to likelihood values
    weighted_df = likelihood_df %>%
      left_join(time_weight_df, by = "date") %>%
      left_join(metric_weight_df, by = c("canton", "metric")) %>%
      mutate(likelihood = likelihood * time_weight * metric_weight) %>%
      select(-canton, -time_weight, -metric_weight)
    
    # Continuously count total likelihood value
    likelihood = likelihood + sum(weighted_df$likelihood)
  }
  
  # # Evaluate parameter priors and increment
  # if (!is.null(prior_function)) 
  #   likelihood = likelihood + prior_function(param_set, o)
  
  # Output the model input and output as well as associated likelihood value
  result = list(param_set    = param_set, 
                model_input  = p, 
                model_output = m, 
                likelihood   = -likelihood)
  
  # NOTE: Returning a positive value to be minimised as this is the only way I can
  #       see to get the acquisition function of the adaptive sampling process working
  
  return(result)
}

