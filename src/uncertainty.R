###########################################################
# UNCERTAINTY
#
# Sample parameter sets to represent parameter uncertainty.
#
###########################################################

# ---------------------------------------------------------
# Store uncertainty parameters in seperate list, denoted u
# ---------------------------------------------------------
parse_uncertainty = function(o, y, u = NULL) {
  
  # Initiate u list if it doesn't current exist
  if (is.null(u)) u = list()
  
  # All individual parameter items (see auxiliary.R)
  param_items = names(unlist_format(y))
  
  # Which of the these are to be sampled for uncertainty
  uncert_idx = grep("\\$uncertainty\\$", param_items)
  
  # Full names of these parameters
  uncert_params = param_items[uncert_idx] %>%
    str_remove("\\$uncertainty\\$.*") %>%
    unique()
  
  # Loop through uncretainty parameters
  for (uncert_param in uncert_params) {
    
    # Store uncertainty details in u list
    eval_str("u[[uncert_param]] = y$", uncert_param, "$uncertainty")
    
    # Remove from y list, which is about to be parsed
    eval_str("y$", uncert_param, " = NULL")
  }
  
  return(list(y, u))
}

# ---------------------------------------------------------
# Sample points across parameter uncertainty space
# ---------------------------------------------------------
sample_uncertainty = function(o) {
  
  # Extract all parameters we want to generate uncertainty for
  u = parse_yaml(o, "all", uncert = "*read*")
  
  # Return out if empty
  if (is_empty(u))
    return()
  
  # Create unit hypercube (one dimension per parameter)
  unit_cube = t(matrix(c(0, 1), nrow = 2, ncol = length(u)))
  
  # Sample unitcube using Latin Hypercube sampling
  lhc_samples = o$n_parameter_sets %>%
    lhs(unit_cube) %>%
    as_named_dt(names(u))
  
  # Now transform to appropriate scale and distribution for each parameter...
  
  # Convert function names (by adding 'q' for qnorm, qunif, etc)
  u = lapply(u, function(x) list_modify(x, fn = paste0("q", x$fn)))
  
  # Loop through parameters
  uncert_list = NULL
  for (param in names(u)) {
    
    # Append values sampled from unitcube to function call list
    param_fn = list.prepend(u[[param]], p = pull(lhc_samples, param))
    
    # Apply distribution of choice to unitcube values
    uncert_list[[param]] = parse_fn(param_fn)
  }
  
  # Convert parameter set list to long datatable
  uncert_df = uncert_list %>%
    as.data.table() %>% 
    mutate(param_set = 1 : o$n_parameter_sets) %>%
    pivot_longer(cols = -param_set, 
                 names_to = "param") %>%
    arrange(param_set, param) %>%
    setDT()
  
  # Save to file for use on cluster
  saveRDS(uncert_df, file = paste0(o$pth$uncertainty, "uncertainty.rds"))
  
  return(uncert_df) 
}

# ---------------------------------------------------------
# Sample single average value from uncertainty distribution(s)
# ---------------------------------------------------------
sample_average = function(o, average_fn = "median") {
  
  # Extract all parameters for which uncertainty has been defined
  u = parse_yaml(o, "all", uncert = "*read*")
  
  # Return out if empty
  if (is_empty(u))
    return()
  
  # We want to sample a load of points and take the average...
  
  # Initiate output
  uncert_list = NULL
  
  # Convert function names (by adding 'r' for rnorm, runif, etc)
  u = lapply(u, function(x) list_modify(x, fn = paste0("r", x$fn)))
  
  # Loop through parameters
  for (param in names(u)) {
    
    # Append number of values to randomly sample from this distribution
    param_fn = list.prepend(u[[param]], n = 1e5)
    
    # Parse the list and take the average
    uncert_list[[param]] = get(average_fn)(parse_fn(param_fn))
  }
  
  return(uncert_list) 
}

# ---------------------------------------------------------
# Apply selected parameter value to yaml list for parsing
# ---------------------------------------------------------
apply_uncertainty = function(y, uncert) {
  
  # Loop through parameters to alter
  for (param in names(uncert)) {
    
    # # Check whether parameter we want to overwrite already exists
    # eval_str("param_exist = is.numeric(y$", param, ")")
    # 
    # # Check parameter is valid - throw error if not
    # if (!param_valid)
    #   stop("Attempting to alter invalid parameter '", param, "' during model fitting")
    
    # Re-define parameter value (using eval to enable change of listed items)
    eval_str("y$", param, " = uncert[[param]]")
  }
  
  return(y)
}

