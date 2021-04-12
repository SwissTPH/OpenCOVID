###########################################################
# STRATEGY
#
# Define intervention 'strategies'. These are essentially more
# elaborate versions of 'scenarios'.
#
# Two stategy defintions required: a 'roll_out' of vaccines and an 
# 'npi_exit' plan. Options for each of this defintions is as follows.
#
# ---- roll_out ----
#
# REQUIRED PARAMETERS:
#  coverage := Vaccine coverage acheived for the priority group
#  start    := Start date of vaccine roll out for the priority group
#  duration := Number of days until coverage is acheived for the priority group
#  vaccine  := Which vaccine is used for the priority group (see vaccine_properties.csv)
#
# NOTE: See 'vaccine_priority' in options.R for defintion of priority groups.
#
# ---- npi_exit ----
#
# REQUIRED PARAMETERS:
#  start := Date measures begin to relax
#
# THEN ONE OF:
#  new_normal_rel := Relative NPI level after delay period
#  new_normal_abs := Absolute NPI after delay period
#  new_normal_osi := Set level of OSI directly (NPI scaling handled internally)
#
# OPTIONAL PARAMETERS:
#  delay := Number of days after measures are relaxed to get back to 'new normal'
#
# NOTE: For extrapolating current NPI levels, use 'new_normal_rel = 1'
#
###########################################################

# Example strategy
n1a = function() {
  
  # NPI exit strategy
  npi_exit = list(
    npi1 = list(name   = "No change in NPI policy", 
                phase1 = list(start = "2021-04-26", new_normal_rel = 1)), 
    npi2 = list(name   = "Example NPI ploicy", 
                phase1 = list(start = "2021-04-26", new_normal_osi = 0.55), 
                phase2 = list(start = "2021-05-24", new_normal_osi = 0.55)))
  
  # Vaccine rollout method 1: coverage targets per vaccine priority group
  # roll_out = list(
  #   fastvax = list(name   = "Example vaccine rollout", 
  #                  group1 = list(coverage = 0.75, start = "2021-01-01", duration = 60, vaccine = "pfizer"),
  #                  group2 = list(coverage = 0.75, start = "2021-03-01", duration = 60, vaccine = "pfizer"),
  #                  group3 = list(coverage = 0.75, start = "2021-06-01", duration = 90, vaccine = "oxford")))
  
  # Vaccine rollout method 2: number vaccinated per day (in order of priority)
  alter_params = c("p$number_vaccines[o$dates_all > max(o$dates_data) + 1] = 20000")
  roll_out = list(defvax = list(name = "Example NPI ploicy", 
                                group1 = list(coverage = 0, start = "2021-01-01", duration = 0)))
  
  # Format strategy details into single list
  details = list(roll_out = roll_out, npi_exit = npi_exit, alter_params = alter_params)
  
  return(details)
}

# ---------------------------------------------------------
# Parse strategy and manipulate parameters accordingly
# ---------------------------------------------------------
parse_strategy = function(o, p, details, this_combination) {
  
  # ---- Vaccine roll out to each priority group ----
  
  # Extract roll out strategy
  this_roll_out = details$roll_out[[this_combination$roll_out]]
  
  # Initiate named vector for which vaccine to be used for each group
  vaccine_type   = rep(o$vaccine_default, length(o$vaccine_priority))
  p$vaccine_type = setNames(vaccine_type, names(o$vaccine_priority))
  
  # Iterate through the groups to vaccinate
  for (group in names(this_roll_out)) {
    this_group = this_roll_out[[group]]
    
    # Index of priority group (see o$vaccine_priority)
    group_idx = which(group == names(o$vaccine_priority))
    
    # Scale up vector (assuming we have all time points available)
    scale_up = seq(0, this_group$coverage, length.out = this_group$duration + 1)
    
    # Get the proposed start and end dates for roll out to this priority group
    start_date = format_date(this_group$start)
    end_date   = start_date + this_group$duration
    
    # Check we are at least starting at a valid time point
    if (start_date > max(o$dates_all))
      stop("You are attempting to start vaccine deployment too far into the future", 
           " - try increasing 'n_future_days'")
    
    # Indices of model-valid time points to implement vaccine roll out
    time_idx = future_time_index(o, from = start_date, to = end_date) 
    
    # Linearly increase target coverage over duration of roll out
    p$vaccine_coverage[time_idx, group_idx] = scale_up[1 : length(time_idx)]
    
    # Vaccine coverage then stays stable at all later points
    p$vaccine_coverage[future_time_index(o, from = end_date + 1), group_idx] = this_group$coverage
    
    # If defined, set vaccine to be used (overwriting default)
    if (!is.null(this_group$vaccine))
      p$vaccine_type[group_idx] = this_group$vaccine
  }
  
  # Check that all vaccine names are well defined
  valid_type = p$vaccine_type %in% names(o$vaccine_dict)
  if (!all(valid_type))
    stop("Invalid vaccine name: ", paste(p$vaccine_type[!valid_type], collapse = ", "))
  
  # ---- NPI exit strategy ----
  
  # Extract NPI exit strategy
  this_npi_exit = details$npi_exit[[this_combination$npi_exit]]
  
  # Iterate through the release phases
  for (phase in names(this_npi_exit)) {
    this_phase = this_npi_exit[[phase]]
    
    # Extract number of days needed tp reach new normal (defaults to immediate)
    this_delay = max(this_phase$delay, 1)
    
    # Get the proposed start and end dates for this NPI release phase
    start_date = format_date(this_phase$start)
    end_date   = start_date + this_delay
    
    # Check we are at least starting at a valid time point
    if (start_date <= max(o$dates_data))
      stop("You are attempting to start an NPI release phase in the past")
    
    # Check we are at least starting at a valid time point
    if (start_date > max(o$dates_all))
      stop("You are attempting to start an NPI release phase too far into the future", 
           " - try increasing 'n_future_days'")
    
    # Level of NPI at starting point
    init_value = p$npi_level[which(o$dates_all == start_date)]
    
    # The 'new normal' can be defined in absolute or relative terms, or by defining new OSI level
    if (!is.null(this_phase$new_normal_rel)) new_normal = this_phase$new_normal_rel * init_value
    if (!is.null(this_phase$new_normal_abs)) new_normal = this_phase$new_normal_abs
    if (!is.null(this_phase$new_normal_osi)) new_normal = this_phase$new_normal_osi * p$npi_scaler
    
    # Check we are at least starting at a valid time point
    if (!exists("new_normal"))
      stop("Neither 'new_normal_abs' or 'new_normal_rel' defined for NPI exit strategy")
    
    # Scale down vector (assuming we have all time points available)
    scale_down = seq(init_value, new_normal, length.out = this_delay + 1)
    
    # Indices of model-valid time points to implement vaccine roll out
    time_idx = future_time_index(o, from = start_date, to = end_date) 
    
    # Linearly increase target coverage over duration of roll out
    p$npi_level[time_idx] = scale_down[1 : length(time_idx)]
    
    # Contact reduction is then the new normal for all later points
    p$npi_level[future_time_index(o, from = end_date + 1)] = new_normal
  }
  
  # ---- Alter other model parameters ----
  
  # NOTE: These changes are applied for ALL combinations in the strategy
  
  # Loop through and evaluate string to apply condition
  for (i in seq_along(details$alter_params))
    eval(parse(text = details$alter_params[i]))
  
  return(p)
}

