###########################################################
# SCENARIOS
#
# Alters model parameters for some predefined past or future 
# scenario.
#
# This function acts on the output of get_parameters 
# (see parameters.R).
#
###########################################################

# ---------------------------------------------------------
# Alter parameters to represent alternative scenarios.
# ---------------------------------------------------------
alter_model_params = function(o, p, scenario) {
  
  # Psuedo 'today' - last date of data
  today = max(o$dates_data)
  
  # ---- Past scenario: No historical interventions ----
  
  # Assume we had no interventions in the past
  if (scenario == "past_evaluation")
    p$npi_level[] = 0
  
  # ---- Past scenario: Earlier response (by 2 weeks) ----
  
  # Start the same response 14 days before actual implementation
  
  # Shift response timing
  time_diff = 14
  
  # Shift timing of past interventions
  if (scenario == "earlier_npi") {
    
    # Cut values to start, add values from end
    cut_dates = p$npi_level[(time_diff + 1) : o$n_dates]
    add_dates = rep(p$npi_level[o$n_dates], time_diff)
    
    # Concatenate for new contact reduction timing
    p$npi_level = c(cut_dates, add_dates)
  }
  
  # ---- Past scenario: Later response (by 2 weeks) ----
  
  # Shift response timing
  time_diff = 14
  
  # Shift timing of past interventions
  if (scenario == "later_npi") {
    
    # Add values to start, cut values from end
    add_dates = rep(0, time_diff)
    cut_dates = p$npi_level[1 : (o$n_dates - time_diff)]
    
    # Concatenate for new contact reduction timing
    p$npi_level = c(add_dates, cut_dates)
  }
  
  # ---- Future scenarios: Relax all interventions in two weeks ----
  
  # Alter parameters from two weeks time until end of analysis
  time_idx = future_time_index(o, from = today + 14) 
  
  # Reduce npi_level effect by 20%
  if (scenario == "relax_npi")
    p$npi_level[time_idx] = p$npi_level[time_idx] * (1 - 0.2)
  
  # ---- Future scenario: Totally remove all interventions in two weeks ----
  
  # NOTE: People revert back to 'normal' behaviour
  
  # Alter parameters from two weeks time until end of analysis
  time_idx = future_time_index(o, from = today + 14) 
  
  # Totally remove npi_level effect
  if (scenario == "remove_npi")
    p$npi_level[time_idx] = 0
  
  # ---- Example vaccine modelling scenarios ----
  
  # Vaccinate 10k people everyday from tomorrow
  if (scenario == "vaccine_10k")
    p$number_vaccines[future_time_index(o)] = 10000
  
  # Vaccinate 25k people everyday from tomorrow
  if (scenario == "vaccine_25k")
    p$number_vaccines[future_time_index(o)] = 25000
  
  # ---- Sanity checks ----
  
  # Sanity check on length of npi_level vector
  if (length(p$npi_level) != o$n_dates)
    stop("We seem to have gained/lost time")
  
  return(p)
}

# ---------------------------------------------------------
# Get indices of time points for altering parameters.
# ---------------------------------------------------------
future_time_index = function(o, from = max(o$dates_data) + 1, to = NULL) {
  
  # Default end point as final date of analysis
  if (is.null(to)) to = max(o$dates_all)
  
  # Ensure any specified end point is not after last point of simulation
  to = min(to, max(o$dates_all))
  
  # Consider trivial case where first point is after end point of analysis
  if (from > max(o$dates_all)) time_idx = NULL
  
  else # All indices between 'to' and 'from' dates
    time_idx = which(o$dates_all == from) : which(o$dates_all == to)
  
  return(time_idx) 
}

# ---------------------------------------------------------
# Define descriptive names associated with scenarios
# ---------------------------------------------------------
scenario_names = function(o) {
  
  # Named vector of scenario IDs and descriptive names to display on plots
  scenario_dict = qc(past_evaluation = "No previous interventions", 
                     earlier_npi     = "NPIs two weeks earlier", 
                     later_npi       = "Response measures two weeks later", 
                     relax_npi       = "20% decrease in response measures",  
                     remove_npi      = "Remove all interventions tomorrow", 
                     vaccine_10k     = "Vaccinate 10k people each day", 
                     vaccine_25k     = "Vaccinate 25k people each day")
  
  return(scenario_dict)
}

