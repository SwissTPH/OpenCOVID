###########################################################
# MODEL
#
# Simulate transmission model for a given parameter set
# (defined by p). Outputs a long datatable, denoted m, and
# several other model properties of interest.
#
# Verbose options: 
#  date - print each date as we iterate
#  bar  - report using progress bar
#  none - no progress reporting
#
###########################################################

# ---------------------------------------------------------
# Main function: The guts of the model
# ---------------------------------------------------------
model = function(o, scenario, seed = NA, fit = NULL, uncert = NULL, do_plot = FALSE, verbose = "date") {
  
  # Initiate model timer
  tic("model")
  
  # Set random seed generator if input defined
  if (!is.na(seed)) set.seed(seed)
  
  # ---- Generate model parameters ----
  
  if (verbose != "none") message(" - Parsing input")
  
  # Load and parse user-defined inputs for this scenario
  yaml = parse_yaml(o, scenario, fit = fit, uncert = uncert, read_array = TRUE)
  
  # Shorthand for model parameters
  p = yaml$parsed
  
  # ---- Model set up ----
  
  if (verbose != "none") message(" - Running model")
  
  # Generate people datatable
  ppl_df = create_ppl(p, verbose = verbose)
  
  # Generate full contact network outside of main model loop (see networks.R)
  network = create_network(p, ppl_df, do_plot, verbose)
  
  # Initialise infections/vaccinations
  ppl_df = initiate_epidemic(p, ppl_df)
  
  # Calculate number of people in each vaccine prioirty group outside of main loop
  n_priority_group = ppl_df %>%
    count(priority_group) %>%
    left_join(p$priority_groups[, .(id, priority)], 
              by = c("priority_group" = "id")) %>% 
    arrange(priority)
  
  # Group certain states for easy access and readability
  states = list(infected   = qc(latent, presym, asym, mild, severe, crit), 
                infectious = qc(presym, asym, mild, severe, crit), 
                in_care    = qc(hospital, icu, rehosp))
  
  # Also append vectors for all disease and care states
  for (state in p$model_states$all)
    states[[state]] = rep(0, p$n_days)
  
  # Preallocate model output datatable (denoted 'm' for short)
  m = preallocate_output(p)
  
  # Preallocate list for collecting epidemiological statistics
  epi_stats = list(incidence  = rep(0, p$n_days), 
                   serial_int = vector("list", p$n_days))  # For calculation of Re
  
  # Preallocate list for healthcare and intervention costs
  cost_units = list(vaccine   = rep(0, p$n_days), 
                    prep      = rep(0, p$n_days),
                    treatment = rep(0, p$n_days), 
                    hospital  = rep(0, p$n_days), 
                    icu       = rep(0, p$n_days))
  
  # ---- Main model loop ----
  
  if (verbose != "none")
    message("  > Simulating for ", p$n_days, " days")
  
  # Initiate progress bar if desired
  if (verbose == "bar")
    pb = start_progress_bar(p$n_days)
  
  # Iterate through time steps - starting from first imported case(s)
  for (i in 1 : p$n_days) {
    
    # Report progress by printing date (useful when running on the cluster)
    if (verbose == "date") 
      message("   - Simulating day ", i)
    
    # Copy people datatable for the next time step
    ppl = data.table::copy(ppl_df)
    
    # Calculate basic epidemiological indicators
    list[m, epi_stats] = fn_epi_stats(m, p, ppl, epi_stats, i)
    
    # ---- New infections ----
    
    # Calculate which infected-susceptible contacts lead to transmission
    list[m, ppl, id] = fn_transmission(m, p, ppl, network, epi_stats, i)
    
    # Import infections acquired elsewhere
    list[m, ppl, id] = fn_importation(m, p, ppl, epi_stats, id, i)
    
    # New viral variant emergence
    list[m, ppl, id] = fn_emergence(m, p, ppl, epi_stats, id, i)
    
    # Count incidence and serial intervals for all new infections
    list[m, epi_stats] = fn_incidence(m, p, ppl, epi_stats, id, i)
    
    # ---- Testing, diagnosis, and isolation ----
    
    # Determine test and diagnosis probability for newly symptomatic people
    list[m, ppl] = fn_test_diagnose(m, p, ppl, i)
    
    # Isolate newly diagnosed individuals where appropriate
    list[m, ppl] = fn_isolate(m, p, ppl, i)
    
    # ---- Vaccination, PrEP, and treatment ----
    
    # Vaccinate eligible priority groups according
    list[m, ppl, cost_units] = fn_vaccinate(m, p, ppl, i, n_priority_group, cost_units)
    
    # Pre-exposure prophylaxis for those eligible
    list[m, ppl, cost_units] = fn_prep(m, p, ppl, i, cost_units)
    
    # Treat eligible infected cases
    list[m, ppl, cost_units] = fn_treatment(m, p, ppl, i, cost_units)
    
    # ---- Epidemiological and population updates ----
    
    # Update disease progression for the next time step and count disease states
    list[m, ppl, cost_units, states] = fn_epi_update(m, p, ppl, i, cost_units, states)
    
    # Non-COVID-related 'natural' deaths
    list[m, ppl] = fn_natural_death(m, p, ppl, i)
    
    # Replace deceased (COVID or otherwise) with newborns
    list[m, ppl] = fn_birth(m, p, ppl, i)
    
    # Age individuals each year
    list[m, ppl] = fn_ageing(m, p, ppl, i)
    
    # Update ppl datatable
    ppl_df = data.table::copy(ppl)
    
    # Update progress bar
    if (verbose == "bar") 
      setTxtProgressBar(pb, i)
  }
  
  # ---- Close up ----
  
  # Close progress bar
  if (verbose == "bar")
    close(pb)
  
  # Produce plot of disease and care states over time
  if (do_plot == TRUE)
    plot_disease_state(o, p, "Population disease and care states", states)
  
  # Calculate healthcare and intervention costs
  m = health_economics(m, p, cost_units)
  
  # Calculate effective reproduction number from incidence and serial interval
  si = unlist(epi_stats$serial_int)
  Re_info = calculate_Re(epi_stats$incidence, si, p$Re_window)
  
  # Prepare final model output
  results = format_output(m, p, Re_info, network, ppl, yaml, seed)
  
  # Model runtime (in seconds)
  time_clock = toc(quiet = TRUE)
  time_taken = round(time_clock$toc - time_clock$tic)
  
  # Store the time taken
  results$time_taken = seconds_to_period(time_taken)
  
  # Display model run time if appropriate
  if (verbose != "none")
    message("  > Model runtime: ", results$time_taken)
  
  return(results)
}

# ---------------------------------------------------------
# Generate ppl datatable by creating new individuals
# ---------------------------------------------------------
create_ppl = function(p, n = NULL, init = TRUE, verbose = "none") {
  
  # Number of new people to create (on initial call use population_size)
  if (init == TRUE) n = p$population_size
  
  # Throw an error if trying to create trivial number of people
  if (is.null(n) || n == 0)
    stop("Attemping to create zero people")
  
  # Display (or not) have many people we are creating
  if (verbose != "none")
    message("  > Initiating population of ", thou_sep(n))
  
  # Initiate new datatable
  ppl = data.table()
  
  # Preallocate variables with default values of appropriate class (see config/model_metrics.yaml)
  lapply(p$model_vars, function(n, x) {
    ppl[, (x$name) := get(x$class)(n)]
    ppl[, (x$name) := x$value]}, n = n)
  
  # ---- Set initial values ----
  
  # Individual IDs are essentially row numbers
  ppl[, id := 1 : n]
  
  # Sample ages from some distribution and a birthday index
  sample_ages = sample_vec(x = p$ages,  size = n, replace = TRUE, prob = p$demography)
  sample_bday = sample_vec(x = 0 : 364, size = n, replace = TRUE)
  
  # Apply age (use zero for newborns)
  if (init == TRUE)  ppl[, age := sample_ages]
  if (init == FALSE) ppl[, age := 0]
  
  # Apply birthday index - do not use 0 for newborns to stagger ageing
  ppl[, birthday := sample_bday]
  
  # ---- Assign risk groups ----
  
  # Sanity check that all defined risk groups are able to modelled
  if (!all(names(p$risk_groups) %in% names(ppl)))
    stop("Unrecognised risk groups")
  
  # Loop through the user-defined risk groups
  for (this_group in names(p$risk_groups)) {
    
    # Probability of acceptance per person (dependent on age)
    risk_prop = p$risk_groups[[this_group]]$probability[ppl[, age] + 1]
    
    # Sample TRUE and FALSE for all in this age group accordingly 
    ppl[, (this_group) := runif(n) < risk_prop]
  }
  
  # ---- Set vaccination priority policy and acceptance ----
  
  # Iterate through the different priority groups
  # 
  # NOTE: We do last at first so any people in multiple groups are placed in the highest possible group
  for (group_id in rev(p$priority_groups$id)) {
    this_group = p$priority_groups[id == group_id, ]
    
    # Evaluate the condition and set priority for any that satisfy
    ppl[eval_str(this_group$condition), priority_group := this_group$id]
  }
  
  # Remove all those that are no able to receive vaccination (but can receive PrEP)
  ppl[vax_unsuitable == TRUE, priority_group := "none"]
  
  # Again iterate through priority groups - this time to set vaccine acceptance probabilities
  #
  # NOTE: We do this after all groups are formalised to prevent double counting
  for (group_id in p$priority_groups$id) {
    
    # Proportion of this group that would accept a vaccine (ie max coverage)
    acceptance = p$vaccine$details[id == group_id, coverage]
    booster    = p$booster_details[id == group_id, probability]
    
    # Sample logical and populate vaccine_accept variable
    ppl[priority_group == group_id, 
        vaccine_accept := sample_vec(c(TRUE, FALSE), size = .N, replace = TRUE, 
                                     prob = c(acceptance, 1 - acceptance))]
    
    # Do similar for booster_accept variable, but only for those that accept vaccine in general
    ppl[priority_group == group_id & vaccine_accept == TRUE, 
        booster_accept := sample_vec(c(TRUE, FALSE), size = .N, replace = TRUE, 
                                     prob = c(booster, 1 - booster))]
  }
  
  # ---- Mass testing acceptance ----
  
  # Number of modelled people willing to be regularly mass tested
  n_mass_test_accept = round(n * p$testing$mass_testing$probability)
  
  # Probability of acceptance per person (dependent on age)
  p_mass_test_accept = p$testing$mass_testing$age_prob[ppl[, age] + 1]
  
  # Sample mass_test_acceptance of the population with age-related probabilities
  accept_id = sample_int_crank(n, n_mass_test_accept, p_mass_test_accept)
  
  # Assign these people to accept mass testing
  ppl[accept_id, mass_test_accept := TRUE]
  
  # ---- Tidy up ----

  # Remove all attributes
  ppl = ppl[, lapply(.SD, as.vector)] 
  
  return(ppl)
}

# ---------------------------------------------------------
# Initially infect a subset of people to kick start epidemic
# ---------------------------------------------------------
initiate_epidemic = function(p, ppl) {
  
  # ---- Initiate previously vaccinated ---
  
  # Iterate through the vaccine and treatment priority groups
  for (priority_id in p$priority_groups$id) {
    group = p$vaccine$details[id == priority_id, ]
    
    # Skip this if we have zero coverage for this group
    if (group$coverage > 1e-6) {
      
      # IDs of all those in this priority group who will accept the vaccine
      group_id = ppl[priority_group == priority_id & vaccine_accept == TRUE, id]
      
      # Absolute number vaccinated in the past
      n_vaccinate = round(length(group_id) * p$vaccine$coverage_init[[group$id]])
      
      # Sample IDs of those previously vaccinated - uniformly for all in this priority group
      vaccinate_id = sample_vec(group_id, size = n_vaccinate)
      
      # Sample number of days each person has been vaccinated for
      ppl[vaccinate_id, days_vaccinated := -sample_vec(group$start : min(group$end, -1), 
                                                       size = .N, replace = TRUE)]
      
      # Assign the 'initial' vaccine to these individuals
      ppl[vaccinate_id, vaccine_type := "init"]
    }
  }
  
  # Count primary vaccine doses
  for (dose_day in c(0, p$vaccine$subsequent_dose_days))
    ppl[days_vaccinated > dose_day, vaccine_doses := vaccine_doses + 1L]
  
  # ---- Set booster dose info for those already vaccinated ----
  
  # Number of doses in primary vaccine schedule
  primary_doses = length(p$vaccine$subsequent_dose_days) + 1
  
  # Booster delivery details for all those who will recieve booster dose(s)
  booster_df = ppl %>%
    filter(vaccine_doses >= primary_doses, 
           booster_accept == TRUE) %>% 
    select(id, priority_group, days_vaccinated, days_booster, booster_due) %>% 
    left_join(y  = p$booster_details[, -"probability"], 
              by = c("priority_group" = "id")) %>% 
    mutate(n_booster_doses = 0)
  
  # Assign first booster delivery date, bounded below if force-starting
  booster_df[, booster_due := pmin(pmax(cycle_period - days_vaccinated, start), force_start, na.rm = TRUE)]
  
  # Count the number of boosters received up until initialise date
  booster_df[booster_due < 0, n_booster_doses := -booster_due %/% cycle_period]
  
  # If first booster is due in the past, determine most recent dose and assign date of next dose 
  booster_df[booster_due < 0, days_booster := -booster_due - n_booster_doses * cycle_period]
  booster_df[booster_due < 0, booster_due := cycle_period - days_booster]
  
  # Apply these details to the main ppl datatable (avoid starting on day zero)
  ppl[id %in% booster_df$id, days_booster := pmax(booster_df$days_booster, 1, na.rm = FALSE)]
  ppl[id %in% booster_df$id, booster_due  := pmax(booster_df$booster_due, 1)]
  
  # Increment counter for each already-received booster dose
  ppl[id %in% booster_df$id, vaccine_doses := vaccine_doses + booster_df$n_booster_doses]
  
  # ---- Initiate previously infected ----
  
  # Number of people infected, and IDs of those selected (uniform selection)
  n_prev_infected  = round(p$population_size * p$previously_infected)
  prev_infected_id = sort(sample.int(p$population_size, size = n_prev_infected))
  
  # Begin infections for these people
  ppl = fn_begin_infection(p, ppl, prev_infected_id, p$variants$id[1])
  
  # Sample day of infection from 1 (first possible day) to p$n_days_init (day 0 of simulation)
  day_infected = sample.int(n = -p$n_days_init, size = n_prev_infected,
                            prob = rev(p$seasonal_init), replace = TRUE)
  
  # Append temporary days_to_allocate variable to determine current state for those previously infected
  ppl[prev_infected_id, days_to_allocate := -p$n_days_init - day_infected + 1]
  
  # With everyone in the latent phase, begin to count number of days infected
  ppl[prev_infected_id, days_infected := pmin(days_next_event, days_to_allocate)]
  
  # Loop through phases of natual history until we get to recovery, death, or current state of infection
  while (sum(ppl$days_to_allocate, na.rm = TRUE) > 0) {
    
    # Number of days in 'next' stage, whatever that may be
    next_event = ppl[, days_next_event]
    
    # Take this period away from the durations we are allocating
    ppl[, days_next_event := pmax(0, days_next_event - days_to_allocate)]
    ppl[, days_to_allocate := pmax(0, days_to_allocate - next_event, na.rm = TRUE)]
    
    # Update disease and care state then sample duration of next phase
    list[ppl, , ] = update_state(p, ppl, 1)
    
    # Add duration of this phase to number of days infected and infectious
    ppl[prev_infected_id, days_infected   := days_infected   + pmin(days_next_event, days_to_allocate)]
    ppl[prev_infected_id, days_infectious := days_infectious + pmin(days_next_event, days_to_allocate)]
    
    # If we've got to recovery stage, allocate the remaining time in this stage
    ppl[days_recovered == 0, days_recovered := days_to_allocate]
    ppl[days_recovered >= 0, days_to_allocate := 0]
  }
  
  # Now remove temporary 'days_to_allocate' variable
  ppl[, days_to_allocate := NULL]
  
  # For ease, do not consider those already with symptoms for testing & diagnosis
  ppl[, test_date := NA]
  
  # ---- Only for those still currently infected ----
  
  # Quantify initial viral load - used to quantify infectiousness
  ppl = fn_viral_load(p, ppl)
  
  # ---- Only for those who have recovered ----
  
  # Variant prevalence is trivial at analysis start time - 100% 'primary variant'
  variant_prevalence = c(1, rep(0, nrow(p$variants) - 1))
  
  # Provide immunity for those recovered (noting they may have already been vaccinated)
  list[ppl, ] = fn_immunity_infection(p, ppl, variant_prevalence)
  
  # Sanity check for numerical vaccine_effect
  if (any(is.na(ppl$vaccine_effect)))
    stop("Non-numerical vaccine effect")
  
  # ---- Tidy up ----
  
  # Remove all attributes
  ppl = ppl[, lapply(.SD, as.vector)] 
  
  return(ppl)
}

# ---------------------------------------------------------
# Construct model output dataframe 
# ---------------------------------------------------------
preallocate_output = function(p) {
  
  # Preallocate output
  m = list()
  
  # Easy access metric info datatable
  metrics = p$metrics$df
  
  # Loop through all types of grouping
  for (this_grouping in unique(p$metrics$groupings)) {
    
    # Stratification of this grouping (can be altered by user)
    stratify_grouping = p$count[[this_grouping]]
    
    # Use NA to denote no grouping
    if (is.null(stratify_grouping))
      stratify_grouping = NA
    
    # All metrics for which we will provide this grouping
    metrics_idx = grepl(this_grouping, metrics$grouping)
    
    # These can be over time or not
    metrics_temporal = metrics$metric[metrics_idx & metrics$temporal]
    metrics_single   = metrics$metric[metrics_idx & !metrics$temporal]
    
    # Full factorial expansion of metrics that use this grouping over time
    m_id = paste0(this_grouping, "_temporal")
    m[[m_id]] = expand_grid(metric   = metrics_temporal, 
                            date     = 1 : p$n_days, 
                            group    = stratify_grouping, 
                            grouping = this_grouping)
    
    # Full factorial expansion of metrics that use this grouping (no time index)
    m_id = paste0(this_grouping, "_single")
    m[[m_id]] = expand_grid(metric   = metrics_single, 
                            date     = NA, 
                            group    = stratify_grouping, 
                            grouping = this_grouping)
  }
  
  # Bind into single dataframe
  m = rbindlist(m) %>%
    mutate(metric = factor(metric, levels = metrics$metric)) %>%
    arrange(metric, grouping, date, group)
  
  # Preallocate column for model ouput
  m$value = 0
  
  return(m)
}

# ---------------------------------------------------------
# Update model output datatable with counts of outcomes
# ---------------------------------------------------------
update_output = function(m, p, date_idx, update_metric, to_count, denom = 1) {
  
  # All groupings to disaggregate by for this metric
  metric_groupings = p$metrics$groupings[names(p$metrics$groupings) == update_metric]
  
  # Iterate through the groupings
  for (this_grouping in metric_groupings) {
    
    # For metrics with no possible disaggregation...
    if (this_grouping == "na") {
      
      # ... append single value
      metric_count = to_count / denom
      m[metric == update_metric & date == date_idx & is.na(group), value := metric_count]
      
      # For metrics we have chosen not to disaggregate...
    } else if (this_grouping == "none") {
      
      # ... count all values
      metric_count = to_count[, .N] / as.data.table(denom)[, .N]
      m[metric == update_metric & date == date_idx & is.na(group), value := metric_count]
      
    } else {  # Otherwise we want to count each grouping
      
      # Count number of occurrences of each group
      metric_count = group_count = table(to_count[, ..this_grouping])
      
      # Also count group occurrences in denominator (if provided)
      if (is.data.frame(denom)) {
        denom_count = table(denom[, ..this_grouping])
        
        # Divide denominator through (only for non-trivial groups)
        metric_count = group_count / denom_count[names(group_count)]
        
      } else { # If denominator is not a dataframe...
        
        # ... ensure it is nothing other than trivial
        if (!identical(denom, 1))
          stop("Denominator must be a dataframe if counting metrics by group")
      }
      
      # Insert these grouped values into output datatable
      m[metric == update_metric &          # This metric
          date == date_idx &               # This time step
          grouping == this_grouping &      # This grouping
          group %in% names(metric_count),  # This subset of 'groups'
        value := metric_count]
    }
  }
  
  # Check the update: m[metric == update_metric & date == date_idx, ]
  
  return(m)
}

# ---------------------------------------------------------
# Prepare final model output
# ---------------------------------------------------------
format_output = function(m, p, Re_info, network, ppl, yaml, seed) {
  
  # Insert NPI effect and seasonality profile
  #
  # NOTE: These are actually inputs, not outputs, but are helpful for visualisation
  m[metric == "contact_reduction", value := p$npi_contact_reduction]
  m[metric == "seasonality",       value := 100 * p$seasonality[1 : p$n_days]]
  
  # Convert number of cases per variant to variant prevalence
  m[metric == "variant_prevalence", value := 100 * value / max(sum(value), 1), by = "date"]
  
  # Convert total vaccine doses into cumulative measures
  m[metric == "total_doses", value := cumsum(value), by = .(grouping, group)]
  
  # Convert total treatments into cumulative measures
  m[metric == "total_treat",         value := cumsum(value), by = .(grouping, group)]
  m[metric == "total_treat_success", value := cumsum(value), by = .(grouping, group)]
  
  # Store effective reproduction number (calculated in postprocess.R)
  Re_info$Re_df$value[1 : (p$Re_fit$start_day - 1)] = NA  # Remove values before fitting period
  m[metric == "Re", value := Re_info$Re_df$value]
  
  # Count number of infections per person, capped by max_infection_count
  n_infections = table(ppl[, num_infections])
  m[metric == "n_infections" & group %in% names(n_infections), value := n_infections]
  
  # Append scenario and seed columns to output datatable
  m$scenario = p$.id
  m$seed = seed
  
  # Append ages of contacts to edge list
  network_age = network %>% 
    mutate(from_age = ppl[from, age],
           to_age   = ppl[to, age], .before = layer)
  
  # Combine inputs and outputs into one list
  results = list(yaml    = yaml$raw, 
                 input   = yaml$parsed, 
                 network = network_age,
                 Re_info = Re_info, 
                 output  = m)
  
  return(results)
}

# ---------------------------------------------------------
# Update disease and care states for all relevant individuals
# ---------------------------------------------------------
update_state = function(p, ppl, date_idx, m = NULL, cost_units = NULL) {
  
  # IDs of people that need to be updated
  update_id = ppl[days_next_event == 0, id]
  
  # Skip this process if no one to update
  if (length(update_id) > 0) {
    
    # Join with model flows dataframe to determine next state
    update_df = ppl[update_id, ] %>%
      left_join(y  = p$model_flows, 
                by = qc(disease_state, care_state, prognosis_state))
    
    # If m datatable is provided we'll want to count metrics
    if (!is.null(m)) {
      
      # Count and store values of each metric in model output
      for (update_metric in unique(na.omit(update_df$metric)))
        m = update_output(m, p, date_idx, update_metric, update_df[metric == update_metric, ])
    }
    
    # Store numbers entering hospital/ICU if required
    if (!is.null(cost_units)) {
      cost_units$hospital[date_idx] = update_df[metric == "hospital_admissions", .N]
      cost_units$icu[date_idx]      = update_df[metric == "icu_admissions", .N]
    }
    
    # Update disease and care state according to model_flows
    ppl[update_id, disease_state := update_df$next_disease]
    ppl[update_id, care_state    := update_df$next_care]
    
    # ---- Reset recovered and dead people ----
    
    # IDs of those no longer infected (recovery or death)
    reset_id = update_df[is.na(next_duration), id]
    if (length(reset_id) > 0) {
      
      # First things first, trivialise time to next event
      ppl[reset_id, days_next_event := NA]
      
      # Reset counters, viral load, variant, diagnosis and prognosis
      #
      # NOTE: Do not reset variant infected with - needed for immuno-escape considerations
      ppl[reset_id, days_infected    := NA]
      ppl[reset_id, days_infectious  := NA]
      ppl[reset_id, diagnosis_date   := NA]
      ppl[reset_id, treatment_date   := NA]
      ppl[reset_id, viral_load       := NA]
      ppl[reset_id, prognosis_state  := "none"]
      
      # Some extra things required for those who have recovered...
      
      # IDs of all who have recovered in this time step
      recover_id = update_df[metric == "recovered", id]
      
      # Signal start of recovery period by starting 'days recovered' count
      ppl[recover_id, days_recovered := 0]  # Incremented at end of daily loop
      
      # IDs of all who have died in this time step
      death_id = update_df[next_disease == "dead", id]
      
      # Check that we're not picking up anyone other than recovery or death
      if (!identical(reset_id, sort(c(recover_id, death_id))))
        stop("We should only 'reset' at recovery or death")
    }
    
    # Remove those we've reset from further updates
    update_df = update_df[!(id %in% reset_id)]
    
    # ---- Sample duration until next event ----
    
    # Loop through the different duration functions
    for (duration in unique(na.omit(update_df$next_duration))) {
      
      # Call duration function for the relevant individuals
      duration_id  = update_df[next_duration == duration, id]
      duration_val = p$durations[[duration]](duration_id)
      
      # Store the time until next event
      ppl[duration_id, days_next_event := duration_val]
    }
    
    # ---- Initaite infectious period and possible test & treat dates ----
    
    # Newly infectious people are those going INTO the pre-symptomatic phase
    new_infectious_id = update_df[next_disease == "presym", id]
    
    # Signal start of infectious period by starting 'days infectious' count
    ppl[new_infectious_id, days_infectious := 0]  # Incremented at end of daily loop
    
    # Newly symptomatic people are those going FROM the pre-symptomatic phase
    new_symptoms_id = update_df[disease_state == "presym", id]
    
    # Date of potential test and diagnosis - see fn_test_diagnose for probability
    ppl[new_symptoms_id, test_date := date_idx + p$diagnosis_delay]
  }
  
  return(list(ppl, m, cost_units))
}

# ---------------------------------------------------------
# Calculate healthcare and intervention costs
# ---------------------------------------------------------
health_economics = function(m, p, cost_units) {
  
  # Total cost of all interventions
  intervention_costs = 
    cost_units$vaccine   * p$intervention_costs$vaccine + 
    cost_units$prep      * p$intervention_costs$prep + 
    cost_units$treatment * p$intervention_costs$treatment
  
  # Total healthcare costs
  healthcare_costs = 
    cost_units$hospital * p$healthcare_costs$hospital + 
    cost_units$icu      * p$healthcare_costs$icu
  
  # Store these costs
  m[metric == "intervention_costs", value := intervention_costs]
  m[metric == "healthcare_costs",   value := healthcare_costs]
  
  # Also store overal costs (cost of interventions + healthcare costs)
  m[metric == "overall_costs", value := intervention_costs + healthcare_costs]
  
  return(m)
}

# ---------------------------------------------------------
# Calculate basic epidemiological indicators
# ---------------------------------------------------------
fn_epi_stats = function(m, p, ppl, epi_stats, date_idx) {
  
  # Total number of people currently infected and infectious
  all_infected   = ppl[days_infected > 0, id]
  all_infectious = ppl[days_infectious > 0, id]  # Including those in isolation
  
  # Easy access total number of people currently infected
  n_infected = length(all_infected)
  
  # Store total number of people infected and infectious
  m = update_output(m, p, date_idx, "currently_infected",   ppl[all_infected, ])
  m = update_output(m, p, date_idx, "currently_infectious", ppl[all_infectious, ])
  
  # Total number of people currently in isolation
  currently_isolated = ppl[care_state %in% "iso", id]
  m = update_output(m, p, date_idx, "currently_isolated", ppl[currently_isolated, ])
  
  # All infectious individuals not in isolation - these contributing to new infections
  epi_stats$all_infecting = setdiff(all_infectious, currently_isolated)
  
  # Total number of people currently infected with symptoms (mild, severe, or critical)
  currently_symptomatic = ppl[disease_state %in% c("mild", "severe", "crit"), id]
  m = update_output(m, p, date_idx, "currently_symptomatic", ppl[currently_symptomatic, ])
  
  # All susceptible people and their level of immunity (from previous infection and/or vaccination)
  epi_stats$all_susceptible      = ppl[disease_state == "susc", id]
  epi_stats$susceptible_immunity = ppl[epi_stats$all_susceptible, immune_state]
  
  # For population susceptibility, weight by immunity state and sum
  pop_susceptibility = 100 * sum(1 - epi_stats$susceptible_immunity) / p$population_size
  m = update_output(m, p, date_idx, "pop_susceptibility", pop_susceptibility)
  
  # Prevalence across whole population (proportion currently infected)
  pop_prevalence = 100 * n_infected / p$population_size
  m = update_output(m, p, date_idx, "pop_prevalence", pop_prevalence)
  
  # Sero-prevalence across whole population (proportion that have been infected thus far)
  pop_seroprevalence = 100 * ppl[num_infections > 0, .N] / p$population_size
  m = update_output(m, p, date_idx, "seroprevalence", pop_seroprevalence)
  
  # Virus variant of each infected individual
  all_variants = factor(ppl[all_infected, variant], levels = p$variants$id)
  m = update_output(m, p, date_idx, "variant_prevalence", ppl[all_infected, ])
  
  # Store variant prevalence as simple unnamed vector for use when importing cases
  if (n_infected > 0) {
    epi_stats$variant_prevalence = as.numeric(table(all_variants) / n_infected)
    
  } else {  # Deal with edge case of no new infections
    epi_stats$variant_prevalence = c(1, rep(0, nrow(p$variants) - 1))
  }
  
  return(list(m, epi_stats))
}

# ---------------------------------------------------------
# Update disease progression for the next time step
# ---------------------------------------------------------
fn_epi_update = function(m, p, ppl, date_idx, cost_units, states) {
  
  # Move one day closer to next 'event' (see model_flow config file for details)
  ppl[!is.na(days_next_event), days_next_event := days_next_event - 1]
  
  # Sanity check: throw an error if days to next event has become negative
  if (ppl[!is.na(days_next_event) & days_next_event < 0, .N] > 0)
    stop("Model error: negative days until next event")
  
  # Count number of people currently in hospital and ICU before updating states
  #
  # NOTE: Metrics for which a person is only counted once (eg admissions) are handled in update_state
  m = update_output(m, p, date_idx, "hospital_beds", ppl[care_state %in% c("hospital", "rehosp"), ])
  m = update_output(m, p, date_idx, "icu_beds",      ppl[care_state == "icu", ])
  
  # Update disease and care state then sample duration of next phase
  list[ppl, m, cost_units] = update_state(p, ppl, date_idx, m = m, cost_units = cost_units)
  
  # Increment the number of days infected and infectious, or recovered
  ppl[!is.na(days_infected),   days_infected   := days_infected   + 1]
  ppl[!is.na(days_infectious), days_infectious := days_infectious + 1]
  
  # Update viral load in next time step for all infected people
  ppl = fn_viral_load(p, ppl)
  
  # Also increment the number of days since recovery and/or vaccination
  ppl[!is.na(days_recovered),  days_recovered  := days_recovered  + 1]
  ppl[!is.na(days_vaccinated), days_vaccinated := days_vaccinated + 1]
  ppl[!is.na(days_booster),    days_booster    := days_booster    + 1]
  ppl[!is.na(days_prep),       days_prep       := days_prep       + 1]
  
  # Calculate number in each disease state
  for (state in p$model_states$disease)
    states[[state]][date_idx] = length(ppl[disease_state == state, id])
  
  # Calculate number in each care state
  for (state in p$model_states$care)
    states[[state]][date_idx] = length(ppl[care_state == state, id])
  
  return(list(m, ppl, cost_units, states))
}

# ---------------------------------------------------------
# Calculate new local infections
# ---------------------------------------------------------
fn_transmission = function(m, p, ppl, network, epi_stats, date_idx) {
  
  # ---- Sample contacts considering NPIs ----
  
  # NPI effect on at this time step
  npi_effect = p$npi_contact_reduction[date_idx]
  
  # Check that NPI is not trivial
  if (npi_effect > 1e-6) {
    
    # Number of effective contacts considering any NPIs in place
    n_contacts  = network[, .N]
    n_effective = n_contacts * (1 - npi_effect)
    
    # Sample indices of effective contacts 
    effective_idx = sample_int_crank(n    = n_contacts,
                                     size = n_effective,
                                     prob = rep(1, n_contacts))
    
    # All infected-susceptible 
    exposures = network[effective_idx][from %in% epi_stats$all_infecting & 
                                         to %in% epi_stats$all_susceptible]
    
  } else {  # No sampling needed if there is no NPI effect
    exposures = network[from %in% epi_stats$all_infecting & 
                          to %in% epi_stats$all_susceptible]
  }
  
  # ---- Per-contact transmission probability ----
  
  # Variant index of infectious individuals (can impact both infectiousness and susceptibility)
  variant_idx = match(ppl[exposures$from, variant], p$variants$id)
  
  # Infectiouness is a multiplier of viral load, seasonality, and viral variant
  beta = ppl[exposures$from, viral_load] * 
    p$beta * 
    p$seasonality[date_idx] * 
    p$air_pollution_factor$susceptibility * 
    p$variants$infectivity[variant_idx]
  
  # Calculate immunity for each exposure (depends on variant exposed to)
  list[ppl, exposures] = fn_immunity_infection(p, ppl, epi_stats$variant_prevalence, exposures)
  
  # ---- Sample random numbers to determine transmission ----
  
  # Draw random numbers to determine which contacts lead to transmission
  transmission_df = exposures %>% 
    mutate(probability = pmin(beta * (1 - immunity), 1),  # Probability of transmission in each contact
           rand = runif(n()), 
           transmission = rand < probability) %>% 
    filter(transmission == TRUE) %>%  # Only interested in transmission events
    group_by(to) %>%  # Uniqueness then guaranteed
    summarise(from = first(from)) %>%  # Only consider first transmission if multiple
    ungroup() %>%
    mutate(variant = ppl[from, variant]) %>%  # Variant is passed down from infector
    setDT()
  
  # NOTE: If desirable, we could here store number of people infected by each individual - could 
  #       be useful if looking at the effects on removing key individuals from transmission chain
  
  # Extract the IDs and variants of new local transmissions
  local_id    = transmission_df$to 
  infector_id = transmission_df$from
  
  # Begin infections for these people
  ppl = fn_begin_infection(p, ppl, local_id, transmission_df$variant)
  
  # Store number of new local infections among susceptibles
  m = update_output(m, p, date_idx, "new_local_infections", ppl[local_id, ])
  
  # Output IDs of newly infected as a list
  id = list(local = local_id, infector = infector_id)
  
  return(list(m, ppl, id))
}

# ---------------------------------------------------------
# Import infections acquired elsewhere
# ---------------------------------------------------------
fn_importation = function(m, p, ppl, epi_stats, id, date_idx) {
  
  # Number of infection imports this time step
  #
  # NOTE: Rounding up here as we always want at least 1 and sample_vec rounds down by default
  n_imports = ceiling(p$import_constant * length(epi_stats$all_susceptible))
  
  # Sample indices - random for now, but may want to consider other factors later
  id$import = sample_vec(x    = epi_stats$all_susceptible, 
                         size = n_imports, 
                         prob = 1 - epi_stats$susceptible_immunity)
  
  # Select virus variants for these individuals based on the proportion currently in circulation
  apply_variant = sample_vec(x    = p$variants$id, 
                             size = n_imports, 
                             prob = epi_stats$variant_prevalence, 
                             replace = TRUE)
  
  # Begin infections for these people
  ppl = fn_begin_infection(p, ppl, id$import, apply_variant)
  
  # Store number of new imported infections
  m = update_output(m, p, date_idx, "new_importations", ppl[id$import, ])
  
  return(list(m, ppl, id))
}

# ---------------------------------------------------------
# New viral variant emergence
# ---------------------------------------------------------
fn_emergence = function(m, p, ppl, epi_stats, id, date_idx) {
  
  # All newly infected people (including imports) we could use to initiate new mutation
  newly_infected_id = c(id$local, id$import)
  
  # Check whether any variants are due to be imported in this time step
  if (date_idx %in% p$variants$import_day) {
    
    # Details of the variant to be imported
    this_variant = p$variants[p$variants$import_day == date_idx, ]
    
    # Do we have enough new infections to assign the new variant to
    n_imports_remain = max(0, this_variant$import_number - length(newly_infected_id))
    
    # Deal with case where we need to introduce elsewhere
    if (n_imports_remain > 0) {
      
      # Can reassign variant for the already infected OR force infect susceptibles - we do the latter
      newly_forced_id = sample_vec(x    = epi_stats$all_susceptible, 
                                   size = n_imports_remain, 
                                   prob = 1 - susceptible_immunity)
      
      # Begin infections for these people
      ppl = fn_begin_infection(p, ppl, newly_forced_id, this_variant$id)
      
      # Extend the people we can choose from to 'reassign' the new variant to
      newly_infected_id = c(newly_infected_id, newly_forced_id)
    }
    
    # Sample this many people from those recently infected 
    id$new_variant = sample_vec(newly_infected_id, this_variant$import_number)
    
    # Assign the new mutation variant to these individuals
    ppl[id$new_variant, variant := this_variant$id]
    
  } else { # No viral emergence on this day...
    
    # Output trivial IDs
    id$new_variant = integer(0)
  }
  
  return(list(m, ppl, id))
}

# ---------------------------------------------------------
# Count incidence and serial intervals for all new infections
# ---------------------------------------------------------
fn_incidence = function(m, p, ppl, epi_stats, id, date_idx) {
  
  # All newly infected people (including all imports)
  all_new_infections = c(id$local, id$import, id$new_variant)
  new_infections     = c(id$local, id$import)
  
  # Store total number of new infections
  m = update_output(m, p, date_idx, "all_new_infections", ppl[all_new_infections, ])
  
  # Store incidence in this time step (including imported cases)
  epi_stats$incidence[date_idx] = length(new_infections)
  
  # IDs of those who infected others in this time step
  infector_id = id$infector
  
  # Calculate serial interval for each of these infector individuals
  epi_stats$serial_int[[date_idx]] = ppl[id %in% infector_id, days_infected]
  
  return(list(m, epi_stats))
}

# ---------------------------------------------------------
# Initiate new infections for newly infected individuals
# ---------------------------------------------------------
fn_begin_infection = function(p, ppl, id, apply_variant) {
  
  # Sanity check: we need a variant defined for each person
  if (!(length(apply_variant) %in% c(length(id), 1)))
    stop("Inconsistent number of variants defined")
  
  # Apply the precalculated variant for this individial
  ppl[id, variant := apply_variant]
  
  # Sample a prognosis for each infected individual
  ppl[id, prognosis_state := fn_prognosis(p, ppl[id, ])]
  
  # All newly infected people go to the latent phase
  ppl[id, disease_state := "latent"]
  
  # Define how long this latency phase will last
  ppl[id, days_next_event := pmax(1, round(p$duration$latency(id)))]
  
  # We signal this new infection by starting 'days infected' count
  ppl[id, days_infected := 0]  # Incremented at end of daily loop
  
  # Increment number of infections experienced thus far
  # 
  # NOTE: We do this after sampling prognosis to not count the current
  #       infection in the exposure-response disease severity effect
  ppl[id, num_infections := pmin(num_infections + 1L, p$max_infection_count)]
  
  # Reset counter for days recovered if necessary
  ppl[id, days_recovered := NA]
  
  return(ppl)
}

# ---------------------------------------------------------
# Generate prognosis state for newly infected individuals
# ---------------------------------------------------------
fn_prognosis = function(p, new_infected) {
  
  # Disease blocking effects of vaccination and/or PrEP
  disease_blocking_df = fn_immunity_disease(p, new_infected)
  
  # Per-person severity factors considering all influencing variables
  severity_df = new_infected %>% 
    select(id, age, variant, vaccine_doses, num_infections) %>%
    left_join(y  = disease_blocking_df,  # Vaccine and PrEP effect
              by = "id") %>%
    left_join(y  = p$prognosis$severity,  # Variant and dose/exposure effect
              by = c("variant", "vaccine_doses", "num_infections")) %>%
    mutate(value = value * p$air_pollution_factor$severity,   # Air pollution effect
           value = value * (1 - disease_blocking)) %>%
    select(id, state, value)
  
  # Mulitply severity factors by baseline age-related prognosis probabilities
  prognosis_df = new_infected %>% 
    select(id, age, comorbidities) %>% 
    mutate(age = pmin(age + comorbidities * 10, max(p$ages))) %>%  # Comorbidity effect
    left_join(y  = p$prognosis$age,  # Base age-related prognosis
              by = "age") %>%
    select(id, state, value) %>%
    bind_rows(severity_df) %>%
    group_by(id, state) %>%
    summarise(value = prod(value)) %>%  # Multiply all factors through
    mutate(state = factor(state, levels = p$model_states$prognosis)) %>%
    arrange(id, state) %>%
    mutate(value = cumsum(value), 
           value = value / max(value)) %>%  # Normalise to 100% max
    ungroup() %>% 
    setDT()
  
  # Generate random number for every individual
  rand_df = new_infected[, .(id, rand = runif(.N))]
  
  # Use this random number to select a prognosis for each person
  prognosis = prognosis_df %>%
    left_join(y  = rand_df, 
              by = "id") %>%
    filter(value > rand) %>%
    group_by(id) %>%
    slice_head(n = 1) %>%
    pull(state)
  
  return(prognosis)
}

# ---------------------------------------------------------
# Which infected people will get tested AND diagnosed, and when
# ---------------------------------------------------------
fn_test_diagnose = function(m, p, ppl, date_idx) {
  
  # ---- Regular testing ----
  
  # All people who could potentially get tested today
  #
  # NOTE: Occasionally people recover before getting tested, hence days_infectious check
  test_potential   = ppl[test_date == date_idx & !is.na(days_infectious), ]
  n_test_potential = test_potential[, .N]
  
  # Skip this process if no potential testers identified
  if (n_test_potential > 0) {
    
    # Preallocate vector for probability of test AND diagnosis
    dx_probability = rep(1, n_test_potential)
    
    # Indices of all those who will have mild or no symptoms
    asym_idx = test_potential$prognosis_state == "asym"
    mild_idx = test_potential$prognosis_state == "mild"
    
    # Age probabilty vectors of test AND diagnosis in these two cases
    #
    # NOTE: These probabilities already include test sensitivity
    asym_age_prob = p$testing$without_symptoms$age_prob
    mild_age_prob = p$testing$with_symptoms$age_prob
    
    # Apply these probabilities based on age
    dx_probability[asym_idx] = asym_age_prob[test_potential[asym_idx, age] + 1]
    dx_probability[mild_idx] = mild_age_prob[test_potential[mild_idx, age] + 1]
    
    # Generate random number vector
    random_number = runif(n_test_potential)
    
    # ID of all those to be tested AND diagnosed
    dx_id = test_potential[which(random_number < dx_probability), id]
    
    # Assign the date that this diagnosis is to occur
    ppl[dx_id, diagnosis_date := date_idx]
  }
  
  # ---- Mass testing ----
  
  # Easy access mass testing details
  mass = p$testing$mass_testing
  
  # Are we due to mass test today?
  if (date_idx %in% p$testing$mass_testing$when) {
    
    # Number of people to mass test on this date
    #
    # NOTE: This is essentially 'effective coverage' as it considers test sensitivity
    n_mass_test = round(p$population_size * mass$coverage * mass$sensitivity)
    
    # IDs of those to be mass tested today
    mass_test_id = sample_vec(ppl[mass_test_accept == TRUE, id], size = n_mass_test)
    
    # All people to diagnose through mass testing today
    mass_diagnose_id = ppl[id %in% mass_test_id & viral_load > 0, id]
    
    # Remove any IDs in mass_test_id that have already been diagnosed
    mass_diagnose_id = setdiff(mass_diagnose_id, ppl[diagnosis_date <= date_idx, id])
    
    # Assign diagnosis date and remove future test date
    ppl[mass_diagnose_id, diagnosis_date := date_idx]
    ppl[mass_diagnose_id, test_date := NA]
  }
  
  # Remove reference to all expired test dates
  ppl[test_date <= date_idx, test_date := NA]
  
  # Store total number of those diagnosed this time step
  m = update_output(m, p, date_idx, "confirmed", ppl[diagnosis_date == date_idx, ])
  
  return(list(m, ppl))
}

# ---------------------------------------------------------
# Isolate newly diagnosed individuals where appropriate
# ---------------------------------------------------------
fn_isolate = function(m, p, ppl, date_idx) {
  
  # Once diagnosed, only those with certain prognoses can go into isolation (see model_flows)
  #
  # INTERPRETATION: You would not isolate if you have severe disease
  can_isolate_id = ppl[diagnosis_date == date_idx & 
                         prognosis_state %in% p$model_states$iso, id]
  
  # Sample number of people due to isolate (according to isolation_probability)
  n_isolate  = length(can_isolate_id) * p$isolation$probability
  isolate_id = sample_vec(can_isolate_id, round(n_isolate))
  
  # Mark these people as isolated and keep them there for isolation_duration days
  ppl[isolate_id, care_state := "iso"]
  ppl[isolate_id, days_next_event := round(p$isolation$duration)]
  
  return(list(m, ppl))
}

# ---------------------------------------------------------
# Vaccinate eligible people
# ---------------------------------------------------------
fn_vaccinate = function(m, p, ppl, date_idx, n_priority_group, cost_units) {
  
  # Extract vaccine type currently being distributed
  current_vaccine = p$vaccine_update %>%
    filter(release_day <= date_idx) %>%
    slice_max(release_day) %>%
    pull(vaccine_type)
  
  # Coverage difference between yesterday and today
  #
  # NOTE: Binding initial coverage so we consistently index a 2D matrix
  coverage_matrix = rbind(p$vaccine$coverage_init, p$vaccine$coverage)
  coverage_diff   = coverage_matrix[date_idx : (date_idx + 1), ] %>%
    matrix(nrow = 2, dimnames = list(1 : 2, colnames(coverage_matrix)))
  
  # Iterate through the vaccine priority groups
  for (group in p$priority_groups$id) {
    
    # Only continue if we have people to vaccinate at this time step
    if (diff(coverage_diff[, group]) > 0) {
      
      # All people eligible for vaccination in this group at this time step 
      eligible = ppl[priority_group == group & vaccine_accept == TRUE & is.na(days_vaccinated), ]
      
      # Number of people we want vaccinated in this group in this time step
      n_group   = n_priority_group[priority_group == group, n]
      n_desired = round(p$vaccine$coverage[date_idx, group] * n_group)
      
      # Number of people currently vaccinated in this group
      n_current = ppl[priority_group == group & !is.na(days_vaccinated), .N]
      
      # How many we will vacinate this time step - bounded above by number available
      n_vaccinate = min(max(n_desired - n_current, 0), nrow(eligible))
      
      # Vaccinate the first chunk of people - randomness is inherent as table is not ordered
      vaccinate_id = eligible[seq_len(n_vaccinate), id]
      
      # Start these people off on their vaccination journey
      ppl[vaccinate_id, days_vaccinated := 0]
      ppl[vaccinate_id, vaccine_type := current_vaccine]
    }
  }
  
  # ---- Apply booster dose ----
  
  # IDs of those to receive a booster today (limited to booster_doses per day)
  booster_id = ppl %>%
    filter(booster_due <= date_idx) %>%
    left_join(y  = n_priority_group, 
              by = "priority_group") %>%
    select(id, booster_due, priority) %>%
    arrange(booster_due, priority) %>%
    slice_head(n = p$booster_doses) %>%
    pull(id)
  
  # For all receiving a booster, initiate days since most recent dose received
  ppl[booster_id, days_booster := 0]
  ppl[booster_id, vaccine_type := current_vaccine]
  
  # ---- Due date of next booster dose ----
  
  # NOTE: Booster cycle initiates on day of final dose UNLESS using force_start, in which
  #       case first booster given on day force_start. Cycle then continues from that point. 
  
  # Day of final dose of initial vaccine schedule
  final_dose = max(p$vaccine$subsequent_dose_days)
  
  # Set date of first booster dose (for those recieving final dose of initial schedule today)
  booster_1_df = ppl %>%
    filter(days_vaccinated == final_dose, 
           booster_accept == TRUE) %>% 
    select(id, priority_group, days_vaccinated, days_booster, booster_due) %>% 
    left_join(y  = p$booster_details[, -"probability"], 
              by = c("priority_group" = "id")) %>%
    mutate(booster_due = pmin(pmax(cycle_period + date_idx, start), 
                              force_start, na.rm = TRUE), 
           booster_due = ifelse(booster_due <= force_end, booster_due, NA))
  
  # Set date of subsequent booster dose (for those receiving a booster today)
  booster_n_df = ppl %>%
    filter(id %in% booster_id) %>% 
    select(id, priority_group, days_vaccinated, days_booster, booster_due) %>% 
    left_join(p$booster_details[, -"probability"], 
              by = c("priority_group" = "id")) %>%
    mutate(booster_due = booster_due + cycle_period, 
           booster_due = ifelse(booster_due <= force_end, booster_due, NA))
  
  # Apply next booster due date calculated above
  ppl[id %in% booster_1_df$id, booster_due := booster_1_df$booster_due]
  ppl[id %in% booster_n_df$id, booster_due := booster_n_df$booster_due]
  
  # ---- Increment number of doses per person ----
  
  # Primary vaccine doses
  ppl[days_vaccinated %in% c(0, p$vaccine$subsequent_dose_days), 
      vaccine_doses := vaccine_doses + 1L]
  
  # Booster doses
  ppl[days_booster == 0, vaccine_doses := vaccine_doses + 1L]
  
  # Stop counting when we get up to max_dose_count
  ppl[, vaccine_doses := pmin(vaccine_doses, p$max_dose_count)]
  
  # Handle all vaccine-related metrics (there are a few!) in seperate function for readability
  list[m, cost_units] = fn_vaccine_count(m, p, ppl, date_idx, cost_units)
  
  return(list(m, ppl, cost_units))
}

# ---------------------------------------------------------
# Handle all vaccine and booster metrics
# ---------------------------------------------------------
fn_vaccine_count = function(m, p, ppl, date_idx, cost_units) {
  
  # Number vaccinated total and historically
  vaccine_today = ppl[days_vaccinated == 0, ]
  vaccine_total = ppl[days_vaccinated >= 0, ]
  
  # Store the number of people vaccinated today and in total
  m = update_output(m, p, date_idx, "n_vaccinated",     vaccine_today)
  m = update_output(m, p, date_idx, "total_vaccinated", vaccine_total)
  
  # People receiving doses: first dose, subsequent doses, and booster doses
  vaccine_dose = ppl[days_vaccinated %in% c(0, p$vaccine$subsequent_dose_days), ]
  booster_dose = ppl[days_booster == 0, ]
  
  # Bind these dataframes together
  all_doses = rbind(vaccine_dose, booster_dose)
  
  # Store number of vaccine doses administered in this time step
  m = update_output(m, p, date_idx, "n_doses",     all_doses)
  m = update_output(m, p, date_idx, "total_doses", all_doses)  # Cum-summed in format_output
  
  # Total number of people eligible for vaccination (used as denominator)
  vaccine_elig = ppl[priority_group != "none", ]
  
  # Vaccine coverage, using both those eligible and total pop as the denominator
  m = update_output(m, p, date_idx, "vaccine_coverage",     vaccine_total, denom = vaccine_elig)
  m = update_output(m, p, date_idx, "vaccine_coverage_pop", vaccine_total, denom = ppl)
  
  # Any number of boosters recieved at any time
  booster_total = ppl[days_booster >= 0, ]
  
  # Booster coverage among those eligible (ie previously vaccinated) and population level
  m = update_output(m, p, date_idx, "booster_coverage",     booster_total, denom = vaccine_total)
  m = update_output(m, p, date_idx, "booster_coverage_pop", booster_total, denom = ppl)
  
  # Boosted in past 365 days
  booster_12m = ppl[days_booster >= 0 & days_booster <= 365, ]
  
  # Coverage of boosters in past 365 days
  m = update_output(m, p, date_idx, "booster_coverage_12m",     booster_12m, denom = vaccine_total)
  m = update_output(m, p, date_idx, "booster_coverage_12m_pop", booster_12m, denom = ppl)
  
  # Update number of vaccine doses in intervention unit counter
  cost_units$vaccine[date_idx] = all_doses[, .N]
  
  return(list(m, cost_units))
}

# ---------------------------------------------------------
# Pre-exposure prophylaxis for those eligible
# ---------------------------------------------------------
fn_prep = function(m, p, ppl, date_idx, cost_units) {
  
  # ---- Distribute PrEP as appropriate ----
  
  # All those eligible that have and haven't yet received PrEP
  prep_yes = ppl[vax_unsuitable == TRUE & !is.na(days_prep), ]
  prep_no  = ppl[vax_unsuitable == TRUE & is.na(days_prep), ]
  
  # Total number of people eligible
  n_total = prep_yes[, .N] + prep_no[, .N]
  
  # Number of people we want to have recieved PrEP in this time step
  n_desired = round(p$prep$coverage[date_idx] * n_total)
  
  # How many we will provide PrEP to this time step - bounded above by number available
  n_prep = min(max(n_desired - prep_yes[, .N], 0), prep_no[, .N])
  
  # Randomly select eligible people to provide PrEP to
  prep_id = sample_vec(x = prep_no[, id], size = n_prep)
  
  # Start these people off on their PrEP journey
  ppl[prep_id, days_prep := 0]
  
  # ---- Store details ----
  
  # Number vaccinated total and historically
  prep_today = ppl[days_prep == 0, ]
  prep_total = ppl[days_prep >= 0, ]
  
  # Store the number of people who have received PrEP (today and in total)
  m = update_output(m, p, date_idx, "n_prep",     prep_today)
  m = update_output(m, p, date_idx, "total_prep", prep_total)
  
  # Total number of people eligible for PrEP (used as denominator)
  prep_elig = ppl[vax_unsuitable == TRUE, ]
  
  # PrEP coverage, using both those eligible and total pop as the denominator
  m = update_output(m, p, date_idx, "prep_coverage",     prep_total, denom = prep_elig)
  m = update_output(m, p, date_idx, "prep_coverage_pop", prep_total, denom = ppl)
  
  # Update number of vaccine doses in intervention unit counter
  cost_units$prep[date_idx] = prep_today[, .N]
  
  return(list(m, ppl, cost_units))
}

# ---------------------------------------------------------
# Treat eligible infected cases
# ---------------------------------------------------------
fn_treatment = function(m, p, ppl, date_idx, cost_units) {
  
  # ---- Treat infected individuals as appropriate ----
  
  # Extract people who could be due treatment today
  treat_elig = ppl %>%
    left_join(y  = p$treat_prob, 
              by = c("priority_group", "disease_state")) %>%
    left_join(y  = p$treat_available, 
              by = "priority_group") %>%
    filter(treat_prob > 0, 
           available <= date_idx,
           diagnosis_date == date_idx - p$treat_delay)
  
  # It remains to check their vaccine status and treatment conditions...
  #
  # NOTE: There is probably a better way of doing this in one step
  
  # Treatments for the unvaccinated, vaccinated, and regardless of vaccine status
  treat_unvax = treat_elig[vaccine_condition == "unvaccinated" & is.na(days_vaccinated), ]
  treat_vax   = treat_elig[vaccine_condition == "vaccinated"   & days_vaccinated > 0, ]
  treat_all   = treat_elig[vaccine_condition == "none", ]
  
  # Probability of treatment and of treatment success
  treat_df = rbind(treat_unvax, treat_vax, treat_all) %>%
    mutate(treat_success = treat_prob * p$treat_efficacy, 
           rand = runif(n()))
  
  # Those with successful treatment, and total treatement
  success_id = treat_df[rand < treat_success, id]
  total_id   = treat_df[rand < treat_prob, id]
  
  # Assign current date index to treatment_date for only those successful
  ppl[success_id, treatment_date := date_idx]
  
  # Reset prognosis, disease, and care values to represent a mild case...
  ppl[success_id, disease_state   := "mild"]
  ppl[success_id, care_state      := "none"]
  ppl[success_id, prognosis_state := "mild"]
  
  # ... which is just about to recover
  ppl[success_id, days_next_event := 1]
  
  # ---- Store details ----
  
  # Store total number of people treated this time step
  m = update_output(m, p, date_idx, "n_treat",     ppl[total_id, ])
  m = update_output(m, p, date_idx, "total_treat", ppl[total_id, ]) # Cum-summed in format_output
  
  # Also store the subset that had successful treatment
  m = update_output(m, p, date_idx, "n_treat_success",     ppl[success_id, ])
  m = update_output(m, p, date_idx, "total_treat_success", ppl[success_id, ]) # Cum-summed in format_output
  
  # Update number of vaccine doses in intervention unit counter
  cost_units$treatment[date_idx] = length(total_id)
  
  return(list(m, ppl, cost_units))
}

# ---------------------------------------------------------
# Update viral load for all infected individuals
# ---------------------------------------------------------
fn_viral_load = function(p, ppl) {
  
  # NOTE: We assume viral load in latent period is zero
  
  # ID of all infectious people
  infectious_id = ppl[days_infectious > 0, id]
  
  # Days since these people became infectious
  days_infectious = ppl[infectious_id, days_infectious]
  
  # Index with days since infection
  viral_load_values = p$viral_load_profile[days_infectious]
  
  # If 'infectious' longer than viral load profile, assume zero infectiousness
  viral_load_values[is.na(viral_load_values)] = 0
  
  # Update viral load variable in ppl dataframe
  ppl[infectious_id, viral_load := viral_load_values]
  
  return(ppl)
}

# ---------------------------------------------------------
# Immunity to infection when exposed
# ---------------------------------------------------------
fn_immunity_infection = function(p, ppl, variant_prevalence, exposures = NULL) {
  
  # ---- Immunity per exposure ----
  
  # First job is to calculate current 'effect' of vaccination - depends on day of vaccination and any booster
  ppl[days_vaccinated > 0, vaccine_effect := pmax(p$booster_profile[days_booster], 
                                                  p$vaccine$profile[days_vaccinated], na.rm = TRUE)]
  
  # Skip this step if no exposure details provided
  if (!is.null(exposures)) {
    
    # For each exposure, append variant details and other factors influencing immunity
    immunity_effect_df = exposures %>%
      left_join(y  = ppl[, .(id, from_variant = variant)], 
                by = c("from" = "id")) %>%
      left_join(y  = ppl[, .(id, to_variant = variant, vaccine_type, vaccine_effect, days_recovered)], 
                by = c("to" = "id")) %>%
      mutate(recovery_effect = p$acquired_immunity[pmax(days_recovered, 1)]) %>%
      select(-days_recovered)
    
    # Effect on vaccine-induced and acquired immunity - depends of previous exposure and vaccine type
    immunity_df = immunity_effect_df %>%
      left_join(y  = p$vaccine_update[, .(vaccine_type, target_variant)], 
                by = "vaccine_type") %>%
      mutate(variant_str = paste0("\\*", from_variant, "\\*"), 
             target_of_vaccine = mapply(grepl, variant_str, target_variant)) %>%
      select(-vaccine_type, -target_variant) %>%
      mutate(escape = 1 - p$variants$immuno_escape[match(from_variant, p$variants$id)], 
             escape_vaccine  = ifelse(target_of_vaccine == TRUE,  1, escape), 
             escape_acquired = ifelse(from_variant == to_variant, 1, escape)) %>%
      select(-escape, -target_of_vaccine)
    
    # TODO: We may want to model hybrid immunity to infection
    
    # Now we can compute overall level of immunity for each exposure - the max of both types
    exposures = immunity_df %>%
      mutate(immunity = pmax(escape_vaccine * vaccine_effect * p$vaccine$infection_blocking, 
                             escape_acquired * recovery_effect, na.rm = TRUE)) %>%
      select(from, to, immunity)
  }
  
  # ---- Population-level immunity ----
  
  # NOTE: The population-level calculation is more crude and uses prevalence of each variant in circulation
  
  # Population-level effect of immuno-escaping variants (prevalence * immuno-escape level)
  pop_escape = variant_prevalence * p$variants$immuno_escape
  
  # For all people, calculate effect of all immuno-escape variants on acquired immunity
  immunity_acquired = ppl[, .(id, variant, days_recovered)] %>%
    mutate(escape = 1 - (sum(pop_escape) - pop_escape[match(variant, p$variants$id)]), 
           value  = escape * p$acquired_immunity[pmax(days_recovered, 1)]) %>%
    pull(value) 
  
  # Then take the biggest immunity effect for all people - through infection or vaccination
  ppl[, immune_state := pmax(vaccine_effect * p$vaccine$infection_blocking * (1 - sum(pop_escape)),
                             immunity_acquired, na.rm = TRUE)]
  
  return(list(ppl, exposures))
}

# ---------------------------------------------------------
# Immunity to severe disease once infected
# ---------------------------------------------------------
fn_immunity_disease = function(p, new_infected) {
  
  # Vaccine effect can be reduced if variant infected with is immuno-escaping
  vaccine_effect = new_infected %>%
    select(id, variant, vaccine_effect, vaccine_type) %>%
    left_join(y  = p$vaccine_update[, .(vaccine_type, target_variant)], 
              by = "vaccine_type") %>%
    mutate(variant_str = paste0("\\*", variant, "\\*"), 
           target_of_vaccine = mapply(grepl, variant_str, target_variant)) %>%
    mutate(escape = 1 - p$variants$immuno_escape[match(variant, p$variants$id)], 
           escape = ifelse(target_of_vaccine == TRUE, 1, escape),
           value  = vaccine_effect * escape * (1 - p$vaccine$infection_blocking)) %>%
    select(id, disease_blocking_vaccine = value)
  
  # PrEP effect (considering how many days since injections)
  prep_effect  = new_infected %>%
    mutate(value = p$prep$profile[days_prep], 
           value = pmax(value, 0, na.rm = TRUE)) %>%
    select(id, disease_blocking_prep = value)
  
  # Combine the effects of vaccination and PrEP for overall disease blocking effect
  disease_blocking_df = vaccine_effect %>%
    full_join(y  = prep_effect, 
              by = "id") %>%
    mutate(value = 1 - (1 - disease_blocking_vaccine) * 
             (1 - disease_blocking_prep)) %>%
    select(id, disease_blocking = value)
  
  return(disease_blocking_df)
}

# ---------------------------------------------------------
# Non-COVID-related 'natural' deaths
# ---------------------------------------------------------
fn_natural_death = function(m, p, ppl, date_idx) {
  
  # Check natural death flag
  if (p$natural_deaths == TRUE) {
    
    # IDs of all those that die naturally, regardless of COVID status
    natural_death_id = ppl[, .(id, age)] %>%
      left_join(y  = p$natural_death_age, 
                by = "age") %>%
      mutate(rand  = runif(n()), 
             death = rand < probability) %>% 
      filter(death == TRUE) %>%
      pull(id)
    
    # Assign these individuals an informative temporary disease state
    #
    # NOTE: These individuals will be replaced with newborns regardless of p$replace_deceased flag
    ppl[natural_death_id, disease_state := "natural_dead"]
    
    # Store details of natural deaths in model output
    natural_deaths = ppl[natural_death_id, ]
    m = update_output(m, p, date_idx, "natural_deaths", natural_deaths)
  }
  
  return(list(m, ppl))
}

# ---------------------------------------------------------
# Add newborns to replace deceased
# ---------------------------------------------------------
fn_birth = function(m, p, ppl, date_idx) {
  
  # TODO: On initial call we may want to do better than simply assigning age zero
  #       A possible improvement could be to you p$n_days_init to sample 'newborn' age
  
  # IDs of all recently deceased
  birth_id = ppl[disease_state %in% p$replace_state, id]
  
  # Skip this process of no newborns to create
  if (length(birth_id) > 0) {
    
    # Create susceptible newborns (age 0, randomly assigned birthday index)
    ppl_newborn = create_ppl(p, n = length(birth_id), init = FALSE)
    
    # Overwrite IDs (1-n by default)
    ppl_newborn[, id := birth_id]
    
    # Replace the recently deceased with these newborns
    ppl[birth_id, (names(ppl)) := ppl_newborn]
  }
  
  # Store in number of births in model output - skip first time step
  if (date_idx > 1)
    m = update_output(m, p, date_idx, "births", length(birth_id))
  
  return(list(m, ppl))
}

# ---------------------------------------------------------
# Age individuals annually
# ---------------------------------------------------------
fn_ageing = function(m, p, ppl, date_idx) {
  
  # TODO: Update priority group status as people age into new group
  
  # Check ageing flag
  if (p$annual_ageing == TRUE) {
    
    # Take the modulo to determine day of the year
    day_of_year = date_idx %% 365  # Day index 0 used instead of 365
    
    # Index of individuals with this birthday
    age_idx = which(ppl$birthday == day_of_year)
    
    # Age these individuals by one year (capped by max age)
    ppl[age_idx, age := pmin(age + 1, max(p$ages))]
  }
  
  # ---- Count outcomes ----
  
  # Average population age (of those surviving)
  pop_age = mean(ppl[disease_state != "dead", age])
  
  # Store in model output
  m = update_output(m, p, date_idx, "pop_age", pop_age)
  
  return(list(m, ppl))
}

