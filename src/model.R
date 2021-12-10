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
model = function(o, scenario, seed = NA, fit = NULL, do_plot = FALSE, verbose = "date") {
  
  # Set random seed generator if input defined
  if (!is.na(seed)) set.seed(seed)
  
  # ---- Generate model parameters ----
  
  if (verbose != "none") message(" - Parsing input")
  
  # Load and parse user-defined inputs for this scenario
  yaml = parse_yaml(o, scenario = scenario, read_array = TRUE)
  
  # Apply fitted parameters - also used during calibration process
  yaml = fit_yaml(o, yaml, fit)
  
  # Short hand for model parameters
  p = yaml$parsed
  
  # ---- Model set up ----
  
  if (verbose != "none") message(" - Running model")
  
  # Initiate people datatable
  ppl_df = initiate_ppl(o, p, verbose)
  
  # Generate full contact network outside of main model loop (see networks.R)
  network = create_network(o, p, ppl_df, do_plot, verbose)
  
  # Age-related correction factors from network generation
  p$age_correction = network$age_correction
  p$age_contacts   = network$age_contacts
  
  # Initialise infections/vaccinations
  ppl_df = initiate_epidemic(o, p, ppl_df)
  
  # Calculate number of people in each vaccine prioirty group outside of main loop
  n_vaccine_group = ppl_df[, .N, by = vaccine_group] %>% 
    left_join(p$vaccine$groups[, .(id, priority)], 
              by = c("vaccine_group" = "id")) %>% 
    arrange(priority)
  
  # Group certain states for easy access and readability
  states = list(infected   = qc(latent, presym, asym, mild, severe, crit), 
                infectious = qc(presym, asym, mild, severe, crit), 
                in_care    = qc(hospital, icu, rehosp))
  
  # Also append vectors for all disease and care states
  for (state in p$model_states$all)
    states[[state]] = rep(0, p$n_days)
  
  # Preallocate vector needed to calculate effective reproduction number
  infectious_time = rep(0, p$n_days)
  
  # Preallocate model output datatable (denoted 'm' for short)
  m = preallocate_output(p)
  
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
    
    # ---- Basic epidemiological indicators ----
    
    # Total number of people currently infected and infectious
    all_infected   = ppl[days_infected > 0, id]
    all_infectious = ppl[days_infectious > 0, id]  # Including those in isolation
    
    # Easy access total number of people currently infected
    n_infected = length(all_infected)
    
    # Store total number of people infected and infectious
    m = update_output(m, p, i, "currently_infected",   ppl[all_infected, ])
    m = update_output(m, p, i, "currently_infectious", ppl[all_infectious, ])
    
    # Total number of people currently in isolation
    currently_isolated = ppl[care_state %in% "iso", id]
    m = update_output(m, p, i, "currently_isolated", ppl[currently_isolated, ])
    
    # Total number of people currently infected with symptoms (mild, severe, or critical)
    currently_symptomatic = ppl[disease_state %in% c("mild", "severe", "crit"), id]
    m = update_output(m, p, i, "currently_symptomatic", ppl[currently_symptomatic, ])
    
    # All susceptible people and their level of immunity (from previous infection and/or vaccination)
    all_susceptible      = ppl[disease_state == "susc", id]
    susceptible_immunity = ppl[all_susceptible, immune_state]
    
    # For population susceptibility, weight by immunity state and sum
    pop_susceptibility = 100 * sum(1 - susceptible_immunity) / p$population_size
    m = update_output(m, p, i, "pop_susceptibility", pop_susceptibility)
    
    # Prevalence across whole population (proportion currently infected)
    pop_prevalence = 100 * n_infected / p$population_size
    m = update_output(m, p, i, "pop_prevalence", pop_prevalence)
    
    # Sero-prevalence across whole population (proportion that have been infected thus far)
    pop_seroprevalence = 100 * nrow(ppl[n_infections > 0, ]) / p$population_size
    m = update_output(m, p, i, "seroprevalence", pop_seroprevalence)
    
    # Virus variant of each infected individual
    all_variants = factor(ppl[all_infected, variant], levels = p$variants$id)
    m = update_output(m, p, i, "variant_prevalence", ppl[all_infected, ])
    
    # Store variant prevalence as simple unnamed vector for use when importing cases
    if (n_infected == 0) variant_prevalence = c(1, rep(0, nrow(p$variants) - 1))
    else variant_prevalence = as.numeric(table(all_variants) / n_infected)
    
    # ---- New local infections ----
    
    # Reduction of contacts due to NPIs
    contact_reduction = min(1 - p$npi_effect[i] * p$npi_scaler, 1)
    
    # All contacts between infectious ('from') and susceptible ('to') people 
    #
    # NOTES: 
    #  1) We remove all contacts with isolated people (ie assume they are actually isolating)
    #  2) The slicing here is a basic implementation of social distancing - remove edges at random
    infectious_contacts = network$edge_list %>%
      slice_sample(prop = contact_reduction) %>%
      filter(from %in% setdiff(all_infectious, currently_isolated), 
             to   %in% all_susceptible)
    
    # Variant index of infectious individuals (can impact both infectiousness and susceptibility)
    variant_idx = match(ppl[infectious_contacts$from, variant], p$variants$id)
    
    # TODO: Consider different beta values based on network layer
    
    # Infectiouness is a multiplier of viral load, seasonality, and viral variant
    beta = ppl[infectious_contacts$from, viral_load] * 
      p$beta * p$seasonality[i] * p$variants$infectivity[variant_idx]
    
    # Calculate immunity for each exposure (depends on variant exposed to)
    immunity_output   = fn_immunity(p, ppl, variant_prevalence, infectious_contacts)
    immunity_contacts = immunity_output$immunity_contacts
    
    # Update immunity status in ppl datatable
    ppl = immunity_output$ppl
    
    # Draw random numbers to determine which contacts lead to transmission
    transmission_df = immunity_contacts %>% 
      mutate(probability = pmin(beta * (1 - immunity), 1),  # Probability of transmission in each contact
             rand = runif(n()), 
             transmission = rand < probability) %>% 
      filter(transmission == TRUE) %>%  # Only interested in transmission events
      rename(transmission_to = to) %>%  # Just for readability / ease of interpretation
      group_by(transmission_to) %>%  # Uniqueness then guaranteed
      summarise(transmission_from = first(from)) %>%  # Only consider first transmission if multiple
      mutate(variant = ppl[transmission_from, variant]) %>%  # Variant is passed down from infector
      as.data.table()
    
    # NOTE: If desirable, we could here store number of people infected by each individual - could 
    #       be useful if looking at the effects on removing key individuals from transmission chain
    
    # Extract the IDs and variants of new local transmissions
    local_id      = transmission_df$transmission_to 
    apply_variant = transmission_df$variant
    
    # Begin infections for these people
    ppl = fn_begin_infection(o, p, ppl, local_id, apply_variant)
    
    # Store number of new local infections among susceptibles
    m = update_output(m, p, i, "new_local_infections", ppl[local_id, ])
    
    # ---- New imported infections ----
    
    # Number of infection imports this time step
    #
    # NOTE: Rounding up here as we always want at least 1 and sample_vec rounds down by default
    n_imports = ceiling(p$import_constant * length(all_susceptible))
    
    # Sample indices - random for now, but may want to consider other factors later
    import_id = sample_vec(all_susceptible, n_imports, prob = 1 - susceptible_immunity)
    
    # Select virus variants for these individuals based on the proportion currently in circulation
    apply_variant = sample_vec(p$variants$id, n_imports, 
                               prob = variant_prevalence, replace = TRUE)
    
    # Begin infections for these people
    ppl = fn_begin_infection(o, p, ppl, import_id, apply_variant)
    
    # Store number of new imported infections
    m = update_output(m, p, i, "new_importations", ppl[import_id, ])
    
    # ---- Introduce new variant ----
    
    # All newly infected people (including imports) we could use to initiate new mutation
    newly_infected_id = c(local_id, import_id)
    new_variant_id    = integer(0)
    
    # Check whether any variants are due to be imported in this time step
    if (i %in% p$variants$import_day) {
      
      # Details of the variant to be imported
      this_variant = p$variants[p$variants$import_day == i, ]
      
      # Do we have enough new infections to assign the new variant to
      n_imports_remain = max(0, this_variant$import_number - length(newly_infected_id))
      
      # Deal with case where we need to introduce elsewhere
      if (n_imports_remain > 0) {
        
        # Can reassign variant for the already infected OR force infect susceptibles - we do the latter
        newly_forced_id = sample_vec(all_susceptible, n_imports_remain, 
                                     prob = 1 - susceptible_immunity)
        
        # Begin infections for these people
        ppl = fn_begin_infection(o, p, ppl, newly_forced_id, this_variant$id)
        
        # Extend the people we can choose from to 'reassign' the new variant to
        newly_infected_id = c(newly_infected_id, newly_forced_id)
      }
      
      # Sample this many people from those recently infected 
      new_variant_id = sample_vec(newly_infected_id, this_variant$import_number)
      
      # Assign the new mutation variant to these individuals
      ppl[new_variant_id, variant := this_variant$id]
    }
    
    # ---- Total number of new infections ----
    
    # All newly infected people (including all imports)
    all_new_infections = c(local_id, import_id, new_variant_id)
    m = update_output(m, p, i, "all_new_infections", ppl[all_new_infections, ])
    
    # ---- Effective reproduction number ----
    
    # INTERPRETATION: New onwards infections per person over the course of one's infection
    
    # Currently infectious in this time step
    infectious_time[i] = length(all_infectious)
    
    # Avoid divide by 0 (occurs if on one is currently infectious)
    if (i > p$infectious_period) {
      
      # Number of people infectious infectious_period days ago
      prevously_infectious = infectious_time[i - p$infectious_period]
      
      # Calculate effective reproductive number using number newly infected this time step
      R_effective = length(newly_infected_id) / (prevously_infectious / p$infectious_period)
      
      # Store in model output datatable
      m = update_output(m, p, i, "R_effective", R_effective)
    }
    
    # ---- Hospital metrics ----
    
    # Number of people currently in hospital and ICU
    #
    # NOTE: Metrics for which a person is only counted once (eg admissions) are handled in update_state
    m = update_output(m, p, i, "hospital_beds", ppl[care_state %in% c("hospital", "rehosp"), ])
    m = update_output(m, p, i, "icu_beds",      ppl[care_state == "icu", ])
    
    # ---- Testing, diagnosis, and isolation ----
    
    # Determine test and diagnosis probability for newly symptomatic people
    ppl = fn_test_diagnose(p, ppl, i)
    
    # Store total number of those diagnosed this time step
    m = update_output(m, p, i, "confirmed", ppl[diagnosis_date == i, ])
    
    # Isolate newly diagnosed individuals where appropriate
    ppl = fn_isolate(p, ppl, i)
    
    # ---- Vaccination ----
    
    # Vaccinate priority groups according to p$vaccine_coverage
    ppl = fn_vaccinate(p, ppl, n_vaccine_group, i)
    
    # Store the number of people vaccinated today and in total
    m = update_output(m, p, i, "n_vaccinated",     ppl[days_vaccinated == 0, ])
    m = update_output(m, p, i, "total_vaccinated", ppl[days_vaccinated >= 0, ])
    
    # People receiving doses: first dose, subsequent doses, and booster doses
    vaccine_dose = ppl[days_vaccinated %in% c(0, p$vaccine$subsequent_dose_days), ]
    booster_dose = ppl[days_booster == 0, ]
    
    # Store number of vaccine doses administered in this time step
    m = update_output(m, p, i, "n_doses", rbind(vaccine_dose, booster_dose))
    
    # ---- Other updates ----
    
    # Move one day closer to next 'event' (see model_flow config file for details)
    ppl[!is.na(days_next_event), days_next_event := days_next_event - 1]
    
    # Sanity check: throw an error if days to next event has become negative
    if (nrow(ppl[!is.na(days_next_event) & days_next_event < 0, ]) > 0)
      stop("Model error: negative days until next event")
    
    # Update disease and care state then sample duration of next phase
    updates = update_state(p, ppl, i, m = m)
    
    # Apply updates to ppl array and metric list, m
    ppl = updates$ppl
    m   = updates$m
    
    # Increment the number of days infected and infectious, or recovered
    ppl[!is.na(days_infected),   days_infected   := days_infected   + 1]
    ppl[!is.na(days_infectious), days_infectious := days_infectious + 1]
    
    # Update viral load in next time step for all infected people
    ppl = fn_viral_load(p, ppl)
    
    # Also increment the number of days since recovery and/or vaccination
    ppl[!is.na(days_recovered),  days_recovered  := days_recovered  + 1]
    ppl[!is.na(days_vaccinated), days_vaccinated := days_vaccinated + 1]
    ppl[!is.na(days_booster),    days_booster    := days_booster    + 1]
    
    # Calculate number in each disease and care state
    for (state in p$model_states$disease) states[[state]][i] = length(ppl[disease_state == state, id])
    for (state in p$model_states$care)    states[[state]][i] = length(ppl[care_state    == state, id])
    
    # TODO: Replace dead people? If yes, do that here
    if (p$replace_dead == TRUE)
      stop("Haven't yet implemented the 'replace_dead' feature")
    
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
  if (do_plot) plot_disease_state(o, p, states)
  
  # Prepare final model output
  results = format_output(o, p, m, network, ppl, yaml, seed)
  
  return(results)
}

# ---------------------------------------------------------
# Set initial conditionals for people array
# ---------------------------------------------------------
initiate_ppl = function(o, p, verbose) {
  
  if (verbose != "none")
    message("  > Initiating population")
  
  # Initiate people datatable	
  ppl = data.table(	
    id                = integer(p$population_size),    # ID of individual, can potentially be replaced by rownumber	
    age               = integer(p$population_size),    # Integer age, will map onto age classes	
    comorbidities     = logical(p$population_size),    # Does this person have relevant comorbidities
    healthcare_worker = logical(p$population_size),    # Is this individual working in healthcare
    disease_state     = character(p$population_size),  # Disease state of infected individual pre-symptomatic, asymptomatic, mild, severe, or critical	
    care_state        = character(p$population_size),  # Is this individual in a hospital, care home or other?	
    prognosis_state   = character(p$population_size),  # Prognosis state - used as alternative to viral load	
    days_infected     = numeric(p$population_size),    # Days since infection
    days_infectious   = numeric(p$population_size),    # Days since becoming infectious (after latent period)	
    n_infections      = integer(p$population_size),    # Total number of infections experienced per person	
    viral_load        = numeric(p$population_size),    # Viral load for infected individuals	
    variant           = character(p$population_size),  # Variant infected with (not modelling mutations, so this requires importation for each variant)	
    mass_test_accept  = logical(p$population_size),    # Whether this individual would be willing to take part in mass testing
    test_date         = numeric(p$population_size),    # Potential data for a test (all infected individuals are assigned this)
    diagnosis_date    = numeric(p$population_size),    # Date of diagnosis (ie positive test result)
    days_next_event   = numeric(p$population_size),    # Days until next event occurs	
    days_recovered    = numeric(p$population_size),    # Number of days since recovery - used to define immunity profile	
    immune_state      = numeric(p$population_size),    # Captures immune state of individual and other relevant comorbidities like blood group, Il6 to Il10 ratio, and neanderthal recombinant genome	
    vaccine_group     = character(p$population_size),  # Priority group for vaccination (group 1 highest priority)	
    vaccine_accept    = logical(p$population_size),    # Whether this individual with accept vaccination
    booster_accept    = logical(p$population_size),    # Whether this individual with accept vaccination booster dose
    booster_due       = numeric(p$population_size),    # Date next booster dose is due to be administered
    days_vaccinated   = numeric(p$population_size),    # Number of days since vaccination - used to determine increasing then waning effect
    days_booster      = numeric(p$population_size),    # Number of days since vaccination booster dose
    vaccine_effect    = numeric(p$population_size))    # Current level of vaccine effect based on which vaccine and when vaccinated
  
  # ---- Set initial values ----
  
  # Sample ages from some distribution
  init_age = sample_vec(x    = p$age$all, 
                        size = p$population_size, 
                        prob = p$demography, 
                        replace = TRUE)
  
  # Initate basic values in people datatable
  ppl[, id  := 1 : p$population_size]
  ppl[, age := init_age]
  
  # Start by assuming everyone is susceptible
  ppl[, comorbidities     := FALSE]   # Assigned below
  ppl[, healthcare_worker := FALSE]   # Assigned below
  ppl[, disease_state     := "susc"]  # Assume everyone starts susceptible
  ppl[, care_state        := "none"]
  ppl[, prognosis_state   := "none"]
  ppl[, days_infected     := NA]
  ppl[, days_infectious   := NA]
  ppl[, n_infections      := 0L]
  ppl[, viral_load        := NA]
  ppl[, variant           := "none"]
  ppl[, mass_test_accept  := FALSE]   # Assigned below
  ppl[, test_date         := NA]
  ppl[, diagnosis_date    := NA]
  ppl[, days_next_event   := NA]
  ppl[, days_recovered    := NA]
  ppl[, immune_state      := 0]       # Assume everyone starts with no immunity
  ppl[, vaccine_group     := "none"]  # Assigned below
  ppl[, vaccine_accept    := FALSE]   # Assigned below
  ppl[, booster_accept    := FALSE]   # Assigned below
  ppl[, booster_due       := NA]
  ppl[, days_vaccinated   := NA]
  ppl[, days_booster      := NA]
  ppl[, vaccine_effect    := 0]
  
  # ---- Assign risk groups ----
  
  # Sanity check that all defined risk groups are able to modelled
  if (!all(p$risk_groups$name %in% names(ppl)))
    stop("Unrecognised risk groups")
  
  # Loop through the user-defined risk groups
  for (i in seq_len(nrow(p$risk_groups))) {
    this_group = p$risk_groups[i, ]
    
    # Age-corrected probability of being in this risk group (with compliment)
    risk_prop = c(this_group$probability, 1 - this_group$probability)
    
    # Sample TRUE and FALSE for all in this age group accordingly 
    ppl[age >= this_group$age_lower & age <= this_group$age_upper, 
        (this_group$name) := sample_vec(c(TRUE, FALSE), 
                                        size = .N,
                                        prob = risk_prop, 
                                        replace = TRUE)]
  }
  
  # ---- Set vaccination priority policy and acceptance ----
  
  # Iterate through the different priority groups
  # 
  # NOTE: We do last at first so any people in multiple groups are placed in the highest possible group
  for (i in 1 : p$vaccine$n_groups) {
    this_group = p$vaccine$groups[p$vaccine$n_groups - i + 1, ]
    
    # Evaluate the condition and set priority for any that satisfy
    ppl[eval(parse(text = this_group$condition)), vaccine_group := this_group$id]
  }
  
  # Again iterate through priority groups - this time to set vaccine acceptance probabilities
  #
  # NOTE: We do this after all groups are formalised to prevent double counting
  for (group_id in p$vaccine$groups$id) {
    
    # Proportion of this group that would accept a vaccine (ie max coverage)
    acceptance = p$vaccine$groups[id == group_id, max_coverage]
    booster    = p$vaccine$booster_groups[vaccine_group == group_id, probability]
    
    # Sample logical and populate vaccine_accept variable
    ppl[vaccine_group == group_id, 
        vaccine_accept := sample_vec(c(TRUE, FALSE), size = .N, replace = TRUE, 
                                     prob = c(acceptance, 1 - acceptance))]
    
    # Do similar for booster_accept variable, but only for those that accept vaccine in general
    ppl[vaccine_group == group_id & vaccine_accept == TRUE, 
        booster_accept := sample_vec(c(TRUE, FALSE), size = .N, replace = TRUE, 
                                     prob = c(booster, 1 - booster))]
  }
  
  # ---- Mass testing acceptance ----
  
  # Number of modelled people willing to be regularly mass tested
  n_mass_test_accept = round(p$population_size * p$testing$mass_testing$probability)
  
  # Probability of acceptance per person (dependent on age)
  p_mass_test_accept = p$testing$mass_testing$age_prob[ppl[, age] + 1]
  
  # Sample mass_test_acceptance of the population with age-related probabilities
  accept_id = wrswoR::sample_int_crank(n = p$population_size,
                                       size = n_mass_test_accept,
                                       prob = p_mass_test_accept)
  
  # Assign these people to accept mass testing
  ppl[accept_id, mass_test_accept := TRUE]
  
  return(ppl)
}

# ---------------------------------------------------------
# Initially infect a subset of people to kick start epidemic
# ---------------------------------------------------------
initiate_epidemic = function(o, p, ppl) {
  
  # ---- Initiate previously vaccinated ---
  
  # Iterate through the vaccine priority groups
  for (group_name in p$vaccine$groups$id) {
    group = p$vaccine$groups[id == group_name, ]
    
    # Skip this if we have zero max coverage for this group
    if (group$max_coverage > 1e-6) {
      
      # IDs of all those in this priority group who will accept the vaccine
      group_id = ppl[vaccine_group == group_name & vaccine_accept == TRUE, id]
      
      # Proportion of group and absolute number previously vaccinated
      p_vaccinate = group$init_coverage / group$max_coverage
      n_vaccinate = round(length(group_id) * p_vaccinate)
      
      # Sample IDs of those previously vaccinated - uniformly for all in this priority group
      vaccinate_id = sample_vec(group_id, size = n_vaccinate)
      
      # Sample number of days each person has been vaccinated for
      ppl[vaccinate_id, days_vaccinated := -sample_vec(group$init_start : min(group$init_end, -1), 
                                                       size = .N, replace = TRUE)]
    }
  }
  
  # ---- Set booster dose info for those already vaccinated ----
  
  # Booster delivery details for all those who will recieve booster dose(s)
  booster_df = ppl[days_vaccinated > 0 & booster_accept == TRUE] %>% 
    select(id, vaccine_group, days_vaccinated, days_booster, booster_due) %>% 
    left_join(p$vaccine$booster_groups[, -"probability"], by = "vaccine_group")
  
  # Assign first booster delivery date, bounded below if force-starting
  booster_df[, booster_due := pmin(pmax(cycle_period - days_vaccinated, start), force_start, na.rm = TRUE)]
  
  # If first booster is due in the past, determine most recent dose and assign date of next dose 
  booster_df[booster_due < 0, days_booster := -booster_due - (-booster_due %/% cycle_period) * cycle_period]
  booster_df[booster_due < 0, booster_due := cycle_period - days_booster]
  
  # Apply these details to the main ppl datatable (avoid starting on day zero)
  ppl[id %in% booster_df$id, days_booster := pmax(booster_df$days_booster, 1, na.rm = TRUE)]
  ppl[id %in% booster_df$id, booster_due  := pmax(booster_df$booster_due, 1)]
  
  # ---- Initiate previously infected ----
  
  # Number of people infected, and IDs of those selected (uniform selection)
  n_prev_infected  = round(p$population_size * p$previously_infected)
  prev_infected_id = sort(sample.int(p$population_size, size = n_prev_infected))
  
  # Begin infections for these people
  ppl = fn_begin_infection(o, p, ppl, prev_infected_id, p$variants$id[1])
  
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
    ppl = update_state(p, ppl, 1)$ppl
    
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
  ppl = fn_immunity(p, ppl, variant_prevalence)$ppl
  
  # Sanity check for numerical vaccine_effect
  if (any(is.na(ppl$vaccine_effect)))
    stop("Non-numerical vaccine effect")
  
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
update_output = function(m, p, date_idx, update_metric, to_count) {
  
  # Loop through each metric defined within update_metric (comma seperated)
  update_metrics = unlist(str_split(str_remove_all(update_metric, " "), ","))
  for (update_metric in update_metrics) {
    
    # Loop through all groupings defined for this metric
    metric_groupings = p$metrics$groupings[names(p$metrics$groupings) == update_metric]
    for (grouping in metric_groupings) {
      
      # For metrics with no possible disaggregation, append single value
      if (grouping == "na") {
        m[metric == update_metric & date == date_idx, value := to_count]
        
        # For metrics we have chosen not to disaggregate, count all values
      } else if (grouping == "none") {
        m[metric == update_metric & date == date_idx, value := nrow(to_count)]
        
      } else {  # Otherwise we want to count each grouping...
        
        # Count number of occurrences of each group
        group_count = table(to_count[, ..grouping])
        
        # Insert these values into output datatable
        m[metric == update_metric &         # This metric
            date == date_idx &              # This time step
            group %in% names(group_count),  # This subset of 'groups'
          value := group_count]
      }
    }
  }
  
  # Check the update: m[metric == update_metric & date == date_idx, ]
  
  return(m)
}

# ---------------------------------------------------------
# Prepare final model output
# ---------------------------------------------------------
format_output = function(o, p, m, network, ppl, yaml, seed) {
  
  # Insert seasonality profile - this is actually input, not output, but is helpful for visualisation
  m[metric == "seasonality", value := p$seasonality[1 : p$n_days]]
  
  # Similar for NPI effect on contact reduction
  m[metric == "contact_reduction", value := p$npi_effect * p$npi_scaler]
  
  # Convert number of cases per variant to variant prevalence
  m[metric == "variant_prevalence", value := 100 * value / max(sum(value), 1), by = "date"]
  
  # Trivialise first R_eff values - can only be calculated after infectious_period days
  m[metric == "R_effective" & date <= p$infectious_period, value := NA]
  
  # Append scenario and seed columns to output datatable
  m$scenario = p$.id
  m$seed = seed
  
  # Append ages of contacts to edge list
  network = network$edge_list %>% 
    mutate(from_age = ppl[from, age],
           to_age   = ppl[to, age], .before = layer)
  
  # Combine inputs and outputs into one list
  results = list(yaml    = yaml$raw, 
                 input   = yaml$parsed, 
                 network = network,
                 output  = m)
  
  return(results)
}

# ---------------------------------------------------------
# Update disease and care states for all relevant individuals
# ---------------------------------------------------------
update_state = function(p, ppl, date_idx, m = NULL) {
  
  # IDs of people that need to be updated
  update_id = ppl[days_next_event == 0, id]
  
  # Skip this process if no one to update
  if (length(update_id) > 0) {
    
    # Join with model flows dataframe to determine next state
    update_df = left_join(ppl[update_id], p$model_flows, 
                          by = qc(disease_state, care_state, prognosis_state))
    
    # If m datatable is provided we'll want to count metrics
    if (!is.null(m)) {
      
      # Count and store values of each metric in model output
      for (update_metric in unique(na.omit(update_df$metric)))
        m = update_output(m, p, date_idx, update_metric, update_df[metric == update_metric, ])
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
    
    # ---- Initaite infectious period and possible test date ----
    
    # Newly infectious people are those going INTO the pre-symptomatic phase
    new_infectious_id = update_df[next_disease == "presym", id]
    
    # Signal start of infectious period by starting 'days infectious' count
    ppl[new_infectious_id, days_infectious := 0]  # Incremented at end of daily loop
    
    # Newly symptomatic people are those going FROM the pre-symptomatic phase
    new_symptoms_id = update_df[disease_state == "presym", id]
    
    # Date of potential test and diagnosis - see fn_test_diagnose for probability
    ppl[new_symptoms_id, test_date := date_idx + p$diagnosis_delay]
  }
  
  return(list(ppl = ppl, m = m))
}

# ---------------------------------------------------------
# Initiate new infections for newly infected individuals
# ---------------------------------------------------------
fn_begin_infection = function(o, p, ppl, id, apply_variant) {
  
  # We signal this new infection by starting 'days infected' count
  ppl[id, days_infected := 0]  # Incremented at end of daily loop
  
  # Increment number of infections experienced thus far
  ppl[id, n_infections := n_infections + 1L]
  
  # Reset counter for days recovered if necessary
  ppl[id, days_recovered := NA]
  
  # Sanity check: we need a variant defined for each person
  if (!(length(apply_variant) %in% c(length(id), 1)))
    stop("Inconsistent number of variants defined")
  
  # Apply the precalculated variant for this individial
  ppl[id, variant := apply_variant]
  
  # All newly infected people go to the latent phase
  ppl[id, disease_state := "latent"]
  
  # Define how log this latency phase will last
  ppl[id, days_next_event := pmax(1, round(p$duration$latency(id)))]
  
  # Sample a prognosis for each infected individual
  ppl[id, prognosis_state := fn_prognosis(o, p, ppl[id, ])]
  
  return(ppl)
}

# ---------------------------------------------------------
# Generate prognosis state for newly infected individuals
# ---------------------------------------------------------
fn_prognosis = function(o, p, new_infected) {
  
  # ---- Construct cumulative probability matrix ----
  
  # Proportion of all cases that are asymptomatic
  p_asym = rep(p$proportion_asymptomatic, length(p$age$groups))
  
  # Proportion of symptomatic cases that are severe and require hospitalisation (by age group)
  # 
  # NOTE: This considers (potentially) less contacts in older ages and a calibration scaling factor
  symptom_severe = pmin(p$severe_symptom_age * p$age_correction, 1)
  
  # Proportion of mild cases: symptomatic, non-severe
  p_mild   = (1 - p_asym) * (1 - symptom_severe)
  p_severe = (1 - p_asym) * symptom_severe
  
  # Proportion of cases that are severe that will / will not become critical
  p_critical     = p_severe * pmin(p$critical_severe_age, 1)
  p_non_critical = p_severe * (1 - pmin(p$critical_severe_age, 1))
  
  # Split these by hospital-seeking behaviour 
  p_seek_S   = p_non_critical * p$seek_hospital
  p_noseek_S = p_non_critical * (1 - p$seek_hospital)
  
  # Of the critical cases admitted to ICU, who will recover and who will die
  p_seek_C = p_critical * p$seek_hospital * (1 - p$death_critical_icu * p$death_critical_age)
  p_seek_D = p_critical * p$seek_hospital * p$death_critical_icu * p$death_critical_age
  
  # Of the critical cases outside of the care system, who will recover and who will die
  p_noseek_C = p_critical * (1 - p$seek_hospital) * (1 - p$death_critical_non_icu * p$death_critical_age)
  p_noseek_D = p_critical * (1 - p$seek_hospital) * p$death_critical_non_icu * p$death_critical_age
  
  # Re-aggregate these probabilities to check we haven't made silly mistakes
  p_seek   = p_seek_S   + p_seek_C   + p_seek_D
  p_noseek = p_noseek_S + p_noseek_C + p_noseek_D
  
  # Sanity check: Severe people should either seek or not seek hospital care
  if (any(abs(p_seek + p_noseek - p_severe) > 1e-6))
    stop("Hospital seeking probabilities must align for all age groups")
  
  # Sanity check: Check that these sum to one and nothing weird is going on
  if (any(abs(p_asym + p_mild + p_seek + p_noseek - 1) > 1e-6))
    stop("Prognosis probabilities must sum to 1 for all age groups")
  
  # Combine probabilities into dataframe (age group x prognosis state)
  prognosis_df = data.table(asym        = p_asym, 
                            mild        = p_mild, 
                            severe_care = p_seek_S, 
                            crit_care   = p_seek_C, 
                            dead_care   = p_seek_D, 
                            severe_home = p_noseek_S, 
                            crit_home   = p_noseek_C,
                            dead_home   = p_noseek_D)
  
  # Convert to a matrix
  prognosis_names  = names(prognosis_df)
  prognosis_matrix = as.matrix(prognosis_df)
  
  # Number of prognoses possible - used to index matrices
  n_prognoses = length(prognosis_names)
  
  # ---- Sample to determine age-related prognosis ----
  
  # Preallocate vector of prognoses
  prognosis = rep(NA, nrow(new_infected))
  
  # NOTE: new_infected is just a slice of ppl and include all variables
  
  # Age and comorbidity status of newly infected individuals
  age     = new_infected[, age]
  variant = new_infected[, variant]
  
  # Represent increase risk of comorbidities by shifting up one age class
  age[new_infected$comorbidities] = age[new_infected$comorbidities] + 10
  age = pmin(age, max(p$age$all))
  
  # Convert age to age group index to reference
  age_bin = floor(age / 10) + 1 
  
  # Probability of severe disease
  prob_severe     = rowSums(prognosis_matrix[, p$prognosis_idx$severe])
  prob_non_severe = rowSums(prognosis_matrix[, p$prognosis_idx$non_severe])
  
  # Repeat calculation for all currently circulating variants
  for (this_variant in unique(variant)) {
    
    # Reset prognosis matrix
    prognosis_variant = prognosis_matrix
    
    # Increase in severity for this variant
    variant_severity = p$variants[id == this_variant, severity]
    
    # This needs to be bounded above otherwise we can get negative non-severe probabilities
    variant_severity = pmin(prob_severe * variant_severity, 1) / prob_severe
    
    # Convert into matrix for element-by-element multiplication for subset of prognosis matrix
    severe_factor = matrix(data = variant_severity, 
                           nrow = length(p$age$groups), 
                           ncol = sum(p$prognosis_idx$severe))
    
    # Each severe case is variant_severity times more likely
    prognosis_variant[, p$prognosis_idx$severe] = 
      prognosis_variant[, p$prognosis_idx$severe] * severe_factor
    
    # Scale factor needed for non-severe prognoses to acheive probability of 1 for each age
    prob_severe_variant = rowSums(prognosis_variant[, p$prognosis_idx$severe])
    non_severe_factor   = (1 - prob_severe_variant) / prob_non_severe
    
    # Multiply through with non-severe factor to obtain prognosis matrix for this variant
    prognosis_variant[, p$prognosis_idx$non_severe] = 
      prognosis_variant[, p$prognosis_idx$non_severe] * non_severe_factor
    
    # We should not have any negative values here
    if (any(prognosis_variant < -1e-6))
      stop("Negative values when calculating disease severity for variant: ", this_variant)
    
    # Cumulative sum the probabilities to generate a matrix that can be sampled
    prognosis_variant = pmax(prognosis_variant, 0) # Bound below in case of rounding errors
    prognosis_cumsum  = rowCumsums(prognosis_variant)
    
    # Sanity check that probabilities still sum to 1
    if (any(abs(prognosis_cumsum[, n_prognoses] - 1) > 1e-6))
      stop("Prognosis probabilities must sum to 1 for all age groups")
    
    # Subset of newly infecteds with this variant
    new_infected_variant = new_infected[variant == this_variant]
    
    # Index of these individuals within the new_infected datatable
    new_infected_idx = match(new_infected_variant[, id], new_infected[, id])
    
    # Generate random number vector
    random_number = runif(nrow(new_infected_variant))
    
    # Select appropriate probability index from each row based on random number
    prognosis_bool = (random_number < prognosis_cumsum[age_bin[new_infected_idx], ]) * 1
    
    # A bit of matrix manipulation to extract prognosis index - maybe theres a nicer way to do this
    prognosis_idx = rowCumsums(matrix(prognosis_bool, ncol = n_prognoses))
    prognosis_idx = n_prognoses - prognosis_idx[, n_prognoses] + 1
    
    # Convert index to prognosis state name for readability
    prognosis[new_infected_idx] = prognosis_names[prognosis_idx]
  }
  
  # ---- Potentially avoid disease if vaccinated ----
  
  # Immunity to severe disease can be reduced if variant infected with is immuno-escaping
  immunity_disease = new_infected[, .(id, variant, vaccine_effect)] %>%
    mutate(escape = 1 - p$variants$immuno_escape[match(variant, p$variants$id)], 
           value  = vaccine_effect * escape) %>%
    pull(value)
  
  # Overall probability of being asymptomatic (considering standard proportion_asymptomatic)
  #
  # NOTE: This is scaled by the current 'effect' of the vaccine
  probability_asym = immunity_disease / p$vaccine$efficacy * 
    (1 - (1 - p$vaccine$disease_prevention) * (1 - p$proportion_asymptomatic))
  
  # Skip if none of these people have been vaccinated
  if (sum(probability_asym) > 0) {
    
    # Draw random number to determine if this person should be asymptomatic due to vaccine effect
    asym_vaccine_effect = runif(nrow(new_infected)) < probability_asym
    
    # Reassign prognosis to asymptomatic disease
    prognosis[asym_vaccine_effect] = "asym"
  }
  
  return(prognosis)
}

# ---------------------------------------------------------
# Which infected people will get tested AND diagnosed, and when
# ---------------------------------------------------------
fn_test_diagnose = function(p, ppl, date_idx) {
  
  # ---- Regular testing ----
  
  # All people who could potentially get tested today
  #
  # NOTE: Occasionally people recover before getting tested, hence days_infectious check
  test_potential = ppl[test_date == date_idx & !is.na(days_infectious), ]
  
  # Skip this process if no potential testers identified
  if (nrow(test_potential) > 0) {
    
    # Preallocate vector for probability of test AND diagnosis
    dx_probability = rep(1, nrow(test_potential))
    
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
    random_number = runif(nrow(test_potential))
    
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
  
  return(ppl)
}

# ---------------------------------------------------------
# Isolate newly diagnosed individuals where appropriate
# ---------------------------------------------------------
fn_isolate = function(p, ppl, date_idx) {
  
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
  
  return(ppl)
}

# ---------------------------------------------------------
# Vaccinate eligible people (according to number_vaccines or vaccine_coverage)
# ---------------------------------------------------------
fn_vaccinate = function(p, ppl, n_vaccine_group, date_idx) {
  
  # Coverage difference between yesterday and today
  #
  # NOTE: Binding initial coverage so we consistently index a 2D matrix
  coverage_matrix = rbind(p$vaccine$groups$init_coverage, p$vaccine$coverage)
  coverage_diff   = coverage_matrix[date_idx : (date_idx + 1), ] %>%
    matrix(nrow = 2, dimnames = list(1 : 2, colnames(coverage_matrix)))
  
  # Iterate through the vaccine priority groups
  for (group in p$vaccine$groups$id) {
    
    # Only continue if we have people to vaccinate at this time step
    if (diff(coverage_diff[, group]) > 0) {
      
      # All people eligible for vaccination in this group at this time step 
      eligible = ppl[vaccine_group == group & vaccine_accept == TRUE & is.na(days_vaccinated), ]
      
      # Number of people we want vaccinated in this group in this time step
      n_group   = n_vaccine_group[vaccine_group == group, N]
      n_desired = round(p$vaccine$coverage[date_idx, group] * n_group)
      
      # Number of people currently vaccinated in this group
      n_current = nrow(ppl[vaccine_group == group & !is.na(days_vaccinated), ])
      
      # How many we will vacinate this time step - bounded above by number available
      n_vaccinate = min(max(n_desired - n_current, 0), nrow(eligible))
      
      # Vaccinate the first chunk of people - randomness is inherent as table is not ordered
      vaccinate_id = eligible[seq_len(n_vaccinate), id]
      
      # Start these people off on their vaccination journey
      ppl[vaccinate_id, days_vaccinated := 0]
    }
  }
  
  # ---- Set next booster dose ----
  
  # NOTE: Booster cycle initiates on day of vaccination UNLESS using force_start, in which case
  #       the first booster is given on day force_start and the cycle continues from that point. 
  
  # Set date of first booster dose (for those just vaccinated initial vaccination)
  booster_1_df = ppl[days_vaccinated == 0 & booster_accept == TRUE] %>% 
    select(id, vaccine_group, days_vaccinated, days_booster, booster_due) %>% 
    left_join(p$vaccine$booster_groups[, -"probability"], by = "vaccine_group") %>%
    mutate(booster_due = pmin(pmax(cycle_period + date_idx, start), force_start, na.rm = TRUE))
  
  # Set date of subsequent booster dose (for those who have just recieved a booster)
  booster_n_df = ppl[booster_due == date_idx] %>% 
    select(id, vaccine_group, days_vaccinated, days_booster, booster_due) %>% 
    left_join(p$vaccine$booster_groups[, -"probability"], by = "vaccine_group") %>%
    mutate(booster_due = cycle_period + date_idx)
  
  # ---- Apply booster dose ----
  
  # If due a booster, initiate days since most recent booster received
  ppl[booster_due == date_idx, days_booster := 0]
  
  # Apply next booster due date calculated above
  ppl[id %in% booster_1_df$id, booster_due := booster_1_df$booster_due]
  ppl[id %in% booster_n_df$id, booster_due := booster_n_df$booster_due]
  
  # Sanity check that no-one has been missed
  if (any(ppl$booster_due <= date_idx, na.rm = TRUE))
    stop("Booster doses have been missed - investigation needed")
  
  return(ppl)
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
# Update immunity state for all relevant individuals
# ---------------------------------------------------------
fn_immunity = function(p, ppl, variant_prevalence, infectious_contacts = NULL) {
  
  # ---- Immunity per exposure ----
  
  # First job is to calculate current 'effect' of vaccination - depends on day of vaccination and any booster
  ppl[days_vaccinated > 0, vaccine_effect := pmax(p$vaccine$booster_profile[days_booster], 
                                                  p$vaccine$profile[days_vaccinated], na.rm = TRUE)]
  
  # Skip this step if no exposure details provided
  if (is.null(infectious_contacts)) immunity_contacts = NULL else {
    
    # For each exposure, append variant details and other factors influencing immunity
    immunity_df = infectious_contacts %>%
      left_join(ppl[, .(id, from_variant = variant)], by = c("from" = "id")) %>%
      left_join(ppl[, .(id, to_variant = variant, vaccine_effect, days_recovered)], by = c("to" = "id")) %>%
      mutate(recovery_effect = p$acquired_immunity[pmax(days_recovered, 1)]) %>%
      select(-days_recovered)
    
    # Level of immuno-escaping depends on variant exposed to
    immunity_df[, escape := 1 - p$variants$immuno_escape[match(from_variant, p$variants$id)]]
    
    # Resest this if susceptible person has been previously infected with the same variant
    immunity_df[from_variant == to_variant, escape := 1]
    
    # Now we can compute level of immunity for each exposure
    immunity_contacts = immunity_df %>%
      mutate(immunity = pmax(escape * vaccine_effect * p$vaccine$transmission_blocking, 
                             escape * recovery_effect, na.rm = TRUE)) %>%
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
  ppl[, immune_state := pmax(vaccine_effect * p$vaccine$transmission_blocking * (1 - sum(pop_escape)),
                             immunity_acquired, na.rm = TRUE)]
  
  # TODO: If number of infections is above a certain cap, set immune_state to 1 with no waning?
  
  return(list(ppl = ppl, immunity_contacts = immunity_contacts))
}

