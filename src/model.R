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

# TODO: 
# 1) Change in disease severity - either through vaccination or due to infection with variant
# 2) Variable 'days_infected' now seems surpless to requirements - think this can be removed

# ---------------------------------------------------------
# Main function: The guts of the model
# ---------------------------------------------------------
model = function(o, p, seed = NA, short_run = FALSE, do_plot = FALSE, verbose = "date") {
  
  if (verbose != "none")
    message(" - Running model")
  
  # Set random seed generator if input defined
  if (!is.na(seed)) set.seed(seed)
  
  # ---- Model set up ----
  
  # Number of days to run the model for (less when calibrating)
  if (short_run) run_time = length(o$dates_data) else run_time = o$n_dates
  
  # Initiate people array
  ppl = initiate_ppl(o, p)
  
  # Generate full contact network outside of main model loop
  network = create_network(o, p, ppl, do_plot, verbose)
  
  # Calculate number of people in each vaccine prioirty group outside of main loop
  n_vaccine_group = ppl[, .N, by = vaccine_priority]
  
  # Group certain states for easy access and readability
  states = list(infected   = qc(latent, presym, asym, mild, severe, crit), 
                infectious = qc(presym, asym, mild, severe, crit), 
                in_care    = qc(hospital, icu, rehosp))
  
  # Also append vectors for all disease and care states
  for (state in o$all_states)
    states[[state]] = rep(0, run_time)
  
  # Preallocate model output datatable (denoted 'm' for short)
  m = preallocate_output(o, run_time)
  
  # Scale interventions by population scaler
  #
  # NOTE: Done here and not in get_parameters as these values may be aletered by the user
  p$number_diagnosed = round(p$number_diagnosed / p$population_scaler)
  p$number_vaccines  = round(p$number_vaccines  / p$population_scaler)
  
  # ---- Main model loop ----
  
  if (verbose != "none")
    message("  > Simulating daily time steps")
  
  # Initiate progress bar if desired
  if (verbose == "bar")
    pb = start_progress_bar(run_time - p$import_date_initial + 1)
  
  # Initially infect a subset of people to kick start epidemic
  #
  # NOTE: Must be called after creating network as we need age-correction factor for prognoses
  p$age_correction = network$age_correction
  p$age_contacts   = network$age_contacts
  ppl = initiate_epidemic(o, p, ppl)  
  
  # Iterate through time steps - starting from first imported case(s)
  for (i in p$import_date_initial : run_time) {
    
    # Report progress by printing date (useful when running on the cluster)
    if (verbose == "date") 
      message("   - Simulating day ", i, ": ", o$dates_all[i])
    
    # Copy people datatable for the next time step
    ppl_next = data.table::copy(ppl)
    
    # ---- Basic epidemiological indicators ----
    
    # Total number of people currently infected and infectious
    all_infected   = ppl[disease_state %in% states$infected,   id]
    all_infectious = ppl[disease_state %in% states$infectious, id]  # Including those in isolation
    
    # Easy access total number of people currently infected
    n_infected = length(all_infected)
    
    # Store total number of people infected and infectious
    m = update_output(o, p, m, i, "currently_infected",   ppl[all_infected, ])
    m = update_output(o, p, m, i, "currently_infectious", ppl[all_infectious, ])
    
    # Total number of people currently in isolation
    currently_isolated = ppl[care_state %in% "iso", id]
    m = update_output(o, p, m, i, "currently_isolated", ppl[currently_isolated, ])
    
    # Total number of people currently infected with symptoms (mild, severe, or critical)
    currently_symptomatic = ppl[disease_state %in% c("mild", "severe", "crit"), id]
    m = update_output(o, p, m, i, "currently_symptomatic", ppl[currently_symptomatic, ])
    
    # All susceptible people and their level of immunity (from previous infection and/or vaccination)
    all_susceptible      = ppl[disease_state == "susc", id]
    susceptible_immunity = ppl[all_susceptible, immune_state]
    
    # For population susceptibility, weight by immunity state and sum
    population_susceptibility = sum(1 - susceptible_immunity)
    m = update_output(o, p, m, i, "pop_susceptibility", population_susceptibility)
    
    # Prevalence across whole population (proportion currently infected)
    population_prevalence = 100 * n_infected / p$population_size
    m = update_output(o, p, m, i, "pop_prevalence", population_prevalence)
    
    # Sero-prevalence across whole population (proportion that have been infected thus far)
    population_seroprevalence = 100 * nrow(ppl[n_infections > 0, ]) / p$population_size
    m = update_output(o, p, m, i, "seroprevalence", population_seroprevalence)
    
    # Virus variant of each infected individual
    all_variants = factor(ppl[all_infected, variant], levels = o$variant_names)
    m = update_output(o, p, m, i, "variant_prevalence", ppl[all_infected, ])
    
    # Store variant prevalence as simple unnamed vector for use when importing cases
    if (n_infected == 0) variant_prevalence = c(1, rep(0, length(o$variants$mutation)))
    else variant_prevalence = as.numeric(table(all_variants) / n_infected)
    
    # ---- New local infections ----
    
    # TODO: Reduce contacts for mild and severe cases - still needed with decaying viral load over time?
    # TODO: This slicing should be symmetric for social distancing programs, and non-symmetric for facemasks
    # TODO: Improve speed by using specialised datatable functions instead of dplyr functions here
    
    # Level of contact reduction: OSI x calibrated sclaer x NPI adherence
    contact_reduction = p$npi_level[i] * p$npi_adherence[i]
    
    # All contacts infectious ('from') and susceptible ('to') people
    #
    # NOTES: 
    #  1) Remove all contacts with isolated people (ie assume they are actually isolating)
    #  2) The slicing here is a very basic implementation of social distancing - remove edges at random
    infectious_contacts = network$all_contacts %>%
      slice_sample(prop = 1 - contact_reduction) %>%
      filter(from %in% setdiff(all_infectious, currently_isolated), 
             to   %in% all_susceptible)
    
    # Infectiouness is a multiplier of viral load and seasonality
    #
    # NOTES: 
    #  1) Facemasks now considered in the exponent (ie within p$npi_level)
    #  2) Interpret seasonality factor as decrease in transmission probability per contact when 
    #     people mix outside (ie in warmer months) rather than inside (ie in colder months)
    beta = ppl[infectious_contacts$from, viral_load] * p$beta * p$seasonality[i]
    
    # Multiply through by variant infectivity factor (trivial if all have primary variant)
    beta = beta * unname(p$variant_infectivity[ppl[infectious_contacts$from, variant]])
    
    # Susceptibility scaler for each susceptible person
    susceptibility = 1 - ppl[infectious_contacts$to, immune_state]
    
    # Draw random numbers to determine which contacts lead to transmission
    transmission_df = infectious_contacts %>% 
      mutate(probability = pmin(1, susceptibility * beta),  # Probability of transmission in each contact
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
    ppl_next = fn_begin_infection(o, p, ppl_next, local_id, apply_variant, i)
    
    # Store number of new local infections among susceptibles
    m = update_output(o, p, m, i, "new_local_infections", ppl[local_id, ])
    
    # ---- New imported infections ----
    
    # Number of infection imports this time step
    n_imports = p$import_constant * length(all_susceptible)
    
    # Sample indices - random for now, but may want to consider other factors later
    #
    # NOTE: The base sample function rounds n_imports down by default - may want something nicer here
    import_id = sample_vec(all_susceptible, n_imports, prob = 1 - susceptible_immunity)
    
    # Select virus variants for these individuals based on the proportion currently in circulation
    apply_variant = sample_vec(names(o$variants_dict), n_imports, 
                               prob = variant_prevalence, replace = TRUE)
    
    # Begin infections for these people
    ppl_next = fn_begin_infection(o, p, ppl_next, import_id, apply_variant, i)
    
    # Store number of new imported infections
    m = update_output(o, p, m, i, "new_importations", ppl[import_id, ])
    
    # ---- Introduce new variant ----
    
    # All newly infected people (including imports) we could use to initiate new mutation
    newly_infected_id = c(local_id, import_id)
    
    # Check whether variants are due to be imported in this time step
    if (p$import_date_variant[i] == TRUE) {
      
      # The mutation due to be imported in this time step
      mutation_idx  = which(o$variants$import_date == o$dates_all[i])
      mutation_name = o$variants$mutation[mutation_idx]
      
      # Number of people we need to import (predefined in variant_properties file)
      n_imports_variant = ceiling(p$import_variant[mutation_idx])
      
      # Do we have enough new infections to assign the new variant to
      n_imports_remain = max(0, n_imports_variant - length(newly_infected_id))
      
      # Deal with case where we need to introduce elsewhere
      if (n_imports_remain > 0) {
        
        # Can reassign variant for the already infected OR force infect susceptibles - we do the latter
        newly_forced_id = sample_vec(all_susceptible, n_imports_remain, 
                                     prob = 1 - susceptible_immunity)
        
        # Begin infections for these people
        ppl_next = fn_begin_infection(o, p, ppl_next, newly_forced_id, mutation_name, i)
        
        # Update number of new imported infections
        # m$new_importations[i] = length(import_id) + length(newly_forced_id)
        
        # Extend the people we can choose from to 'reassign' the new variant to
        newly_infected_id = c(newly_infected_id, newly_forced_id)
      }
      
      # Sample this many people from those recently infected 
      new_variant_id = sample_vec(newly_infected_id, n_imports_variant)
      
      # Assign the new mutation variant to these individuals
      ppl_next[new_variant_id, variant := mutation_name]
    }
    
    # ---- Effective reproduction number ----
    
    # INTERPRETATION: new onwards infections per person over the course of one's infection
    
    # Avoid divide by 0 (occurs if on one is currently infectious)
    if (length(all_infectious) > 0) {
      
      # Calculate effective reproductive number using number newly infected this time step
      R_effective = length(newly_infected_id) / (length(all_infectious) / p$infectious_period)
      
      # Store in model output datatable
      m = update_output(o, p, m, i, "R_effective", R_effective)
    }
    
    # ---- Hospital metrics ----
    
    # Number of people currently in hospital and ICU
    #
    # NOTE: Metrics for which a person is only counted once (eg admissions) are handled in update_state
    m = update_output(o, p, m, i, "hospital_beds", ppl[care_state %in% c("hospital", "rehosp"), ])
    m = update_output(o, p, m, i, "icu_beds",      ppl[care_state == "icu", ])
    
    # ---- Testing, diagnosis, and isolation ----
    
    # NOTE: Four routes for being tested & diagnosed: 
    #  1) all severe cases admitted to hospital
    #  2) all severe cases that do not seek hospital care
    #  3) some non-severe cases (depending on number of tests available / number diagnoses)
    #  4) some asymptomatic & mild cases through contact tracing
    
    # Define which people to test and diagnose
    to_diagnose = fn_test_diagnose(p, ppl, n_infected, i)
    
    # Calculate severe to non-severe diagnosis ratio - used for projecting future diagnoses
    p = fn_diagnosis_ratio(o, p, to_diagnose, i)
    
    # Store total number of those diagnosed this time step
    m = update_output(o, p, m, i, "confirmed", ppl[to_diagnose$id, ])
    
    # Isolate newly diagnosed individuals where appropriate
    ppl_next = fn_isolate(p, ppl_next, to_diagnose$id, i)
    
    # ---- Vaccination ----
    
    # Vaccinate priority groups according to p$number_vaccines or p$vaccine_coverage
    ppl_next = fn_vaccinate(p, ppl_next, n_vaccine_group, i)
    
    # Store the number of people vaccinated today and in total
    m = update_output(o, p, m, i, "n_vaccinated",     ppl_next[days_vaccinated == 0, ])
    m = update_output(o, p, m, i, "total_vaccinated", ppl_next[days_vaccinated >= 0, ])
    
    # TODO: Consider time to second dose here...
    
    # Store number of vaccine doses required to achieve this
    n_doses = sum(p$n_doses[ppl_next[days_vaccinated == 0, vaccine_name]])
    m = update_output(o, p, m, i, "n_doses", n_doses)
    
    # ---- Other updates ----
    
    # Move one day closer to next 'event' (see model_flow config file for details)
    ppl_next[!is.na(days_next_event), days_next_event := days_next_event - 1]
    
    # Sanity check: throw an error if days to next event has become negative
    if (nrow(ppl_next[!is.na(days_next_event) & days_next_event < 0, ]) > 0)
      stop("Model error: negative days until next event")
    
    # Update disease and care state then sample duration of next phase
    updates = update_state(o, p, ppl_next, i, m = m)
    
    # Apply updates to ppl array and metric list, m
    ppl_next = updates$ppl
    m        = updates$m
    
    # Increment the number of days infected and infectious, or recovered
    ppl_next[!is.na(days_infected),   days_infected   := days_infected   + 1]
    ppl_next[!is.na(days_infectious), days_infectious := days_infectious + 1]
    
    # Update viral load in next time step for all infected people
    ppl_next = fn_viral_load(p, ppl_next)
    
    # TODO: Not using days_recovered yet, this can be used to quantify waning immunity...
    
    # Also increment the number of days since recovery and/or vaccination
    ppl_next[!is.na(days_recovered),  days_recovered  := days_recovered  + 1]
    ppl_next[!is.na(days_vaccinated), days_vaccinated := days_vaccinated + 1]
    
    # Update immunity state in next time step for all non-infected people
    ppl_next = fn_vaccine_effect(p, ppl_next)
    ppl_next = fn_immunity(p, ppl_next)
    
    # Calculate number in each disease and care state
    for (state in o$disease_states) states[[state]][i] = length(ppl[disease_state == state, id])
    for (state in o$care_states)    states[[state]][i] = length(ppl[care_state    == state, id])
    
    # TODO: Replace dead people? If yes, do that here
    if (o$replace_dead == TRUE)
      stop("Haven't yet implemented the 'replace_dead' feature")
    
    # Update ppl array
    ppl = data.table::copy(ppl_next)
    
    # Update progress bar
    if (verbose == "bar") 
      setTxtProgressBar(pb, i - p$import_date_initial + 1)
  }
  
  # ---- Close up ----
  
  # Close progress bar
  if (verbose == "bar")
    close(pb)
  
  # Prepare final model output
  m = format_output(o, p, m, run_time)
  
  # Produce plots if desired
  if (do_plot) {
    
    # Disease and care states over time
    plot_disease_state(o, p, states, run_time)
    
    # Future testing assumptions
    plot_future_testing(o, p)
  }
  
  return(m)
}

# ---------------------------------------------------------
# Set initial conditionals for people array
# ---------------------------------------------------------
initiate_ppl = function(o, p) {
  
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
    test_date         = numeric(p$population_size),    # Date index for getting tested	
    diagnosed         = logical(p$population_size),    # Whether this individual's current infection has been diagnosed	
    n_infections      = integer(p$population_size),    # Total number of infections experienced per person	
    viral_load        = numeric(p$population_size),    # Viral load for infected individuals	
    variant           = character(p$population_size),  # Variant infected with (not modelling mutations, so this requires importation for each variant)	
    days_next_event   = numeric(p$population_size),    # Days until next event occurs	
    days_recovered    = numeric(p$population_size),    # Number of days since recovery - used to define immunity profile	
    immune_state      = numeric(p$population_size),    # Captures immune state of individual and other relevant comorbidities like blood group, Il6 to Il10 ratio, and neanderthal recombinant genome	
    vaccine_priority  = character(p$population_size),  # Priority group for vaccination (group 1 highest priority)	
    vaccine_accept    = logical(p$population_size),    # Whether this individual with accept vaccination
    days_vaccinated   = numeric(p$population_size),    # Number of days since vaccination - used to determine increasing then waning effect
    vaccine_name      = character(p$population_size),  # Name of vaccine recieved (see data/vaccine_properties.csv)
    vaccine_effect    = numeric(p$population_size))    # Current level of vaccine effect based on which vaccine and when vaccinated
  
  # ---- Set initial values ----
  
  # Sample ages from some distribution
  init_age = sample_vec(x    = o$all_ages, 
                        size = p$population_size, 
                        prob = p$demography, 
                        replace = TRUE)
  
  # Initate basic values in people datatable
  ppl[, id  := 1 : p$population_size]
  ppl[, age := init_age]
  
  # Start by assuming everyone is susceptible
  ppl[, comorbidities     := FALSE]  # Assigned below
  ppl[, healthcare_worker := FALSE]  # Assigned below
  ppl[, disease_state    := "susc"]  # Assume everyone starts susceptible
  ppl[, care_state       := "none"]
  ppl[, prognosis_state  := "none"]
  ppl[, days_infected    := NA]
  ppl[, days_infectious  := NA]
  ppl[, test_date        := NA]
  ppl[, diagnosed        := NA]
  ppl[, n_infections     := 0L]
  ppl[, viral_load       := NA]
  ppl[, variant          := "none"]
  ppl[, days_next_event  := NA]
  ppl[, days_recovered   := NA]
  ppl[, immune_state     := 0]  # Assume everyone starts with no immunity
  ppl[, vaccine_priority := "none"]
  ppl[, days_vaccinated  := NA]
  ppl[, vaccine_accept   := FALSE]
  ppl[, vaccine_name     := "none"]
  ppl[, vaccine_effect   := 0]
  
  # ---- Assign risk groups ----
  
  # TODO: Generalise risk groups and remove this hardcoding
  
  # Healthcare workers
  ppl[age >= 18 & age <= 65, 
      healthcare_worker := sample_vec(c(FALSE, TRUE), size = .N, replace = TRUE, 
                                      prob = c(1 - p$risk_group_ratio$healthcare_workers, 
                                               p$risk_group_ratio$healthcare_workers))]
  
  # People with comorbidities (under 65)
  ppl[age >= 18 & age <= 65, 
      comorbidities := sample_vec(c(FALSE, TRUE), size = .N, replace = TRUE, 
                                  prob = c(1 - p$risk_group_ratio$young_comorbidities, 
                                           p$risk_group_ratio$young_comorbidities))]
  
  # People with comorbidities (over 65)
  ppl[age > 65, 
      comorbidities := sample_vec(c(FALSE, TRUE), size = .N, replace = TRUE, 
                                  prob = c(1 - p$risk_group_ratio$old_comorbidities, 
                                           p$risk_group_ratio$old_comorbidities))]
  
  # ---- Set vaccination priority policy ----
  
  # Iterate through the different priority groups
  # 
  # NOTE: We do last at first so any people in multiple groups are placed in the highest possible group
  for (group in rev(names(o$vaccine_priority))) {
    
    # String representing condition(s) to be in this priority group
    group_condition = o$vaccine_priority[[group]]
    
    # Evaluate this string and set priority for any that satisfy
    ppl[eval(parse(text = group_condition)), vaccine_priority := group]
  }
  
  # Again iterate through priority groups - this time to set vaccine acceptance probabilities
  #
  # NOTE: We do this after all groups are formalised to prevent double counting
  for (group in names(o$vaccine_priority)) {
    
    # Proportion of this group that would accept a vaccine
    group_acceptance = as.numeric(p$vaccine_acceptance[group])
    
    # Sample logical and populate vaccine_accept variable
    ppl[vaccine_priority == group, 
        vaccine_accept := sample_vec(c(TRUE, FALSE), size = .N, replace = TRUE, 
                                     prob = c(group_acceptance, 1 - group_acceptance))]
  }
  
  return(ppl)
}

# ---------------------------------------------------------
# Initially infect a subset of people to kick start epidemic
# ---------------------------------------------------------
initiate_epidemic = function(o, p, ppl) {
  
  # Initiate with a few random people being infected (see below code for more nuance here)
  init_infected = base::sample.int(p$population_size, size = p$import_initial)
  
  # Set basic properties for these initially infected people
  ppl[init_infected, n_infections := 1L]             # Count number of infections experienced
  ppl[init_infected, days_infected := 1]             # Initiate days since this infection
  ppl[init_infected, variant := o$variants$primary]  # They all start with the primary variant
  ppl[init_infected, diagnosed := FALSE]             # Start undiagnosed
  
  # Sample a prognosis for these initially infected people
  ppl[init_infected, prognosis_state := fn_prognosis(o, p, ppl[init_infected, ], p$import_date_initial)]
  
  # We'll force initiate with presymptomatic infection by setting a trivial latent phase
  ppl[init_infected, disease_state   := "latent"]
  ppl[init_infected, days_next_event := 0]
  
  # Update disease state and sample duration of next phase
  ppl = update_state(o, p, ppl, p$import_date_initial)$ppl
  
  # Initiate infectiousness immediately to kick start the epidemic
  ppl[init_infected, days_infectious := 1]           
  
  # Quantify initial viral load - used to quantify infectiousness
  ppl = fn_viral_load(p, ppl)
  
}

# ---------------------------------------------------------
# Set initial conditionals for people array
# ---------------------------------------------------------
create_network = function(o, p, ppl, do_plot, verbose) {
  
  # ---- Simple random network ----
  
  # Simple random network with no particular structure
  if (o$network_structure == "random") {
    
    if (verbose != "none")
      message("  > Creating simple random contact network")
    
    # Number of nodes (people) and edges (contacts)
    n_edges = round(p$population_size * p$contacts / 2)
    
    # Generate simple network
    network = play_erdos_renyi(n = p$population_size, m = n_edges, directed = FALSE)
    
    # Format edges into a datatable - this is what we really need
    network_df = network %>% activate(edges) %>% as.data.table()
    
    # Repeat pair-wise so all contacts are double directed
    reverse_df = data.table(from = network_df$to,
                            to   = network_df$from)
    
    # Bind into single datatable
    all_contacts = rbind(network_df, reverse_df)
  }
  
  # ---- Age-structured network ----
  
  # Age structured network informed by socialmixr contact matrices
  if (o$network_structure == "age") {
    
    if (verbose != "none")
      message("  > Creating age-structured contact network")
    
    # Create age groups, depends on initiate_ppl
    age_limits <- sort(unique(ppl$age))
    
    # Create contact matrix from socialmixr
    #
    # NOTE: The contact_matrix function complains when age_limits are not consistent with it's data, however it performs
    #       interpolation between ages to prevent data loss. We therefore justify the use of suppressWarnings here.
    polymod_data = suppressWarnings(contact_matrix(survey     = polymod, 
                                                   countries  = o$contact_matrix_countries, 
                                                   age.limits = age_limits, 
                                                   symmetric  = TRUE, 
                                                   quiet      = TRUE))
    
    # Convert matrix to datatable
    contact_matrix_data = data.table(polymod_data$matrix) 
    
    # Extract the age groups that were created by socialmixr, usually the oldest individuals are lumped together
    contact_age_groups <- 1 : length(names(contact_matrix_data))
    
    # Temporarily append an age_group variable in ppl, to later sample ID's according to age group
    ppl[, age_group := min(which(age_limits == age[1]), length(contact_age_groups)), by = list(age)]
    
    # How many individuals per age group have we created? Using data.tables .N variable and order by age_group
    created_demog <- ppl[order(age_group), .N, by = list(age_group)]$N
    
    # Multiply the contact matrix by the demography counts and unlist to sample later. Unlist works by column
    contact_probs <- unlist(contact_matrix_data[, Map("*", .SD, created_demog)], use.names = FALSE)
    contact_probs[is.na(contact_probs)] <- 0
    
    # Create ego_cell and alter_cell vectors to get cell identities
    ego_cell   <- rep(1:length(contact_age_groups), each  = length(contact_age_groups))
    alter_cell <- rep(1:length(contact_age_groups), times = length(contact_age_groups))
    
    # Pre-allocate edgelist
    n_contacts   = round(p$population_size * p$contacts / 2)
    all_contacts = data.table(ego_age_group   = integer(n_contacts), 
                              alter_age_group = integer(n_contacts), 
                              from            = integer(n_contacts), 
                              to              = integer(n_contacts))
    
    # Sample contacts
    realized_contacts <- sample(x = 1 : length(contact_probs), 
                                size = n_contacts, 
                                prob = contact_probs, 
                                replace = TRUE)
    
    # Assign the age groups per contacts to the all_contacts edgelist
    all_contacts[, `:=` (ego_age_group   = ego_cell[realized_contacts], 
                         alter_age_group = alter_cell[realized_contacts])]
    
    # Now sample from ppl$age_groups, using data.tables by = as well as .N (for outgoing and incoming contacts)
    all_contacts[, from := sample(ppl[age_group == ego_age_group[1],   id], size = .N, replace = TRUE), by = list(ego_age_group)]
    all_contacts[, to   := sample(ppl[age_group == alter_age_group[1], id], size = .N, replace = TRUE), by = list(alter_age_group)]
    
    # Remove self-loops and additional variables
    all_contacts <- all_contacts[from != to, c("from", "to")]
    
    # Mirror to reflect two sidedness of contacts
    all_contacts_reverse <- data.table(from = all_contacts$to, to = all_contacts$from)
    all_contacts <- rbindlist(list(all_contacts, all_contacts_reverse), use.names = FALSE, fill = FALSE)
    
    # Finally, remove temporary 'age groups' variable
    ppl[, age_group := NULL]
  }
  
  # ---- Perform checks and calculate age correction factor ----
  
  # Throw an error if no edge list created - likely because of invalid network_structure
  if (!exists("all_contacts"))
    stop("Network structure '", o$network_structure, "' not recognised")
  
  # Throw an error if average number of contacts is drastically wrong
  mean_contacts = nrow(all_contacts) / p$population_size
  if (abs(mean_contacts - p$contacts) > 0.1)
    stop("Requested an average of ", p$contacts, " contacts but generated ", round(mean_contacts, 3))
  
  # Number of contacts per person 
  count_contact = table(all_contacts$from)
  
  # Seperate into IDs and number of contacts
  id = as.numeric(names(unlist(count_contact)))
  n_contacts = as.numeric(unlist(count_contact))
  
  # Breaks to create age bins
  # 
  # NOTE: The 10 year age bins relates the the 10-year prognosis data we use (see parameters.R)
  age_breaks = seq(0, max(o$all_ages) + 1, by = 10)
  
  # Total number of contacts per age group
  age_df = data.table(id = id, n_contacts = n_contacts) %>%
    full_join(ppl[, .(id, age)], by = "id") %>%
    mutate(age_group = cut(age, age_breaks, include.lowest = TRUE)) %>%
    group_by(age_group) %>%
    summarise(n_contacts = sum(n_contacts, na.rm = TRUE))
  
  # Relative probability of a contact (and therefore of infection - when all other factors are equal)
  age_contacts = age_df$n_contacts / sum(age_df$n_contacts)
  
  # Use this to correct prognosis probabilities to achieve required disease and death age-distributions
  age_correction = (1 - age_df$n_contacts / mean(age_df$n_contacts)) + 1  # See fn_prognosis to see this in action
  
  # ---- Diagnostic plots and outputs ----
  
  # Produce network visualisations if desired
  if (do_plot == TRUE) {
    
    message("   - Plotting contact network properties")
    
    # Append ages of contacts to edge list
    plot_contacts = all_contacts %>% 
      mutate(age1 = ppl[from, age], 
             age2 = ppl[to, age])
    
    # Network figure 1) Series of network-related properties
    plot_network_properties(o, p, plot_contacts)
    
    # Network figure 2) Age matrix of contact density per age (single age bins)
    plot_network_matrix(o, plot_contacts)
  }
  
  # Append output details to a list
  network = list(all_contacts   = all_contacts, 
                 age_contacts   = age_contacts, 
                 age_correction = age_correction)
  
  return(network)
}

# ---------------------------------------------------------
# Construct model output dataframe 
# ---------------------------------------------------------
preallocate_output = function(o, run_time) {
  
  # Preallocate output
  m = list()
  
  # Easy access metric info datatable
  metrics = o$all_metrics
  
  # Loop through all types of grouping
  for (this_grouping in unique(o$metric_groupings)) {
    
    # Stratification of this grouping (can be altered by user)
    stratify_grouping = o[[paste0("count_", this_grouping)]]
    
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
                            date     = 1 : run_time, 
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
update_output = function(o, p, m, date_idx, update_metric, to_count) {
  
  # Loop through all groupings defined for this metric
  metric_groupings = o$metric_groupings[names(o$metric_groupings) == update_metric]
  for (grouping in metric_groupings) {
    
    # Apply reporting delay (if any) for this metric
    report_idx = date_idx + p$reporting_delays[update_metric]
    
    # For metrics with no possible disaggregation, append single value
    if (grouping == "na") {
      m[metric == update_metric & date == report_idx, value := to_count]
      
      # For metrics we have chosen not to disaggregate, count all values
    } else if (grouping == "none") {
      m[metric == update_metric & date == report_idx, value := nrow(to_count)]
      
    } else {  # Otherwise we want to count each grouping...
      
      # Count number of occurrences of each group
      group_count = table(to_count[, ..grouping])
      
      # Insert these values into output datatable
      m[metric == update_metric &         # This metric
          date == report_idx &            # This time step
          group %in% names(group_count),  # This subset of 'groups'
        value := group_count]
    }
  }
  
  # Check the update: m[metric == update_metric & date == report_idx, ]
  
  return(m)
}

# ---------------------------------------------------------
# Prepare final model output
# ---------------------------------------------------------
format_output = function(o, p, m, run_time) {
  
  # Convert number of cases per variant to variant prevalence
  m[metric == "variant_prevalence", value := 100 * value / max(sum(value), 1), by = "date"]
  
  # Append effect of infectiousness per contact due to seasonality
  m[metric == "seasonality", value := p$seasonality[1 : run_time] * p$beta * 100]
  
  # Extract contact reduction level due to NPIs over time
  npi_level = p$npi_level[1 : run_time]
  
  # Append total NPI effect and unscaled Oxford Stringency Index
  m[metric == "npi_effect", value := npi_level * p$npi_adherence[1 : run_time]]
  m[metric == "osi_level",  value := npi_level * 100 / p$npi_scaler]
  
  # Metrics for which population scaling is to be applied
  scale_metrics = o$all_metrics$metric[o$all_metrics$scaled]
  
  # Apply the population scaling for these metrics
  m[metric %in% scale_metrics, value := value * p$population_scaler]
  
  # Convert date indices to dates
  m[, date := o$dates_all[date]]
  
  # Reorder columns to be consistant with data format
  setcolorder(m, qc(date, metric, grouping, group, value))
  
  return(m)
}

# ---------------------------------------------------------
# Update disease and care states for all relevant individuals
# ---------------------------------------------------------
update_state = function(o, p, ppl, date_idx, m = NULL) {
  
  # IDs of people that need to be updated
  update_id = ppl[days_next_event == 0, id]
  
  # Skip this process if no one to update
  if (length(update_id) > 0) {
    
    # Join with model flows dataframe to determine next state
    update_df = left_join(ppl[update_id], p$model_flows, 
                          by = qc(disease_state, care_state, prognosis_state))
    
    # Count and store values of each metric in model output
    for (update_metric in unique(na.omit(update_df$metric)))
      m = update_output(o, p, m, date_idx, update_metric, update_df[metric == update_metric, ])
    
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
      ppl[reset_id, days_infected   := NA]
      ppl[reset_id, days_infectious := NA]
      ppl[reset_id, test_date       := NA]
      ppl[reset_id, diagnosed       := NA]
      ppl[reset_id, viral_load      := NA]
      ppl[reset_id, variant         := "none"]  # We may want to record this for those that die
      ppl[reset_id, prognosis_state := "none"]
      
      # Some extra things required for those who have recovered
      recover_id = update_df[metric == "recovered", id]
      
      # TODO: Waning acquired immunity over time - consider integrating into fn_immunity
      # TODO: If number of infections is above a certain cap, set immune_state to 1 with no waning
      
      # Signal start of recovery period by starting 'days infectious' count
      ppl[recover_id, immune_state := pmax(p$acquired_immunity, immune_state)]
      ppl[recover_id, days_recovered := 0]  # Incremented at end of daily loop
    }
    
    # ---- Sample duration until next event ----
    
    # Loop through the different duration functions
    for (duration in unique(na.omit(update_df$next_duration))) {
      
      # Call duration function for the relevant individuals
      duration_id  = update_df[next_duration == duration, id]
      duration_val = p[[duration]](duration_id)
      
      # Store the time until next event
      ppl[duration_id, days_next_event := duration_val]
    }
    
    # ---- Initaite infectious period and symptom onset ----
    
    # Newly infectious people are those going INTO the pre-symptomatic phase
    new_infectious_id = update_df[next_disease == "presym", id]
    
    # Signal start of infectious period by starting 'days infectious' count
    ppl[new_infectious_id, days_infectious := 0]  # Incremented at end of daily loop
    
    # Newly symptomatic people are those coming FROM pre-symptoms phase
    new_symptoms_id = update_df[disease_state == "presym", id]
    
    # Define date for any potential test after symptom onset
    ppl[new_symptoms_id, test_date := date_idx + round(p$diagnosis_delay)]
  }
  
  return(list(ppl = ppl, m = m))
}

# ---------------------------------------------------------
# Initiate new infections for newly infected individuals
# ---------------------------------------------------------
fn_begin_infection = function(o, p, ppl, id, apply_variant, date_idx) {
  
  # We signal this new infection by starting 'days infected' count
  ppl[id, days_infected := 0]  # Incremented at end of daily loop
  
  # Increment number of infections experienced thus far
  ppl[id, n_infections := n_infections + 1L]
  
  # Sanity check: we need a variant defined for each person
  if (!(length(apply_variant) %in% c(length(id), 1)))
    stop("Inconsistent number of variants defined")
  
  # Apply the precalculated variant for this individial
  ppl[id, variant := apply_variant]
  
  # All newly infected people go to the latent phase
  ppl[id, disease_state := "latent"]
  
  # None of these people start diagnosed
  ppl[id, diagnosed := FALSE]
  
  # Define how log this latency phase will last
  ppl[id, days_next_event := p$duration_latency(id)]
  
  # Sample a prognosis for each infected individual
  ppl[id, prognosis_state := fn_prognosis(o, p, ppl[id, ], date_idx)]
  
  return(ppl)
}

# ---------------------------------------------------------
# Generate prognosis state for newly infected individuals
# ---------------------------------------------------------
fn_prognosis = function(o, p, new_infected, date_idx) {
  
  # ---- Construct cumulative probability matrix ----
  
  # Proportion of all cases that are asymptomatic
  p_asym = rep(p$proportion_asymptomatic, length(o$age_groups))
  
  # Proportion of symptomatic cases that are severe and require hospitalisation (by age group)
  # 
  # NOTE: This considers (potentially) less contacts in older ages and a calibration scaling factor
  symptom_severe = pmin(p$symptom_severe_age * p$severe_factor * p$age_correction, 1)
  
  # plot_df = data.frame(x = 1 : length(p$symptom_severe_age),
  #                      y1 = pmin(p$symptom_severe_age * p$severe_factor, 1),
  #                      y2 = pmin(p$symptom_severe_age * p$severe_factor * p$age_correction, 1))
  # 
  # g = ggplot(plot_df, aes(x = x)) +
  #   geom_line(aes(y = y1), colour = "blue") +
  #   geom_line(aes(y = y2), colour = "green")
  
  # Proportion of mild cases: symptomatic, non-severe
  p_mild   = (1 - p_asym) * (1 - symptom_severe)
  p_severe = (1 - p_asym) * symptom_severe
  
  # Proportion of severe cases that will become critical and require care in ICU (by age group)
  severe_critical = pmin(p$severe_critical_age * p$critical_factor * 
                           p$improved_care_scaler[date_idx], 1)  # With calibration factor(s)
  
  # Proportion of cases that are severe that will / will not become critical
  p_critical     = p_severe * severe_critical
  p_non_critical = p_severe * (1 - severe_critical)
  
  # TEMP: Additional correction factors for the most elderly age group
  elderly_idx = length(o$age_groups)
  
  # TEMP: Those in oldest age group are less likely to seek care
  seek_hospital = rep(p$seek_hospital, elderly_idx)
  seek_hospital[elderly_idx] = seek_hospital[elderly_idx] * p$elderly_seek_factor
  
  # Split these by hospital-seeking behaviour 
  p_seek_S   = p_non_critical * seek_hospital
  p_noseek_S = p_non_critical * (1 - seek_hospital)
  
  # Of the critical cases admitted to ICU, who will recover and who will die
  p_seek_C = p_critical * seek_hospital * (1 - p$critical_death_icu * p$death_critical_age)
  p_seek_D = p_critical * seek_hospital * p$critical_death_icu * p$death_critical_age
  
  # Of the critical cases outside of the care system, who will recover and who will die
  p_noseek_C = p_critical * (1 - seek_hospital) * (1 - p$critical_death_non_icu * p$death_critical_age)
  p_noseek_D = p_critical * (1 - seek_hospital) * p$critical_death_non_icu * p$death_critical_age
  
  # death_dist = (p_seek_D + p_noseek_D) * p$age_contacts
  # death_norm = death_dist  / sum(death_dist)
  # 
  # plot_df = data.frame(x = 1 : length(death_norm), 
  #                      y1 = p$death_age, 
  #                      y2 = death_norm)
  # 
  # g = ggplot(plot_df, aes(x = x)) +
  #   geom_line(aes(y = y1), colour = "black") +
  #   geom_line(aes(y = y2), colour = "green")
  
  # Re-aggregate these probabilities to check we haven't made silly mistakes
  p_seek   = p_seek_S   + p_seek_C   + p_seek_D
  p_noseek = p_noseek_S + p_noseek_C + p_noseek_D
  
  # Sanity check: Severe people should either seek or not seek hospital care
  if (any(abs(p_seek + p_noseek - p_severe) > 1e-6))
    stop("Hospial seeking probabilities must align for all age groups")
  
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
  
  # Indices of death probabilities - may need to adjust death probabilities for certain variants
  death_idx = grepl("^dead", prognosis_names)
  
  # ---- Sample to determine age-related prognosis ----
  
  # Preallocate vector of prognoses
  prognosis = rep(NA, nrow(new_infected))
  
  # NOTE: new_infected is just a slice of ppl and include all variables
  
  # Age and comorbidity status of newly infected individuals
  age     = new_infected[, age]
  variant = new_infected[, variant]
  
  # TODO: This is ok for now, but will want to improve on this approach...
  
  # Represent increase risk of comorbidities by shifting up one age class
  age[new_infected$comorbidities] <- age[new_infected$comorbidities] + 10
  age <- pmin(age, max(o$all_ages))
  
  # Convert age to age group index to reference
  age_bin = floor(age / 10) + 1  # TODO: Avoid this /10 hardcoding
  
  # Repeat calculation for all currently circulating variants
  for (this_variant in unique(variant)) {
    
    # Reset prognosis matrix
    prognosis_variant = prognosis_matrix
    
    # Increase probability of death for particular variants (see p$variant_severity)
    prognosis_variant[, death_idx] = prognosis_variant[, death_idx] * p$variant_severity[[this_variant]]
    
    # Normalise by reducing probability of being asymptomatic
    prognosis_variant[, 1] = 1 - rowsums(prognosis_variant[, -1])
    
    # Cumulative sum the probabilities to generate a matrix that can be sampled
    prognosis_cumsum = rowCumsums(prognosis_variant)
    
    # Check saniy check that probabilities still sum to 1
    if (any(abs(prognosis_cumsum[, n_prognoses] - 1) > 1e6))
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
  
  # This only comes into play if we assume vaccine only prevents disease
  if (p$vaccine_effect %in% c("severity", "mixed")) {
    
    # There is an increased probability (on top of standard proportion_asymptomatic) of being asym
    increased_prob_asym = new_infected[, vaccine_effect] * (1 - p$proportion_asymptomatic)
    
    # Skip if none of these people have been vaccinated
    if (sum(increased_prob_asym) > 0) {
      
      # TODO: Need to do some more tests to check this does what we expect
      
      # Draw random number to determine if this person should be asymptomatic due to vaccine effect
      asym_vaccine_effect = runif(nrow(new_infected)) < increased_prob_asym * p$disease_prevention_scaler
      
      # Reassign prognosis to asymptomatic disease
      prognosis[asym_vaccine_effect] = "asym"
    }
  }
  
  return(prognosis)
}

# ---------------------------------------------------------
# Test and diagnose individuals
# ---------------------------------------------------------
fn_test_diagnose = function(p, ppl, n_infected, date_idx) {
  
  # For all past dates, we know how many people are to be diagnosed from data
  n_diagnose = p$number_diagnosed[date_idx]
  
  # For the future, assume similar proportion of infected cases get diagnosed
  if (is.na(n_diagnose))
    n_diagnose = round(n_infected * p$future_ratio[date_idx])
  
  # There is an issue if we still have an NA
  if (is.na(n_diagnose))
    stop("Error calculating future diagnosis ratio")
  
  # Calculate diagnosis ratio - this is what we'll use to generate future assumptions
  diagnosed_ratio = n_diagnose / n_infected
  diagnosed_ratio[!is.finite(diagnosed_ratio)] = NA
  
  # All those with severe-enough disease that they are guarantee to get diagnosed today
  diagnosed_id = ppl[test_date == date_idx & disease_state %in% p$diagnosis_guarantee, id]
  
  # Number of 'other' people (ie those without severe) to be diagnosed today
  n_diagnose_others = max(0, n_diagnose - length(diagnosed_id))
  
  # Skip this step if there is no one else to diagnose
  if (n_diagnose_others > 0) {
    
    # All others who could potentially be tested and diagnosed today
    diagnose_potential = ppl[test_date == date_idx & !(id %in% diagnosed_id), ]
    
    # Order who could potentially be tested by disease state (severe gets priority)
    priority_order = order(p$diagnosis_priority[diagnose_potential$disease_state])
    
    # Take n_diagnose people from the top of the priority list
    priority_select = priority_order[seq_len(min(n_diagnose, length(priority_order)))]
    
    # Concatenate with any guaranteed diagnoses
    diagnosed_id = sort(c(diagnosed_id, diagnose_potential[priority_select, id]))
  }
  
  # ---- Contact tracing ----
  
  # TODO: Re-implement contact-tracing
  
  # Output IDs of those diagnosed along with mild : severe ratio
  to_diagnose = list(id = diagnosed_id, ratio = diagnosed_ratio)
  
  return(to_diagnose)
}

# ---------------------------------------------------------
# Calculate diagnoses per infected case probability for future projections
# ---------------------------------------------------------
fn_diagnosis_ratio = function(o, p, to_diagnose, date_idx) {
  
  # This function has three jobs:
  #  1) Preallocate vector for diagnoses per infected case
  #  2) Store the ratio in said vector at each time point (ratio calculated in fn_test_diagnose)
  #  3) When the data runs out, use the ratio vector to project future ratio to ensure smooth projections
  
  # On first pass we need to initiate vector
  if (is.null(p$diagnosis_ratio))
    p$diagnosis_ratio = rep(0, o$n_dates)
  
  # Store ratio calculate in fn_test_diagnose
  p$diagnosis_ratio[date_idx] = to_diagnose$ratio
  
  # ---- Fit future ratio based on model outcomes ---
  
  # Number of confirmed cases in the past - taken directly from data
  diagnosed_data = p$number_diagnosed[date_idx + 1]
  
  # When data runs out we need to fit and project - only do this once
  if (is.na(diagnosed_data) && is.null(p$future_ratio)) {
    
    # If vector completely trivial populate with some nominally low value
    if (sum(p$diagnosis_ratio, na.rm = TRUE) == 0)
      p$diagnosis_ratio[1 : date_idx] = 1e-6
    
    # Time horizon for fitting severe : non-severe diagnosis ratio
    time_start = which(p$diagnosis_ratio[1 : date_idx] > 0)[1]
    time_seq   = time_start : date_idx
    
    # Impute missing points through linear interpolation
    diagnosis_ratio = na_interpolation(p$diagnosis_ratio[time_seq])
    
    # Method 1: fit a straight line
    if (o$future_diagnoses == "linear_model") {
      
      # Fit a simple linear model to this data
      linear_model  = stats::lm(formula = diagnosis_ratio ~ time_seq)$coefficients
      extrapolation = linear_model[1] + linear_model[2] * (1 : o$n_dates)
    }
    
    # Method 2: project a constant (last 30 days)
    if (o$future_diagnoses == "constant")
      extrapolation = rep(mean(tail(diagnosis_ratio, 14)), o$n_dates)
    
    # We'll use this for ratio for forward projections
    p$future_ratio = pmax(extrapolation, 0)
    p$future_ratio[1 : date_idx] = NA
  }
  
  return(p) 
}

# ---------------------------------------------------------
# Isolate newly diagnosed individuals where appropriate
# ---------------------------------------------------------
fn_isolate = function(p, ppl, diagnosed_id, date_idx) {
  
  # TODO: Do we want a delay before isolation?
  
  # First things first, mark individuals as diagnosed
  ppl[diagnosed_id, diagnosed := TRUE]
  
  # Reset any now redundant test_date values
  ppl[test_date == date_idx, test_date := NA]
  
  # Once diagnosed, only those with certain prognoses can go into isolation (see model_flows)
  #
  # INTERPRETATION: You would not isolate if you have severe disease
  can_isolate_id = ppl[id %in% diagnosed_id & prognosis_state %in% p$iso_prognosis, id]
  
  # Sample number of people due to isolate (according to isolation_probability)
  n_isolate  = length(can_isolate_id) * p$isolation_probability
  isolate_id = sample_vec(can_isolate_id, round(n_isolate))
  
  # Mark these people as isolated and keep them there for isolation_duration days
  ppl[isolate_id, care_state := "iso"]
  ppl[isolate_id, days_next_event := round(p$isolation_duration)]
  
  return(ppl)
}

# ---------------------------------------------------------
# Vaccinate eligible people (according to number_vaccines or vaccine_coverage)
# ---------------------------------------------------------
fn_vaccinate = function(p, ppl, n_vaccine_group, date_idx) {
  
  # Are there any people to vaccinate (by number or coverage)
  to_vaccinate_number   = sum(p$number_vaccines[date_idx])
  to_vaccinate_coverage = sum(p$vaccine_coverage[date_idx, ])
  
  # Skip this step is there is no one to vaccinate
  if (to_vaccinate_number >= 1 || to_vaccinate_coverage > 0) {
    
    # All people eligible for vaccination at this time step (regardless of priority group)
    vaccine_eligible = ppl[eval(parse(text = p$vaccine_criteria)), ]
    
    # Exclude any unaccepting of a vaccine and order by vaccine priority group
    vaccine_queue = vaccine_eligible[vaccine_accept == TRUE][order(vaccine_priority)]
    
    # ---- Distribute vaccines by number ----
    
    # Check for vaccination by number
    if (to_vaccinate_number >= 1) {
      
      # Take the top n from vaccine queue
      vaccinate_id = vaccine_queue[seq_len(to_vaccinate_number), id]
      
      # Start these people off on their vaccination journey
      #
      # TODO: Currently just using default vaccine here, should be able to choose
      ppl[vaccinate_id, days_vaccinated := 0]
      ppl[vaccinate_id, vaccine_name := p$vaccine_default]
    }
    
    # ---- Distribute vaccines by group coverage ----
    
    # Check for vaccination by group coverage
    if (to_vaccinate_coverage > 0) {
      
      # Iterate through the vaccine priority groups
      for (group in colnames(p$vaccine_coverage)) {
        
        # Define which vaccine this group are due to recieve
        vaccine_type = as.character(p$vaccine_type[group])
        
        # Number of people we want vaccinated in this group in this time step
        n_group   = n_vaccine_group[vaccine_priority == group, N]
        n_desired = round(p$vaccine_coverage[date_idx, group] * n_group)
        
        # Number of people currently vaccinated in this group
        n_current = nrow(ppl[vaccine_priority == group & !is.na(days_vaccinated), ])
        
        # All eligible for vaccination in this goup
        group_queue = vaccine_queue[vaccine_priority == group]
        
        # How many we will vacinate this time step - bounded above by number available
        n_vaccinate = min(max(n_desired - n_current, 0), nrow(group_queue))
        
        # Vaccinate the first chunk of people - randomness is inherent as table is not ordered
        vaccinate_id = group_queue[seq_len(n_vaccinate), id]
        
        # Start these people off on their vaccination journey
        ppl[vaccinate_id, days_vaccinated := 0]
        ppl[vaccinate_id, vaccine_name := vaccine_type]
      }
    }
  }
  
  return(ppl)
}

# ---------------------------------------------------------
# Determine level of vaccine effect based on what vaccine and when
# ---------------------------------------------------------
fn_vaccine_effect = function(p, ppl) {
  
  # ID of all vaccinated individuals
  vaccinated_id = ppl[days_vaccinated > 0, id]
  
  # Skip this process if no-one vaccinated
  if (length(vaccinated_id) > 0) {
    
    # Loop through the different vaccines
    vaccine_names = ppl[vaccinated_id, vaccine_name]
    
    # Loop through to represent the unique properties of each vaccine
    for (vaccine_name in unique(vaccine_names)) {
      
      # IDs of people receiving this vaccine
      this_vaccine_id = vaccinated_id[vaccine_names == vaccine_name]
      
      # Days since vaccination for all these people
      this_vaccine_days = ppl[this_vaccine_id, days_vaccinated]
      
      # What this means in terms of 'effect' for this vaccine
      this_vaccine_effect = p$vaccine_profile[[vaccine_name]][this_vaccine_days]
      
      # Update variable in ppl dataframe for these vaccinated individuals
      ppl[this_vaccine_id, vaccine_effect := this_vaccine_effect]
    }
  }
  
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
fn_immunity = function(p, ppl) {
  
  # ---- Immunity due to infection ----
  
  # TODO: Incorporate acquired immunity here too (see days_recovered in ppl)
  
  # ---- Immunity due to vaccination ----
  
  # This only comes into play if we assume vaccine blocks transmission
  if (p$vaccine_effect %in% c("transmission", "mixed")) {
    
    # ID of all vaccinated individuals
    vaccinated_id = ppl[days_vaccinated > 0, id]
    
    # Update immunity variable in ppl dataframe for these vaccinated individuals
    ppl[vaccinated_id, immune_state := 
          pmax(vaccine_effect * p$transmission_blocking_scaler, immune_state)]
  }
  
  return(ppl)
}

