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

# ---- Steering committee meeting: 9th March 2021 ----

# OSI notes:
#  - OSI @ end of Feb: 0.6389
#  - OSI @ 1st March: 0.583

# Vaccine notes:
#  - Feasible to vaccinate 25k people a day (50k doses a day)
#  - Upper limit assumed to be 50k people a day based on BAG feedback (100k doses a day)

# Scenario 0a) No relax after March 1, default vaccine
s0a = function() {
  
  # NPI exit strategy
  npi_exit = list(
    npi0 = list(name   = "No relax after March 1", 
                phase1 = list(start = "2021-03-22", new_normal_rel = 1)))
  
  # Vaccine rollout defined by number_vaccines
  roll_out = list(
    defvax = list(name   = "Proposed vaccine rollout", 
                  group1 = list(coverage = 0.0, start = "2021-01-01", duration = 0)))
  
  # Model parameters to be altered
  alter_params = c(paste0("p$number_vaccines[o$dates_all %in% seq(max(o$dates_data) + 1, format_date('2021-04-01'), by = 'day')] = ", 
                          "round(seq(6000, 25000, length.out = as.numeric(format_date('2021-04-01') - max(o$dates_data))))"), 
                   "p$number_vaccines[o$dates_all > '2021-04-01'] = 25000", 
                   "p$vaccine_effect = 'mixed'", 
                   "p$transmission_blocking_scaler = 0.8",
                   "p$disease_prevention_scaler    = 0.8",
                   "p$variant_infectivity['B117']  = 1.6")
  
  # Format strategy details into single list
  details = list(roll_out = roll_out, npi_exit = npi_exit, alter_params = alter_params)
  
  return(details)
}

# Scenario 0b) No relax after March 1, fast vaccine
s0b = function() {
  
  # NPI exit strategy
  npi_exit = list(
    npi0 = list(name   = "No relax after March 1", 
                phase1 = list(start = "2021-03-22", new_normal_rel = 1)))
  
  # Vaccine rollout defined by number_vaccines
  roll_out = list(
    fastvax = list(name   = "Fast vaccine rollout", 
                   group1 = list(coverage = 0.0, start = "2021-01-01", duration = 0)))
  
  # Model parameters to be altered
  alter_params = c(paste0("p$number_vaccines[o$dates_all %in% seq(max(o$dates_data) + 1, format_date('2021-04-01'), by = 'day')] = ", 
                          "round(seq(6000, 50000, length.out = as.numeric(format_date('2021-04-01') - max(o$dates_data))))"), 
                   "p$number_vaccines[o$dates_all > '2021-04-01'] = 50000", 
                   "p$vaccine_effect = 'mixed'", 
                   "p$transmission_blocking_scaler = 0.8",
                   "p$disease_prevention_scaler    = 0.8",
                   "p$variant_infectivity['B117']  = 1.6")
  
  # Format strategy details into single list
  details = list(roll_out = roll_out, npi_exit = npi_exit, alter_params = alter_params)
  
  return(details)
}

# Scenario 1a) March 22nd only, default vaccine
s1a = function() {
  
  # NPI exit strategy
  npi_exit = list(
    npi1 = list(name   = "Proposed relax March 22", 
                phase1 = list(start = "2021-03-22", new_normal_osi = 0.535)))
  
  # Vaccine rollout defined by number_vaccines
  roll_out = list(
    defvax = list(name   = "Proposed vaccine rollout", 
                  group1 = list(coverage = 0.0, start = "2021-01-01", duration = 0)))
  
  # Model parameters to be altered
  alter_params = c(paste0("p$number_vaccines[o$dates_all %in% seq(max(o$dates_data) + 1, format_date('2021-04-01'), by = 'day')] = ", 
                          "round(seq(6000, 25000, length.out = as.numeric(format_date('2021-04-01') - max(o$dates_data))))"), 
                   "p$number_vaccines[o$dates_all > '2021-04-01'] = 25000", 
                   "p$vaccine_effect = 'mixed'", 
                   "p$transmission_blocking_scaler = 0.8",
                   "p$disease_prevention_scaler    = 0.8",
                   "p$variant_infectivity['B117']  = 1.6")
  
  # Format strategy details into single list
  details = list(roll_out = roll_out, npi_exit = npi_exit, alter_params = alter_params)
  
  return(details)
}

# Scenario 1b) March 22nd only, fast vaccine
s1b = function() {
  
  # NPI exit strategy
  npi_exit = list(
    npi1 = list(name   = "Proposed relax March 22", 
                phase1 = list(start = "2021-03-22", new_normal_osi = 0.535)))
  
  # Vaccine rollout defined by number_vaccines
  roll_out = list(
    fastvax = list(name   = "Fast vaccine rollout", 
                   group1 = list(coverage = 0.0, start = "2021-01-01", duration = 0)))
  
  # Model parameters to be altered
  alter_params = c(paste0("p$number_vaccines[o$dates_all %in% seq(max(o$dates_data) + 1, format_date('2021-04-01'), by = 'day')] = ", 
                          "round(seq(6000, 50000, length.out = as.numeric(format_date('2021-04-01') - max(o$dates_data))))"), 
                   "p$number_vaccines[o$dates_all > '2021-04-01'] = 50000", 
                   "p$vaccine_effect = 'mixed'", 
                   "p$transmission_blocking_scaler = 0.8",
                   "p$disease_prevention_scaler    = 0.8",
                   "p$variant_infectivity['B117']  = 1.6")
  
  # Format strategy details into single list
  details = list(roll_out = roll_out, npi_exit = npi_exit, alter_params = alter_params)
  
  return(details)
}

# Scenario 2a) Keep relaxing, default vaccine
s2a = function() {
  
  # NPI exit strategy
  npi_exit = list(
    npi2 = list(name   = "Further relax April 12 & May 3", 
                phase1 = list(start = "2021-03-22", new_normal_osi = 0.535), 
                phase2 = list(start = "2021-04-12", new_normal_osi = 0.485), 
                phase3 = list(start = "2021-05-03", new_normal_osi = 0.435)))
  
  # Vaccine rollout defined by number_vaccines
  roll_out = list(
    defvax = list(name   = "Proposed vaccine rollout", 
                  group1 = list(coverage = 0.0, start = "2021-01-01", duration = 0)))
  
  # Model parameters to be altered
  alter_params = c(paste0("p$number_vaccines[o$dates_all %in% seq(max(o$dates_data) + 1, format_date('2021-04-01'), by = 'day')] = ", 
                          "round(seq(6000, 25000, length.out = as.numeric(format_date('2021-04-01') - max(o$dates_data))))"), 
                   "p$number_vaccines[o$dates_all > '2021-04-01'] = 25000", 
                   "p$vaccine_effect = 'mixed'", 
                   "p$transmission_blocking_scaler = 0.8",
                   "p$disease_prevention_scaler    = 0.8",
                   "p$variant_infectivity['B117']  = 1.6")
  
  # Format strategy details into single list
  details = list(roll_out = roll_out, npi_exit = npi_exit, alter_params = alter_params)
  
  # Format strategy details into single list
  details = list(roll_out = roll_out, npi_exit = npi_exit, alter_params = alter_params)
  
  return(details)
}

# Scenario 2b) Keep relaxing, fast vaccine
s2b = function() {
  
  # NPI exit strategy
  npi_exit = list(
    npi2 = list(name   = "Further relax April 12 & May 3", 
                phase1 = list(start = "2021-03-22", new_normal_osi = 0.535), 
                phase2 = list(start = "2021-04-12", new_normal_osi = 0.485), 
                phase3 = list(start = "2021-05-03", new_normal_osi = 0.435)))
  
  # Vaccine rollout defined by number_vaccines
  roll_out = list(
    fastvax = list(name   = "Fast vaccine rollout", 
                   group1 = list(coverage = 0.0, start = "2021-01-01", duration = 0)))
  
  # Model parameters to be altered
  alter_params = c(paste0("p$number_vaccines[o$dates_all %in% seq(max(o$dates_data) + 1, format_date('2021-04-01'), by = 'day')] = ", 
                          "round(seq(6000, 50000, length.out = as.numeric(format_date('2021-04-01') - max(o$dates_data))))"), 
                   "p$number_vaccines[o$dates_all > '2021-04-01'] = 50000", 
                   "p$vaccine_effect = 'mixed'", 
                   "p$transmission_blocking_scaler = 0.8",
                   "p$disease_prevention_scaler    = 0.8",
                   "p$variant_infectivity['B117']  = 1.6")
  
  # Format strategy details into single list
  details = list(roll_out = roll_out, npi_exit = npi_exit, alter_params = alter_params)
  
  return(details)
}

# Scenario 3a) Slower relax, default vaccine
s3a = function() {
  
  # NPI exit strategy
  npi_exit = list(
    npi3 = list(name   = "Slower relax until July 5", 
                phase1 = list(start = "2021-03-22", new_normal_osi = 0.559), 
                phase2 = list(start = "2021-04-12", new_normal_osi = 0.535), 
                phase3 = list(start = "2021-05-03", new_normal_osi = 0.510), 
                phase4 = list(start = "2021-05-24", new_normal_osi = 0.485), 
                phase5 = list(start = "2021-06-14", new_normal_osi = 0.460), 
                phase6 = list(start = "2021-07-05", new_normal_osi = 0.435)))
  
  # Vaccine rollout defined by number_vaccines
  roll_out = list(
    defvax = list(name   = "Proposed vaccine rollout", 
                  group1 = list(coverage = 0.0, start = "2021-01-01", duration = 0)))
  
  # Model parameters to be altered
  alter_params = c(paste0("p$number_vaccines[o$dates_all %in% seq(max(o$dates_data) + 1, format_date('2021-04-01'), by = 'day')] = ", 
                          "round(seq(6000, 25000, length.out = as.numeric(format_date('2021-04-01') - max(o$dates_data))))"), 
                   "p$number_vaccines[o$dates_all > '2021-04-01'] = 25000", 
                   "p$vaccine_effect = 'mixed'", 
                   "p$transmission_blocking_scaler = 0.8",
                   "p$disease_prevention_scaler    = 0.8",
                   "p$variant_infectivity['B117']  = 1.6")
  
  # Format strategy details into single list
  details = list(roll_out = roll_out, npi_exit = npi_exit, alter_params = alter_params)
  
  return(details)
}

# Scenario 3a) Slower relax, fast vaccine
s3b = function() {
  
  # NPI exit strategy
  npi_exit = list(
    npi3 = list(name   = "Slower relax until July 5", 
                phase1 = list(start = "2021-03-22", new_normal_osi = 0.559), 
                phase2 = list(start = "2021-04-12", new_normal_osi = 0.535), 
                phase3 = list(start = "2021-05-03", new_normal_osi = 0.510), 
                phase4 = list(start = "2021-05-24", new_normal_osi = 0.485), 
                phase5 = list(start = "2021-06-14", new_normal_osi = 0.460), 
                phase6 = list(start = "2021-07-05", new_normal_osi = 0.435)))
  
  # Vaccine rollout defined by number_vaccines
  roll_out = list(
    fastvax = list(name   = "Fast vaccine rollout", 
                   group1 = list(coverage = 0.0, start = "2021-01-01", duration = 0)))
  
  # Model parameters to be altered
  alter_params = c(paste0("p$number_vaccines[o$dates_all %in% seq(max(o$dates_data) + 1, format_date('2021-04-01'), by = 'day')] = ", 
                          "round(seq(6000, 50000, length.out = as.numeric(format_date('2021-04-01') - max(o$dates_data))))"), 
                   "p$number_vaccines[o$dates_all > '2021-04-01'] = 50000", 
                   "p$vaccine_effect = 'mixed'", 
                   "p$transmission_blocking_scaler = 0.8",
                   "p$disease_prevention_scaler    = 0.8",
                   "p$variant_infectivity['B117']  = 1.6")
  
  # Format strategy details into single list
  details = list(roll_out = roll_out, npi_exit = npi_exit, alter_params = alter_params)
  
  # Format strategy details into single list
  details = list(roll_out = roll_out, npi_exit = npi_exit, alter_params = alter_params)
  
  return(details)
}

# Scenario 4a1) Uncertainty (default vaccine), 50% B117 infectiousness increase
s4a1 = function() {
  
  # NPI exit strategy
  npi_exit = list(
    npi = list(name   = "Slower relax until July 5", 
               phase1 = list(start = "2021-03-22", new_normal_osi = 0.559), 
               phase2 = list(start = "2021-04-12", new_normal_osi = 0.535), 
               phase3 = list(start = "2021-05-03", new_normal_osi = 0.510), 
               phase4 = list(start = "2021-05-24", new_normal_osi = 0.485), 
               phase5 = list(start = "2021-06-14", new_normal_osi = 0.460), 
               phase6 = list(start = "2021-07-05", new_normal_osi = 0.435)))
  
  # Vaccine rollout defined by number_vaccines
  roll_out = list(
    vax = list(name   = "Proposed vaccine rollout", 
               group1 = list(coverage = 0.0, start = "2021-01-01", duration = 0)))
  
  # Model parameters to be altered
  alter_params = c(paste0("p$number_vaccines[o$dates_all %in% seq(max(o$dates_data) + 1, format_date('2021-04-01'), by = 'day')] = ", 
                          "round(seq(6000, 25000, length.out = as.numeric(format_date('2021-04-01') - max(o$dates_data))))"), 
                   "p$number_vaccines[o$dates_all > '2021-04-01'] = 25000", 
                   "p$vaccine_effect = 'mixed'", 
                   "p$transmission_blocking_scaler = 0.8",
                   "p$disease_prevention_scaler    = 0.8",
                   "p$variant_infectivity['B117']  = 1.5")
  
  # Format strategy details into single list
  details = list(roll_out = roll_out, npi_exit = npi_exit, alter_params = alter_params)
  
  return(details)
}

# Scenario 4a2) Uncertainty (default vaccine), 70% B117 infectiousness increase
s4a2 = function() {
  
  # NPI exit strategy
  npi_exit = list(
    npi = list(name   = "Slower relax until July 5", 
               phase1 = list(start = "2021-03-22", new_normal_osi = 0.559), 
               phase2 = list(start = "2021-04-12", new_normal_osi = 0.535), 
               phase3 = list(start = "2021-05-03", new_normal_osi = 0.510), 
               phase4 = list(start = "2021-05-24", new_normal_osi = 0.485), 
               phase5 = list(start = "2021-06-14", new_normal_osi = 0.460), 
               phase6 = list(start = "2021-07-05", new_normal_osi = 0.435)))
  
  # Vaccine rollout defined by number_vaccines
  roll_out = list(
    vax = list(name   = "Proposed vaccine rollout", 
               group1 = list(coverage = 0.0, start = "2021-01-01", duration = 0)))
  
  # Model parameters to be altered
  alter_params = c(paste0("p$number_vaccines[o$dates_all %in% seq(max(o$dates_data) + 1, format_date('2021-04-01'), by = 'day')] = ", 
                          "round(seq(6000, 25000, length.out = as.numeric(format_date('2021-04-01') - max(o$dates_data))))"), 
                   "p$number_vaccines[o$dates_all > '2021-04-01'] = 25000", 
                   "p$vaccine_effect = 'mixed'", 
                   "p$transmission_blocking_scaler = 0.8",
                   "p$disease_prevention_scaler    = 0.8",
                   "p$variant_infectivity['B117']  = 1.7")
  
  # Format strategy details into single list
  details = list(roll_out = roll_out, npi_exit = npi_exit, alter_params = alter_params)
  
  return(details)
}

# Scenario 4b1) Uncertainty (fast vaccine), 50% B117 infectiousness increase
s4b1 = function() {
  
  # NPI exit strategy
  npi_exit = list(
    npi = list(name   = "Slower relax until July 5", 
                phase1 = list(start = "2021-03-22", new_normal_osi = 0.559), 
                phase2 = list(start = "2021-04-12", new_normal_osi = 0.535), 
                phase3 = list(start = "2021-05-03", new_normal_osi = 0.510), 
                phase4 = list(start = "2021-05-24", new_normal_osi = 0.485), 
                phase5 = list(start = "2021-06-14", new_normal_osi = 0.460), 
                phase6 = list(start = "2021-07-05", new_normal_osi = 0.435)))
  
  # Vaccine rollout defined by number_vaccines
  roll_out = list(
    vax = list(name   = "Fast vaccine rollout", 
                       group1 = list(coverage = 0.0, start = "2021-01-01", duration = 0)))
  
  # Model parameters to be altered
  alter_params = c(paste0("p$number_vaccines[o$dates_all %in% seq(max(o$dates_data) + 1, format_date('2021-04-01'), by = 'day')] = ", 
                          "round(seq(6000, 50000, length.out = as.numeric(format_date('2021-04-01') - max(o$dates_data))))"), 
                   "p$number_vaccines[o$dates_all > '2021-04-01'] = 50000", 
                   "p$vaccine_effect = 'mixed'", 
                   "p$transmission_blocking_scaler = 0.8",
                   "p$disease_prevention_scaler    = 0.8",
                   "p$variant_infectivity['B117']  = 1.5")
  
  # Format strategy details into single list
  details = list(roll_out = roll_out, npi_exit = npi_exit, alter_params = alter_params)
  
  # Format strategy details into single list
  details = list(roll_out = roll_out, npi_exit = npi_exit, alter_params = alter_params)
  
  return(details)
}

# Scenario 4b2) Uncertainty (fast vaccine), 70% B117 infectiousness increase
s4b2 = function() {
  
  # NPI exit strategy
  npi_exit = list(
    npi = list(name   = "Slower relax until July 5", 
                phase1 = list(start = "2021-03-22", new_normal_osi = 0.559), 
                phase2 = list(start = "2021-04-12", new_normal_osi = 0.535), 
                phase3 = list(start = "2021-05-03", new_normal_osi = 0.510), 
                phase4 = list(start = "2021-05-24", new_normal_osi = 0.485), 
                phase5 = list(start = "2021-06-14", new_normal_osi = 0.460), 
                phase6 = list(start = "2021-07-05", new_normal_osi = 0.435)))
  
  # Vaccine rollout defined by number_vaccines
  roll_out = list(
    vax = list(name   = "Fast vaccine rollout", 
                       group1 = list(coverage = 0.0, start = "2021-01-01", duration = 0)))
  
  # Model parameters to be altered
  alter_params = c(paste0("p$number_vaccines[o$dates_all %in% seq(max(o$dates_data) + 1, format_date('2021-04-01'), by = 'day')] = ", 
                          "round(seq(6000, 50000, length.out = as.numeric(format_date('2021-04-01') - max(o$dates_data))))"), 
                   "p$number_vaccines[o$dates_all > '2021-04-01'] = 50000", 
                   "p$vaccine_effect = 'mixed'", 
                   "p$transmission_blocking_scaler = 0.8",
                   "p$disease_prevention_scaler    = 0.8",
                   "p$variant_infectivity['B117']  = 1.7")
  
  # Format strategy details into single list
  details = list(roll_out = roll_out, npi_exit = npi_exit, alter_params = alter_params)
  
  # Format strategy details into single list
  details = list(roll_out = roll_out, npi_exit = npi_exit, alter_params = alter_params)
  
  return(details)
}

# Scenario 5a1) Uncertainty (default vaccine), 60% transmission blocking
s5a1 = function() {
  
  # NPI exit strategy
  npi_exit = list(
    npi = list(name   = "Slower relax until July 5", 
               phase1 = list(start = "2021-03-22", new_normal_osi = 0.559), 
               phase2 = list(start = "2021-04-12", new_normal_osi = 0.535), 
               phase3 = list(start = "2021-05-03", new_normal_osi = 0.510), 
               phase4 = list(start = "2021-05-24", new_normal_osi = 0.485), 
               phase5 = list(start = "2021-06-14", new_normal_osi = 0.460), 
               phase6 = list(start = "2021-07-05", new_normal_osi = 0.435)))
  
  # Vaccine rollout defined by number_vaccines
  roll_out = list(
    vax = list(name   = "Proposed vaccine rollout", 
               group1 = list(coverage = 0.0, start = "2021-01-01", duration = 0)))
  
  # Model parameters to be altered
  alter_params = c(paste0("p$number_vaccines[o$dates_all %in% seq(max(o$dates_data) + 1, format_date('2021-04-01'), by = 'day')] = ", 
                          "round(seq(6000, 25000, length.out = as.numeric(format_date('2021-04-01') - max(o$dates_data))))"), 
                   "p$number_vaccines[o$dates_all > '2021-04-01'] = 25000", 
                   "p$vaccine_effect = 'mixed'", 
                   "p$transmission_blocking_scaler = 0.6",
                   "p$disease_prevention_scaler    = 0.88",
                   "p$variant_infectivity['B117']  = 1.6")
  
  # Format strategy details into single list
  details = list(roll_out = roll_out, npi_exit = npi_exit, alter_params = alter_params)
  
  return(details)
}

# Scenario 5a2) Uncertainty (default vaccine), 98% transmission blocking
s5a2 = function() {
  
  # NPI exit strategy
  npi_exit = list(
    npi3 = list(name   = "Slower relax until July 5", 
                phase1 = list(start = "2021-03-22", new_normal_osi = 0.559), 
                phase2 = list(start = "2021-04-12", new_normal_osi = 0.535), 
                phase3 = list(start = "2021-05-03", new_normal_osi = 0.510), 
                phase4 = list(start = "2021-05-24", new_normal_osi = 0.485), 
                phase5 = list(start = "2021-06-14", new_normal_osi = 0.460), 
                phase6 = list(start = "2021-07-05", new_normal_osi = 0.435)))
  
  # Vaccine rollout defined by number_vaccines
  roll_out = list(
    defvax_t98 = list(name   = "Proposed vaccine rollout", 
                      group1 = list(coverage = 0.0, start = "2021-01-01", duration = 0)))
  
  # Model parameters to be altered
  alter_params = c(paste0("p$number_vaccines[o$dates_all %in% seq(max(o$dates_data) + 1, format_date('2021-04-01'), by = 'day')] = ", 
                          "round(seq(6000, 25000, length.out = as.numeric(format_date('2021-04-01') - max(o$dates_data))))"), 
                   "p$number_vaccines[o$dates_all > '2021-04-01'] = 25000", 
                   "p$vaccine_effect = 'mixed'", 
                   "p$transmission_blocking_scaler = 0.98",
                   "p$disease_prevention_scaler    = 0.28",
                   "p$variant_infectivity['B117']  = 1.6")
  
  # Format strategy details into single list
  details = list(roll_out = roll_out, npi_exit = npi_exit, alter_params = alter_params)
  
  return(details)
}

# Scenario 5b1) Uncertainty (fast vaccine), 60% transmission blocking
s5b1 = function() {
  
  # NPI exit strategy
  npi_exit = list(
    npi3 = list(name   = "Slower relax until July 5", 
                phase1 = list(start = "2021-03-22", new_normal_osi = 0.559), 
                phase2 = list(start = "2021-04-12", new_normal_osi = 0.535), 
                phase3 = list(start = "2021-05-03", new_normal_osi = 0.510), 
                phase4 = list(start = "2021-05-24", new_normal_osi = 0.485), 
                phase5 = list(start = "2021-06-14", new_normal_osi = 0.460), 
                phase6 = list(start = "2021-07-05", new_normal_osi = 0.435)))
  
  # Vaccine rollout defined by number_vaccines
  roll_out = list(
    fastvax_t60 = list(name   = "Fast vaccine rollout", 
                       group1 = list(coverage = 0.0, start = "2021-01-01", duration = 0)))
  
  # Model parameters to be altered
  alter_params = c(paste0("p$number_vaccines[o$dates_all %in% seq(max(o$dates_data) + 1, format_date('2021-04-01'), by = 'day')] = ", 
                          "round(seq(6000, 50000, length.out = as.numeric(format_date('2021-04-01') - max(o$dates_data))))"), 
                   "p$number_vaccines[o$dates_all > '2021-04-01'] = 50000", 
                   "p$vaccine_effect = 'mixed'", 
                   "p$transmission_blocking_scaler = 0.6",
                   "p$disease_prevention_scaler    = 0.88",
                   "p$variant_infectivity['B117']  = 1.6")
  
  # Format strategy details into single list
  details = list(roll_out = roll_out, npi_exit = npi_exit, alter_params = alter_params)
  
  # Format strategy details into single list
  details = list(roll_out = roll_out, npi_exit = npi_exit, alter_params = alter_params)
  
  return(details)
}

# Scenario 5b2) Uncertainty (fast vaccine), 98% transmission blocking
s5b2 = function() {
  
  # NPI exit strategy
  npi_exit = list(
    npi3 = list(name   = "Slower relax until July 5", 
                phase1 = list(start = "2021-03-22", new_normal_osi = 0.559), 
                phase2 = list(start = "2021-04-12", new_normal_osi = 0.535), 
                phase3 = list(start = "2021-05-03", new_normal_osi = 0.510), 
                phase4 = list(start = "2021-05-24", new_normal_osi = 0.485), 
                phase5 = list(start = "2021-06-14", new_normal_osi = 0.460), 
                phase6 = list(start = "2021-07-05", new_normal_osi = 0.435)))
  
  # Vaccine rollout defined by number_vaccines
  roll_out = list(
    fastvax_t98 = list(name   = "Fast vaccine rollout", 
                       group1 = list(coverage = 0.0, start = "2021-01-01", duration = 0)))
  
  # Model parameters to be altered
  alter_params = c(paste0("p$number_vaccines[o$dates_all %in% seq(max(o$dates_data) + 1, format_date('2021-04-01'), by = 'day')] = ", 
                          "round(seq(6000, 50000, length.out = as.numeric(format_date('2021-04-01') - max(o$dates_data))))"), 
                   "p$number_vaccines[o$dates_all > '2021-04-01'] = 50000", 
                   "p$vaccine_effect = 'mixed'", 
                   "p$transmission_blocking_scaler = 0.98",
                   "p$disease_prevention_scaler    = 0.28",
                   "p$variant_infectivity['B117']  = 1.6")
  
  # Format strategy details into single list
  details = list(roll_out = roll_out, npi_exit = npi_exit, alter_params = alter_params)
  
  # Format strategy details into single list
  details = list(roll_out = roll_out, npi_exit = npi_exit, alter_params = alter_params)
  
  return(details)
}

# Scenario 6a1) Uncertainty (default vaccine), low vaccine acceptance
s6a1 = function() {
  
  # NPI exit strategy
  npi_exit = list(
    npi3 = list(name   = "Slower relax until July 5", 
                phase1 = list(start = "2021-03-22", new_normal_osi = 0.559), 
                phase2 = list(start = "2021-04-12", new_normal_osi = 0.535), 
                phase3 = list(start = "2021-05-03", new_normal_osi = 0.510), 
                phase4 = list(start = "2021-05-24", new_normal_osi = 0.485), 
                phase5 = list(start = "2021-06-14", new_normal_osi = 0.460), 
                phase6 = list(start = "2021-07-05", new_normal_osi = 0.435)))
  
  # Vaccine rollout defined by number_vaccines
  roll_out = list(
    defvax_low = list(name   = "Proposed vaccine rollout", 
                      group1 = list(coverage = 0.0, start = "2021-01-01", duration = 0)))
  
  # Model parameters to be altered
  alter_params = c(paste0("p$number_vaccines[o$dates_all %in% seq(max(o$dates_data) + 1, format_date('2021-04-01'), by = 'day')] = ", 
                          "round(seq(6000, 25000, length.out = as.numeric(format_date('2021-04-01') - max(o$dates_data))))"), 
                   "p$number_vaccines[o$dates_all > '2021-04-01'] = 25000", 
                   "p$vaccine_effect = 'mixed'", 
                   "p$transmission_blocking_scaler = 0.8",
                   "p$disease_prevention_scaler    = 0.8",
                   "p$variant_infectivity['B117']  = 1.6", 
                   "p$vaccine_acceptance = c(group1=0.75, group2=0.6, group3=0.6, group4=0.6, group5=0.6)")
  
  # Format strategy details into single list
  details = list(roll_out = roll_out, npi_exit = npi_exit, alter_params = alter_params)
  
  return(details)
}

# Scenario 6a2) Uncertainty (default vaccine), high vaccine acceptance
s6a2 = function() {
  
  # NPI exit strategy
  npi_exit = list(
    npi3 = list(name   = "Slower relax until July 5", 
                phase1 = list(start = "2021-03-22", new_normal_osi = 0.559), 
                phase2 = list(start = "2021-04-12", new_normal_osi = 0.535), 
                phase3 = list(start = "2021-05-03", new_normal_osi = 0.510), 
                phase4 = list(start = "2021-05-24", new_normal_osi = 0.485), 
                phase5 = list(start = "2021-06-14", new_normal_osi = 0.460), 
                phase6 = list(start = "2021-07-05", new_normal_osi = 0.435)))
  
  # Vaccine rollout defined by number_vaccines
  roll_out = list(
    defvax_high = list(name   = "Proposed vaccine rollout", 
                       group1 = list(coverage = 0.0, start = "2021-01-01", duration = 0)))
  
  # Model parameters to be altered
  alter_params = c(paste0("p$number_vaccines[o$dates_all %in% seq(max(o$dates_data) + 1, format_date('2021-04-01'), by = 'day')] = ", 
                          "round(seq(6000, 25000, length.out = as.numeric(format_date('2021-04-01') - max(o$dates_data))))"), 
                   "p$number_vaccines[o$dates_all > '2021-04-01'] = 25000", 
                   "p$vaccine_effect = 'mixed'", 
                   "p$transmission_blocking_scaler = 0.8",
                   "p$disease_prevention_scaler    = 0.8",
                   "p$variant_infectivity['B117']  = 1.6", 
                   "p$vaccine_acceptance = c(group1=0.75, group2=0.9, group3=0.9, group4=0.9, group5=0.9)")
  
  # Format strategy details into single list
  details = list(roll_out = roll_out, npi_exit = npi_exit, alter_params = alter_params)
  
  return(details)
}

# Scenario 6b1) Uncertainty (fast vaccine), low vaccine acceptance
s6b1 = function() {
  
  # NPI exit strategy
  npi_exit = list(
    npi3 = list(name   = "Slower relax until July 5", 
                phase1 = list(start = "2021-03-22", new_normal_osi = 0.559), 
                phase2 = list(start = "2021-04-12", new_normal_osi = 0.535), 
                phase3 = list(start = "2021-05-03", new_normal_osi = 0.510), 
                phase4 = list(start = "2021-05-24", new_normal_osi = 0.485), 
                phase5 = list(start = "2021-06-14", new_normal_osi = 0.460), 
                phase6 = list(start = "2021-07-05", new_normal_osi = 0.435)))
  
  # Vaccine rollout defined by number_vaccines
  roll_out = list(
    fastvax_low = list(name   = "Fast vaccine rollout", 
                       group1 = list(coverage = 0.0, start = "2021-01-01", duration = 0)))
  
  # Model parameters to be altered
  alter_params = c(paste0("p$number_vaccines[o$dates_all %in% seq(max(o$dates_data) + 1, format_date('2021-04-01'), by = 'day')] = ", 
                          "round(seq(6000, 50000, length.out = as.numeric(format_date('2021-04-01') - max(o$dates_data))))"), 
                   "p$number_vaccines[o$dates_all > '2021-04-01'] = 50000", 
                   "p$vaccine_effect = 'mixed'", 
                   "p$transmission_blocking_scaler = 0.8",
                   "p$disease_prevention_scaler    = 0.8",
                   "p$variant_infectivity['B117']  = 1.6",
                   "p$vaccine_acceptance = c(group1=0.75, group2=0.6, group3=0.6, group4=0.6, group5=0.6)")
  
  # Format strategy details into single list
  details = list(roll_out = roll_out, npi_exit = npi_exit, alter_params = alter_params)
  
  # Format strategy details into single list
  details = list(roll_out = roll_out, npi_exit = npi_exit, alter_params = alter_params)
  
  return(details)
}

# Scenario 6b2) Uncertainty (fast vaccine), high vaccine acceptance
s6b2 = function() {
  
  # NPI exit strategy
  npi_exit = list(
    npi3 = list(name   = "Slower relax until July 5", 
                phase1 = list(start = "2021-03-22", new_normal_osi = 0.559), 
                phase2 = list(start = "2021-04-12", new_normal_osi = 0.535), 
                phase3 = list(start = "2021-05-03", new_normal_osi = 0.510), 
                phase4 = list(start = "2021-05-24", new_normal_osi = 0.485), 
                phase5 = list(start = "2021-06-14", new_normal_osi = 0.460), 
                phase6 = list(start = "2021-07-05", new_normal_osi = 0.435)))
  
  # Vaccine rollout defined by number_vaccines
  roll_out = list(
    fastvax_high = list(name   = "Fast vaccine rollout", 
                        group1 = list(coverage = 0.0, start = "2021-01-01", duration = 0)))
  
  # Model parameters to be altered
  alter_params = c(paste0("p$number_vaccines[o$dates_all %in% seq(max(o$dates_data) + 1, format_date('2021-04-01'), by = 'day')] = ", 
                          "round(seq(6000, 50000, length.out = as.numeric(format_date('2021-04-01') - max(o$dates_data))))"), 
                   "p$number_vaccines[o$dates_all > '2021-04-01'] = 50000", 
                   "p$vaccine_effect = 'mixed'", 
                   "p$transmission_blocking_scaler = 0.8",
                   "p$disease_prevention_scaler    = 0.8",
                   "p$variant_infectivity['B117']  = 1.6",
                   "p$vaccine_acceptance = c(group1=0.75, group2=0.9, group3=0.9, group4=0.9, group5=0.9)")
  
  # Format strategy details into single list
  details = list(roll_out = roll_out, npi_exit = npi_exit, alter_params = alter_params)
  
  # Format strategy details into single list
  details = list(roll_out = roll_out, npi_exit = npi_exit, alter_params = alter_params)
  
  return(details)
}

# Scenario 7a) Uncertainty (default vaccine), 50% B117 mortality increase
s7a = function() {
  
  # NPI exit strategy
  npi_exit = list(
    npi3 = list(name   = "Slower relax until July 5", 
                phase1 = list(start = "2021-03-22", new_normal_osi = 0.559), 
                phase2 = list(start = "2021-04-12", new_normal_osi = 0.535), 
                phase3 = list(start = "2021-05-03", new_normal_osi = 0.510), 
                phase4 = list(start = "2021-05-24", new_normal_osi = 0.485), 
                phase5 = list(start = "2021-06-14", new_normal_osi = 0.460), 
                phase6 = list(start = "2021-07-05", new_normal_osi = 0.435)))
  
  # Vaccine rollout defined by number_vaccines
  roll_out = list(
    defvax_high = list(name   = "Proposed vaccine rollout", 
                       group1 = list(coverage = 0.0, start = "2021-01-01", duration = 0)))
  
  # Model parameters to be altered
  alter_params = c(paste0("p$number_vaccines[o$dates_all %in% seq(max(o$dates_data) + 1, format_date('2021-04-01'), by = 'day')] = ", 
                          "round(seq(6000, 25000, length.out = as.numeric(format_date('2021-04-01') - max(o$dates_data))))"), 
                   "p$number_vaccines[o$dates_all > '2021-04-01'] = 25000", 
                   "p$vaccine_effect = 'mixed'", 
                   "p$transmission_blocking_scaler = 0.8",
                   "p$disease_prevention_scaler    = 0.8",
                   "p$variant_infectivity['B117']  = 1.6", 
                   "p$variant_severity['B117']     = 1.5")
  
  # Format strategy details into single list
  details = list(roll_out = roll_out, npi_exit = npi_exit, alter_params = alter_params)
  
  return(details)
}

# Scenario 7b) Uncertainty (fast vaccine), 50% B117 mortality increase
s7b = function() {
  
  # NPI exit strategy
  npi_exit = list(
    npi3 = list(name   = "Slower relax until July 5", 
                phase1 = list(start = "2021-03-22", new_normal_osi = 0.559), 
                phase2 = list(start = "2021-04-12", new_normal_osi = 0.535), 
                phase3 = list(start = "2021-05-03", new_normal_osi = 0.510), 
                phase4 = list(start = "2021-05-24", new_normal_osi = 0.485), 
                phase5 = list(start = "2021-06-14", new_normal_osi = 0.460), 
                phase6 = list(start = "2021-07-05", new_normal_osi = 0.435)))
  
  # Vaccine rollout defined by number_vaccines
  roll_out = list(
    fastvax_low = list(name   = "Fast vaccine rollout", 
                       group1 = list(coverage = 0.0, start = "2021-01-01", duration = 0)))
  
  # Model parameters to be altered
  alter_params = c(paste0("p$number_vaccines[o$dates_all %in% seq(max(o$dates_data) + 1, format_date('2021-04-01'), by = 'day')] = ", 
                          "round(seq(6000, 50000, length.out = as.numeric(format_date('2021-04-01') - max(o$dates_data))))"), 
                   "p$number_vaccines[o$dates_all > '2021-04-01'] = 50000", 
                   "p$vaccine_effect = 'mixed'", 
                   "p$transmission_blocking_scaler = 0.8",
                   "p$disease_prevention_scaler    = 0.8",
                   "p$variant_infectivity['B117']  = 1.6", 
                   "p$variant_severity['B117']     = 1.5")
  
  # Format strategy details into single list
  details = list(roll_out = roll_out, npi_exit = npi_exit, alter_params = alter_params)
  
  # Format strategy details into single list
  details = list(roll_out = roll_out, npi_exit = npi_exit, alter_params = alter_params)
  
  return(details)
}

# Scenario 8a) March 22nd only (3w delay), default vaccine
s8a = function() {
  
  # NPI exit strategy
  npi_exit = list(
    npi8 = list(name   = "Relax April 12", 
                phase1 = list(start = "2021-04-12", new_normal_osi = 0.535)))
  
  # Vaccine rollout defined by number_vaccines
  roll_out = list(
    defvax = list(name   = "Proposed vaccine rollout", 
                  group1 = list(coverage = 0.0, start = "2021-01-01", duration = 0)))
  
  # Model parameters to be altered
  alter_params = c(paste0("p$number_vaccines[o$dates_all %in% seq(max(o$dates_data) + 1, format_date('2021-04-01'), by = 'day')] = ", 
                          "round(seq(6000, 25000, length.out = as.numeric(format_date('2021-04-01') - max(o$dates_data))))"), 
                   "p$number_vaccines[o$dates_all > '2021-04-01'] = 25000", 
                   "p$vaccine_effect = 'mixed'", 
                   "p$transmission_blocking_scaler = 0.8",
                   "p$disease_prevention_scaler    = 0.8",
                   "p$variant_infectivity['B117']  = 1.6")
  
  # Format strategy details into single list
  details = list(roll_out = roll_out, npi_exit = npi_exit, alter_params = alter_params)
  
  return(details)
}

# Scenario 8b) March 22nd only (3w delay), fast vaccine
s8b = function() {
  
  # NPI exit strategy
  npi_exit = list(
    npi8 = list(name   = "Relax April 12", 
                phase1 = list(start = "2021-04-12", new_normal_osi = 0.535)))
  
  # Vaccine rollout defined by number_vaccines
  roll_out = list(
    fastvax = list(name   = "Fast vaccine rollout", 
                   group1 = list(coverage = 0.0, start = "2021-01-01", duration = 0)))
  
  # Model parameters to be altered
  alter_params = c(paste0("p$number_vaccines[o$dates_all %in% seq(max(o$dates_data) + 1, format_date('2021-04-01'), by = 'day')] = ", 
                          "round(seq(6000, 50000, length.out = as.numeric(format_date('2021-04-01') - max(o$dates_data))))"), 
                   "p$number_vaccines[o$dates_all > '2021-04-01'] = 50000", 
                   "p$vaccine_effect = 'mixed'", 
                   "p$transmission_blocking_scaler = 0.8",
                   "p$disease_prevention_scaler    = 0.8",
                   "p$variant_infectivity['B117']  = 1.6")
  
  # Format strategy details into single list
  details = list(roll_out = roll_out, npi_exit = npi_exit, alter_params = alter_params)
  
  return(details)
}

# Scenario 9a) March 22nd only (6w delay), default vaccine
s9a = function() {
  
  # NPI exit strategy
  npi_exit = list(
    npi9 = list(name   = "Relax May 3", 
                phase1 = list(start = "2021-05-03", new_normal_osi = 0.535)))
  
  # Vaccine rollout defined by number_vaccines
  roll_out = list(
    defvax = list(name   = "Proposed vaccine rollout", 
                  group1 = list(coverage = 0.0, start = "2021-01-01", duration = 0)))
  
  # Model parameters to be altered
  alter_params = c(paste0("p$number_vaccines[o$dates_all %in% seq(max(o$dates_data) + 1, format_date('2021-04-01'), by = 'day')] = ", 
                          "round(seq(6000, 25000, length.out = as.numeric(format_date('2021-04-01') - max(o$dates_data))))"), 
                   "p$number_vaccines[o$dates_all > '2021-04-01'] = 25000", 
                   "p$vaccine_effect = 'mixed'", 
                   "p$transmission_blocking_scaler = 0.8",
                   "p$disease_prevention_scaler    = 0.8",
                   "p$variant_infectivity['B117']  = 1.6")
  
  # Format strategy details into single list
  details = list(roll_out = roll_out, npi_exit = npi_exit, alter_params = alter_params)
  
  return(details)
}

# Scenario 9b) March 22nd only (6w delay), fast vaccine
s9b = function() {
  
  # NPI exit strategy
  npi_exit = list(
    npi9 = list(name   = "Relax May 3", 
                phase1 = list(start = "2021-05-03", new_normal_osi = 0.535)))
  
  # Vaccine rollout defined by number_vaccines
  roll_out = list(
    fastvax = list(name   = "Fast vaccine rollout", 
                   group1 = list(coverage = 0.0, start = "2021-01-01", duration = 0)))
  
  # Model parameters to be altered
  alter_params = c(paste0("p$number_vaccines[o$dates_all %in% seq(max(o$dates_data) + 1, format_date('2021-04-01'), by = 'day')] = ", 
                          "round(seq(6000, 50000, length.out = as.numeric(format_date('2021-04-01') - max(o$dates_data))))"), 
                   "p$number_vaccines[o$dates_all > '2021-04-01'] = 50000", 
                   "p$vaccine_effect = 'mixed'", 
                   "p$transmission_blocking_scaler = 0.8",
                   "p$disease_prevention_scaler    = 0.8",
                   "p$variant_infectivity['B117']  = 1.6")
  
  # Format strategy details into single list
  details = list(roll_out = roll_out, npi_exit = npi_exit, alter_params = alter_params)
  
  return(details)
}

# Scenario 10a) March 22nd only (9w delay), default vaccine
s10a = function() {
  
  # NPI exit strategy
  npi_exit = list(
    npi10 = list(name   = "Relax May 24", 
                phase1 = list(start = "2021-05-24", new_normal_osi = 0.535)))
  
  # Vaccine rollout defined by number_vaccines
  roll_out = list(
    defvax = list(name   = "Proposed vaccine rollout", 
                  group1 = list(coverage = 0.0, start = "2021-01-01", duration = 0)))
  
  # Model parameters to be altered
  alter_params = c(paste0("p$number_vaccines[o$dates_all %in% seq(max(o$dates_data) + 1, format_date('2021-04-01'), by = 'day')] = ", 
                          "round(seq(6000, 25000, length.out = as.numeric(format_date('2021-04-01') - max(o$dates_data))))"), 
                   "p$number_vaccines[o$dates_all > '2021-04-01'] = 25000", 
                   "p$vaccine_effect = 'mixed'", 
                   "p$transmission_blocking_scaler = 0.8",
                   "p$disease_prevention_scaler    = 0.8",
                   "p$variant_infectivity['B117']  = 1.6")
  
  # Format strategy details into single list
  details = list(roll_out = roll_out, npi_exit = npi_exit, alter_params = alter_params)
  
  return(details)
}

# Scenario 10b) March 22nd only (9w delay), fast vaccine
s10b = function() {
  
  # NPI exit strategy
  npi_exit = list(
    npi10 = list(name   = "Relax May 24", 
                phase1 = list(start = "2021-05-24", new_normal_osi = 0.535)))
  
  # Vaccine rollout defined by number_vaccines
  roll_out = list(
    fastvax = list(name   = "Fast vaccine rollout", 
                   group1 = list(coverage = 0.0, start = "2021-01-01", duration = 0)))
  
  # Model parameters to be altered
  alter_params = c(paste0("p$number_vaccines[o$dates_all %in% seq(max(o$dates_data) + 1, format_date('2021-04-01'), by = 'day')] = ", 
                          "round(seq(6000, 50000, length.out = as.numeric(format_date('2021-04-01') - max(o$dates_data))))"), 
                   "p$number_vaccines[o$dates_all > '2021-04-01'] = 50000", 
                   "p$vaccine_effect = 'mixed'", 
                   "p$transmission_blocking_scaler = 0.8",
                   "p$disease_prevention_scaler    = 0.8",
                   "p$variant_infectivity['B117']  = 1.6")
  
  # Format strategy details into single list
  details = list(roll_out = roll_out, npi_exit = npi_exit, alter_params = alter_params)
  
  return(details)
}

# ---------------------------------------------------------
# Define commonly used policies for easy referencing
# ---------------------------------------------------------
policies = function(...) {
  
  # Initiate list of policies
  policy = this_policy = list()
  
  # Which policies we want to export
  which_policies = as.character(list(...))
  
  # ---- Vaccine roll out ----
  
  # No vaccination campaign
  policy$novax = list(
    novax = list(name   = "No vaccine rollout", 
                 group1 = list(coverage = 0.0, start = "2021-01-01", duration = 0)))
  
  # Default (aka expected) vaccination campaign
  policy$defvax = list(
    defvax = list(name   = "Expected vaccine rollout", 
                  group1 = list(coverage = 0.75, start = "2021-03-01", duration = 30, vaccine = "pfizer"),
                  group2 = list(coverage = 0.75, start = "2021-04-01", duration = 90, vaccine = "pfizer"),
                  group3 = list(coverage = 0.75, start = "2021-05-01", duration = 61, vaccine = "pfizer"),
                  group4 = list(coverage = 0.75, start = "2021-05-15", duration = 30, vaccine = "pfizer"),
                  group5 = list(coverage = 0.75, start = "2021-06-15", duration = 60, vaccine = "pfizer")))
  
  # Faster vaccination roll out campaign (up to 75k)
  policy$fastvax = list(
    fastvax = list(name   = "Faster vaccine rollout", 
                   group1 = list(coverage = 0.75, start = "2021-01-01", duration = 115, vaccine = "pfizer"),
                   group2 = list(coverage = 0.75, start = "2021-03-01", duration = 46,  vaccine = "pfizer"),
                   group3 = list(coverage = 0.75, start = "2021-03-20", duration = 29,  vaccine = "pfizer"),
                   group4 = list(coverage = 0.75, start = "2021-04-19", duration = 2,   vaccine = "pfizer"),
                   group5 = list(coverage = 0.75, start = "2021-04-22", duration = 44,  vaccine = "pfizer")))
  
  # ---- NPI exit strategies ----
  
  # Phased relaxation of NPIs as now mandated by federal council (TF report 5th March 2021)
  policy$fcplan = list(
    npiplanned  = list(name   = "NPI proposal of Federal Council", 
                       phase1 = list(start = "2021-03-01", delay = 0, new_normal_rel = 0.9126),
                       phase2 = list(start = "2021-03-22", delay = 0, new_normal_rel = 0.9172)),
    npi6wlater  = list(name   = "Same measures 6 week delay", 
                       phase1 = list(start = "2021-03-01", delay = 0, new_normal_rel = 0.9126), 
                       phase2 = list(start = "2021-03-22", delay = 0, new_normal_rel = 0.9724), 
                       phase3 = list(start = "2021-04-12", delay = 0, new_normal_rel = 0.9716), 
                       phase4 = list(start = "2021-05-03", delay = 0, new_normal_rel = 0.9708)), 
    npi12wlater = list(name   = "Same measures 12 week delay", 
                       phase1 = list(start = "2021-03-01", delay = 0, new_normal_rel = 0.9126), 
                       phase2 = list(start = "2021-03-22", delay = 0, new_normal_rel = 0.9834), 
                       phase3 = list(start = "2021-04-12", delay = 0, new_normal_rel = 0.9832), 
                       phase4 = list(start = "2021-05-03", delay = 0, new_normal_rel = 0.9829), 
                       phase5 = list(start = "2021-05-24", delay = 0, new_normal_rel = 0.9826), 
                       phase6 = list(start = "2021-06-14", delay = 0, new_normal_rel = 0.9823)), 
    npip1only   = list(name   = "No further relax after Mar 1", 
                       phase1 = list(start = "2021-03-01", delay = 0, new_normal_rel = 0.9126)))
  
  # Continued phased relaxation of NPIs
  policy$fcplanplus = list(
    npiplanned = list(name   = "Continued relax of measures", 
                      phase1 = list(start = "2021-03-01", delay = 0, new_normal_rel = 0.9126),
                      phase2 = list(start = "2021-03-22", delay = 0, new_normal_rel = 0.9172), 
                      phase3 = list(start = "2021-04-12", delay = 0, new_normal_rel = 0.915), 
                      phase4 = list(start = "2021-04-24", delay = 0, new_normal_rel = 0.915)),
    npi6wlater = list(name   = "Same measures 6 week delay", 
                      phase1 = list(start = "2021-03-01", delay = 0, new_normal_rel = 0.9126), 
                      phase2 = list(start = "2021-03-22", delay = 0, new_normal_rel = 0.9536), 
                      phase3 = list(start = "2021-04-12", delay = 0, new_normal_rel = 0.9513), 
                      phase4 = list(start = "2021-05-03", delay = 0, new_normal_rel = 0.9488), 
                      phase5 = list(start = "2021-05-24", delay = 0, new_normal_rel = 0.9461), 
                      phase6 = list(start = "2021-06-14", delay = 0, new_normal_rel = 0.943)), 
    npiless    = list(name   = "More measures with slower relax", 
                      phase1 = list(start = "2021-03-01", delay = 0, new_normal_rel = 0.9126), 
                      phase2 = list(start = "2021-03-22", delay = 0, new_normal_rel = 0.9834), 
                      phase3 = list(start = "2021-04-12", delay = 0, new_normal_rel = 0.9832), 
                      phase4 = list(start = "2021-05-03", delay = 0, new_normal_rel = 0.9829), 
                      phase5 = list(start = "2021-05-24", delay = 0, new_normal_rel = 0.9826), 
                      phase6 = list(start = "2021-06-14", delay = 0, new_normal_rel = 0.9823)))
  
  # Return policy of interest
  for (which_policy in which_policies)
    this_policy = c(this_policy, policy[[which_policy]])
  
  return(this_policy)
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
    #
    # NOTE: Commenting out as we want we want to allow this behaviour for now
    # if (start_date <= max(o$dates_data))
    #   stop("You are attempting to start an NPI release phase in the past")
    
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

