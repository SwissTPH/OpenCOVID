###########################################################
# PROGNOSIS
#
# Calculate base probilities of each prognosis based on age.
# Also calculate severe disease multiplicate factors based 
# on variant infected with, number of past exposures, and 
# number of vaccine doses.
#
# These factors are used within model.R to generate a prognosis
# for each newly infected individual, at which stage waning 
# immunity (for vaccine-induced and acquired immunity) as well 
# as factors such as comorbidities and PrEP are considered. 
#
###########################################################

# ---------------------------------------------------------
# Prognosis probabilities by age, variant, doses, and exposures
# ---------------------------------------------------------
prognosis_probabilities = function(y) {
  
  # ---- Age-related disease probabilities ----
  
  # Parse function: age-related probabiliy of severe disease given person is symptomatic
  severe_symptom = parse_fn(y$severe_symptom_age, along = list(x = y$ages)) %>% pmax(1e-12)
  
  # Parse function: gge-related probabiliy of critical disease given person is severe
  critical_severe = parse_fn(y$critical_severe_age, along = list(x = y$ages))
  
  # Parse function: age-related probabiliy of death given person is critical
  death_critical = parse_fn(y$death_critical_age, along = list(x = y$ages))
  
  # Remove redundant variables
  y[c("severe_symptom_age", "critical_severe_age", "death_critical_age")] = NULL
  
  # ---- Age-related cumulative probability matrix ----
  
  # Proportion of all cases that are asymptomatic
  p_asym = rep(y$proportion_asymptomatic, length(y$ages))
  
  # Proportion of all cases that are mild - the rest require hospitalisation
  p_mild   = (1 - p_asym) * (1 - severe_symptom)
  p_severe = (1 - p_asym) * severe_symptom
  
  # Proportion of cases that are severe that will / will not become critical
  p_critical     = p_severe * pmin(critical_severe, 1)
  p_non_critical = p_severe * (1 - pmin(critical_severe, 1))
  
  # Scale probabilities of death for critical cases in / out of care
  p_death_care = pmin(death_critical * y$death_critical_icu, 1)
  p_death_home = pmin(death_critical * y$death_critical_non_icu, 1)
  
  # Split these by hospital-seeking behaviour
  p_seek_S   = p_non_critical * y$seek_hospital
  p_noseek_S = p_non_critical * (1 - y$seek_hospital)
  
  # Of the critical cases admitted to ICU, who will recover and who will die
  p_seek_C = p_critical * y$seek_hospital * (1 - p_death_care)
  p_seek_D = p_critical * y$seek_hospital * p_death_care
  
  # Of the critical cases outside of the care system, who will recover and who will die
  p_noseek_C = p_critical * (1 - y$seek_hospital) * (1 - p_death_home)
  p_noseek_D = p_critical * (1 - y$seek_hospital) * p_death_home
  
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
  y$prognosis$age = data.table(
    age         = y$ages,
    asym        = p_asym,
    mild        = p_mild,
    severe_care = p_seek_S,
    crit_care   = p_seek_C,
    dead_care   = p_seek_D,
    severe_home = p_noseek_S,
    crit_home   = p_noseek_C,
    dead_home   = p_noseek_D) %>%
    pivot_longer(cols = -age, 
                 names_to = "state") %>%
    setDT()
  
  # Remove redundant variables
  y[c("proportion_asymptomatic", "seek_hospital", 
      "death_critical_icu", "death_critical_non_icu")] = NULL
  
  # ---- Vaccine severity and dose/exposure-response relationships ----
  
  # Extract variant severity factor for all circulating variants
  variant_df = y$variants %>%
    select(variant = id, 
           severe  = severity) %>%
    mutate(crit = severe, 
           dead = severe)
  
  # INTERPRETATION: Reductions in ICU and death for those hospitalised...
  
  # Vaccine dose - response relationship
  dose_df = fn_response(x = y$dose_response, 
                        count = 0 : y$max_dose_count, 
                        name  = "vaccine_doses")
  
  # Expsoure - response relationship
  exposure_df = fn_response(x = y$exposure_response, 
                            count = 0 : y$max_infection_count, 
                            name  = "num_infections")
  
  # Severity risk factors considering variant, exposures, and vaccine doses
  severity_df = 
    expand_grid(variant_df, dose_df, exposure_df, 
                .name_repair = make.names) %>%
    pivot_longer(cols = c(-variant, -vaccine_doses, -num_infections), 
                 names_to = "state") %>%
    group_by(variant, vaccine_doses, num_infections, state) %>%
    summarise(value = prod(value)) %>%
    ungroup() %>%
    setDT()
  
  # Repeat the dataframe for care-seeking and non-care-seeking severe states
  y$prognosis$severity = 
    rbind(severity_df %>% mutate(state = paste0(state, "_care")), 
          severity_df %>% mutate(state = paste0(state, "_home"))) %>%
    arrange(variant, vaccine_doses, num_infections, state)
  
  # Remove redundant variables
  y[c("dose_response", "exposure_response")] = NULL
  
  return(y)
}

# ---------------------------------------------------------
# Format response dose/exposure-response relationships
# ---------------------------------------------------------
fn_response = function(x, count, name) {
  
  # Construct state-response dictionary
  responses = qc(icu_risk_factor, death_risk_factor)
  response_dict = setNames(qc(crit, dead), responses)
  
  # States referred to in response list
  states = response_dict[names(x)]
  
  # Construct datatable of response for each count
  count_df = as_named_dt(x, states) %>% 
    cbind(count) %>% 
    mutate(across(.cols = -count, 
                  .fns  = function(x) x ^ count)) %>% 
    select(!!name := count, all_of(unname(states)))
  
  return(count_df)
}

