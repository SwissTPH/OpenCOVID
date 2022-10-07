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
  
  # ---- Age-related prognosis probabilities ----
  
  # Proportion of all cases that are asymptomatic
  p_asym = rep(y$proportion_asymptomatic, length(y$ages))
  
  # Proportion of all cases that are mild - the rest require hospitalisation
  p_mild = (1 - p_asym) * (1 - severe_symptom)
  p_hosp = (1 - p_asym) * severe_symptom
  
  # Proportion of cases that are severe that will / will not become critical
  p_severe = p_hosp * (1 - pmin(critical_severe, 1))
  p_icu    = p_hosp * pmin(critical_severe, 1)
  
  # Of the critical cases admitted to ICU, who will recover and who will die
  p_crit = p_icu * (1 - death_critical)
  p_dead = p_icu * death_critical
  
  # Sanity check: Check that these sum to one and nothing weird is going on
  if (any(abs(p_asym + p_mild + p_severe + p_crit + p_dead - 1) > 1e-6))
    stop("Prognosis probabilities must sum to 1 for all age groups")
  
  # Combine probabilities into dataframe (age group x prognosis state)
  y$prognosis$age = data.table(age    = y$ages,
                               asym   = p_asym, 
                               mild   = p_mild, 
                               severe = p_severe, 
                               crit   = p_crit, 
                               dead   = p_dead) %>%
    pivot_longer(cols = -age, 
                 names_to = "state") %>%
    setDT()
  
  # Remove redundant variables
  y[c("proportion_asymptomatic")] = NULL
  
  # plot_df = y$prognosis$age %>%
  #   mutate(state = fct_inorder(state))
  # 
  # g = ggplot(plot_df, aes(x = age, y = value, fill = state)) +
  #   geom_area()
  
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
  y$prognosis$severity = 
    expand_grid(variant_df, dose_df, exposure_df, 
                .name_repair = make.names) %>%
    pivot_longer(cols = c(-variant, -vaccine_doses, -num_infections), 
                 names_to = "state") %>%
    group_by(variant, vaccine_doses, num_infections, state) %>%
    summarise(value = prod(value)) %>%
    ungroup() %>%
    setDT()
  
  # Remove redundant variables
  y[c("dose_response", "exposure_response")] = NULL
  
  # plot_df = y$prognosis$severity %>%
  #   mutate(state = factor(state, levels = unique(y$prognosis$age$state)))
  # 
  # g = ggplot(plot_df) +
  #   aes(x = num_infections,
  #       y = vaccine_doses,
  #       colour = value,
  #       size   = value) +
  #   geom_point() +
  #   facet_grid(variant~state)
  
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

