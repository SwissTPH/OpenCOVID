###########################################################
# PARAMETERS
#
# Assign an input vector of parameter values to a list, p, which
# is fed into the model itself (see model.R). Also assign any 
# fixed parameters, some of which are determined from actual data
# (e.g demographic parameters).
#
###########################################################

# ---------------------------------------------------------
# Assign calibrated and fixed parameters ready for model simulation
# ---------------------------------------------------------
get_parameters = function(o, d, calibration_params, sample_fixed = FALSE, do_plot = FALSE) {
  
  # Sanity check: data must be provided
  if (is.null(d))
    stop("No data provided - check function call")
  
  # ---- Preprocess calibration parameters ----
  
  # If nothing supplied, use prior means 
  if (is.null(calibration_params)) {
    
    # Extract prior means and convert to named vector
    calibration_params = o$calibration_df$values
    names(calibration_params) = o$calibration_df$names
  }
  
  # Remove canton name from any non-global parameters
  names(calibration_params) = str_remove(names(calibration_params), "\\.[A-Z]*")
  
  # Shorthand for parameter names (after removing canton reference)
  calibration_names = names(calibration_params)
  
  # ---- Extract all parameters ----
  
  # Initiate parameters list - this is fed into the model
  p = list()
  
  # Then store value of calibrated parameters
  p[calibration_names] = as.numeric(calibration_params)
  
  # Then store the value of fixed parameters
  for (i in seq_len(nrow(o$fixed_df))) {
    fixed  = o$fixed_df[i, ]
    f_name = fixed$names
    
    # Use predefined values when calibrating and for best estimates
    if (!sample_fixed | fixed$prior_sd == 0) 
      p[[f_name]] = fixed$values
    
    else { # When running uncertainty analysis we want to sample from fixed distribution
      
      stop("TODO...") # Simply use normal distribution as standard ...
      
      # Get and apply some distribution
      fixed_dist = get(fixed$prior_dist)
      fixed_val  = fixed_dist(1, fixed$values, fixed$prior_param)
      
      # Ensure this is within the bounds for this parameter
      p[[f_name]] = max(min(fixed_val, fixed$upper_bound), fixed$lower_bound)
    }
  }
  
  # ---- Check that we have what we need ----
  
  # Check that we don't have any repeated values in the input parameter vector
  if (!identical(calibration_names, unique(calibration_names)))
    stop("Non unique parameters provided - investigate 'calibration_params' input")
  
  # Sanity check on inputs provided - should align with loaded calibration table
  if (!identical(sort(names(p)), sort(o$parameter_table$names))) {
    
    # Throw an error if parameters are missing
    # missing_params = setdiff(o$parameter_table$names, names(p))
    # if (length(missing_params) > 0)
    #   stop("Parameters not accounted for: ", paste(missing_params, collapse = ", "))
    
    # Throw an error if any parameters are unrecognised
    unknown_params = setdiff(names(p), o$parameter_table$names)
    if (length(unknown_params) > 0)
      stop("Parameters not recognosied: ", paste(unknown_params, collapse = ", "))
    
    # If neither of the above, something strange is going on
    # stop("Parameters not well defined - investigation needed")
  }
  
  # ---- Model flows ----
  
  # Store model flows in p so we don't need to pass o around in model
  p$model_flows = o$model_flows
  
  # Once diagnosed, only those with certain prognoses will go into isolation
  p$iso_prognosis = unique(filter(p$model_flows, care_state == "iso")$prognosis_state)
  
  # ---- National demographics ----
  
  # Simply append demographic breakdown by age
  p$demography = d$demog
  
  # Determine population scaler required to satisfy maximum population size
  #
  # NOTE: If population_size <= o$max_population_size then population_scaler equals one
  p$population_scaler = ceiling(sum(d$demog) / o$max_population_size)
  
  # Total population size
  p$population_size = round(sum(d$demog) / p$population_scaler)
  
  # ---- Risk groups ----
  
  # TODO: Remove hardcoding and call to load_data_pop here...
  
  # Get Swiss population size
  swiss_pop = 8606033 # sum(load_data_pop(o, force_national = TRUE)$value)
  
  # Calculate proportion of risk group of total swiss pop
  p$risk_group_ratio <- list(young_comorbidities = 100000 / swiss_pop, 
                             old_comorbidities   = 621600 / swiss_pop, 
                             healthcare_workers  = 560000 / swiss_pop)
  
  # ---- Importation (initial, constant, and new variants) ----
  
  # Interpret first import to be # days before epidemic outbreak	
  p$import_date_initial = max(1, d$outbreak_start - p$import_date)
  
  # Logical vector of when to introduce these new variants
  p$import_date_variant = o$dates_all %in% o$variants$import_date
  
  # Scale imports per 100,000
  #
  # NOTE: No need to population-scale import_constant as this occurs organically in the model
  p$import_initial  = p$import_initial / 1e5 * p$population_size
  p$import_constant = p$import_constant / 1e5
  
  # TODO: Need to make this more robust (multiple variants fitted at once, remove calibration param
  #       when all are fixed in variant_properties.csv etc etc)
  
  # Number of people to import can be fixed or calibrated (if no fixed value provided)
  import_variant = o$variants$import_number
  import_variant[is.na(import_variant)] = p$import_variants
  
  # Absolute number of people to import with new variants
  p$import_variant = import_variant / 1e5 * p$population_size
  
  # ---- Viral variant properties ----
  
  # Variant properties defined in input file (see variant_properties.csv)
  p$variant_infectivity = o$variants$infectivity_factor
  p$variant_severity    = o$variants$severity_factor  # Currently interpreted as increased risk of death
  
  # Update variant infectivity only if selected for calibration
  if ("variant_infectivity_factor" %in% calibration_names)
    p$variant_infectivity[o$variants$calibrated] = p$variant_infectivity_factor
  
  # Update variant severity only if selected for calibration
  if ("variant_severity_factor" %in% calibration_names)
    p$variant_severity[o$variants$calibrated] = p$variant_severity_factor
  
  # ---- Viral load profile ----
  
  # TODO: Infectious dose, possibly age, and other stuff should contribute here, 
  #       in effect these properties should alter the distribution parameters...
  
  # Define max possible viral shedding days to be period of infectiousness for average severe case
  days_shed = p$presymptomatic_days + p$infectious_days_severe
  
  # Either a constant daily average (multiplied by beta in model)
  if (o$viral_load_infectivity == FALSE) {
    p$viral_load_profile = rep(1 / days_shed, days_shed)
    
  } else {  # Or a profile represnted by viral load...
    
    # TODO: Move distribution values to parameters config excel file...
    
    # Curve over time assumed to follow a gamma distrubtion, peaks around day of symptom onset
    viral_load_shape = dgamma(1 : days_shed, 3, rate = 0.5)
    
    # Normalise so we have the peak at 1 (multiplied by beta in model)
    p$viral_load_profile = viral_load_shape / max(viral_load_shape)
  }
  
  # ---- Seasonality profile ----
  
  # Extrapolation dates - everything beyond the last data date
  extrap_dates = o$dates_all[o$dates_all > max(o$dates_data)]
  
  # Smoothly interpolate between monthly means over a calendar year
  extrap_year = splineInterpolateMonthlytoDailyforSeveralYears(d$weather_extrap[, "temperature"], 
                                                               start_year = year(min(extrap_dates)))
  
  # Index this vector for future dates and concatenate with daily data
  temperature_future = extrap_year[yday(extrap_dates)]
  temperature_extrap = c(d$weather$temperature, temperature_future)
  
  # TODO: Check nothing goes awry when we have negative values here...
  
  # Inverse temperature (so peak is during coldest period) and normalise
  temperature_shift   = temperature_extrap - min(temperature_extrap)
  temperature_inverse = 1 - temperature_shift / max(temperature_shift)
  
  # Scale seasonality effect based on inverse temprerature and seasonality_scaler
  p$seasonality = 1 + p$seasonality_scaler * (temperature_inverse - 1)
  
  # Store a few vectors - simply for plotting purposes
  p$temperature_extrap  = temperature_extrap
  p$temperature_inverse = temperature_inverse
  
  # ---- Prognosis track probabilities ----
  
  # TODO: This needs an overhaul using Swiss data
  
  # NOTES:
  # - Consider using a single parameter to define these growths/decays rather than all this hardcoding
  # - Severe and critical cases can change prognosis track when hospital capacity has been breached (NOT YET IMPLEMENTED)
  # - Source: https://www.imperial.ac.uk/media/imperial-college/medicine/sph/ide/gida-fellowships/Imperial-College-COVID19-NPI-modelling-16-03-2020.pdf
  
  # Proportion of symptomatic cases that are severe and require hospitalisation (by age group)
  p$symptom_severe_age = c(0.1, 0.3, 1.2, 3.2, 4.9, 10.2, 16.6, 24.3, 27.3) / 100
  
  # Proportion of severe cases that will become critical and require care in ICU (by age group)
  p$severe_critical_age = c(5, 5, 5, 5, 6.3, 12.2, 27.4, 43.2, 70.9) / 100
  
  # TEMP: The very elderly are less likely to seek care
  p$elderly_seek_factor = 0.5
  
  # Age scaler for probability of death once in critical state (both in and out of ICU)
  death_decay = dexp(1 : length(o$age_groups), rate = 1.2)
  p$death_critical_age = rev(death_decay) / max(death_decay)
  
  # TEMP: Total probability of death by age group (from Swiss mortality data - use this for sanity checking only)
  #
  # TODO: Calculate this from loaded data and remove this hardcoding
  p$death_age = c(0.01, 0.01, 0.01, 0.10, 0.27, 1.61, 5.85, 20.08, 72.06) / 100
  
  # Reduction in severe cases becoming criticial cases due to improved care and treament procedures (over time)
  p$improved_care_scaler = rep(1, o$n_dates)
  p$improved_care_scaler[o$dates_all >= o$improved_care_date] = p$improved_care_factor
  
  # ---- Disease and care state durations ----
  
  # TODO: Extracting durations can be generalised and applied in a loop. Surely.
  
  # Duration of latent period: infected not infectious
  # p$duration_latency = function(id) {pmax(1, round(rgamma(length(id), p$latency_days, rate = 1)))}
  p$duration_latency = function(id) {pmax(1, round(rnorm(length(id), p$latency_days, sd = 1)))}
  
  # Duration of pre-symptomatic period: infectious without symptoms
  p$duration_presymptomatic = function(id) {pmax(1, round(rnorm(length(id), p$presymptomatic_days, sd = 1)))}
  
  # Duration of severe and non-severe infectious periods
  p$duration_non_severe = function(id) {pmax(1, round(rnorm(length(id), p$infectious_days_mild,   sd = 1)))}
  p$duration_severe     = function(id) {pmax(1, round(rnorm(length(id), p$infectious_days_severe, sd = 1)))}
  
  # Duration of period between symptom onset and hospitalisation for severe cases
  p$duration_onset_to_hospital = function(id) {pmax(1, round(rnorm(length(id), p$onset_to_hospital_days, sd = 4)))}
  
  # Other hospital flow duration distributions
  p$duration_hospital_stay     = function(id) {pmax(1, round(rnorm(length(id), p$hospital_stay_days,     sd = 1)))}
  p$duration_hospital_to_icu   = function(id) {pmax(1, round(rnorm(length(id), p$hospital_to_icu_days,   sd = 1)))}
  p$duration_icu_stay          = function(id) {pmax(1, round(rnorm(length(id), p$icu_stay_days,          sd = 1)))}
  p$duration_icu_stay_death    = function(id) {pmax(1, round(rnorm(length(id), p$icu_stay_death_days,    sd = 1)))}
  p$duration_hospital_transfer = function(id) {pmax(1, round(rnorm(length(id), p$hospital_transfer_days, sd = 1)))}
  
  # Time until death for critical cases outside of the hospital setting
  p$duration_home_death = function(id) {pmax(1, round(rnorm(length(id), p$home_death_days, sd = 1)))}
  
  # Most common infectious period - used to calculate effective reproduction number
  p$infectious_period = p$presymptomatic_days + p$infectious_days_mild
  
  # ---- Reporting delays ----
  
  # Preallocate vectors for all metrics (although reporting delays do not apply to most)
  p$reporting_delays = setNames(rep(0, nrow(o$all_metrics)), o$all_metrics$metric)
  
  # All the various types of reporting delays to consider
  all_delays = str_remove(na.omit(o$all_metrics$delay), " ")
  all_delays = unique(unlist(str_split(all_delays, ",")))
  
  # Loop through the various types of reporting delays
  for (this_delay in all_delays) {
    
    # TODO: The rounding is not ideal here...
    
    # All metrics associated with this delay, and value of this delay
    this_delay_idx = grepl(this_delay, o$all_metrics$delay)
    this_delay_val = round(p[[paste0(this_delay, "_delay")]])
    
    # Increment delays (some metrics may have multiple delays)
    p$reporting_delays[this_delay_idx] = 
      p$reporting_delays[this_delay_idx] + this_delay_val
  }
  
  # Remove trivial values 
  # p$reporting_delays = p$reporting_delays[p$reporting_delays > 0]
  
  # ---- Testing and diagnosis ----
  
  # Extend number of diagnosis vector - keep NAs in the future as we calculate as we go in model.R
  diagnosed_df = data.table(date = o$dates_all) %>%
    left_join(d$diagnoses, by = "date")
  
  # Store extrapolated number diagnosed in parameter list
  p$number_diagnosed = diagnosed_df$number_diagnosed
  
  # Named vector for testing and diagnosis priority (worst disease state to best)
  p$diagnosis_priority = setNames(1 : length(o$disease_states), rev(o$disease_states))
  
  # We guarantee that the most severe cases will receive a test and diagnosis
  p$diagnosis_guarantee = rev(o$disease_states)[1 : which(rev(o$disease_states) == "severe")]
  
  # ... but how should this interact with those tested via contact tracing?
  p$contact_tracing = rep(FALSE, o$n_dates)  # TEMP: Trivialised for now
  
  # ---- Vaccine properties ----
  
  # TODO: 
  #  1) We may want a baseline extrapolation for number_vaccines
  #  2) It may also be a good idea to do some data smoothing
  #  3) number_vaccines could be a matrix for rows for number of EACH VACCINE to distribute every day
  
  # Extract number of vaccines distributed so far
  vaccine_df = data.table(date = o$dates_all) %>%
    left_join(d$vaccine, by = "date") %>%
    replace_na(list(number_vaccines = 0))
  
  # NOTE: These can be altered into the future for different scenarios ...
  
  # Number of vaccines to distribute
  #
  # NOTE: Population scaling occurs within model.R
  p$number_vaccines = vaccine_df$number_vaccines
  
  # Vaccine coverage per priority group
  p$vaccine_coverage = matrix(0, nrow = o$n_dates, ncol = length(o$vaccine_priority))
  colnames(p$vaccine_coverage) = names(o$vaccine_priority)
  
  # Store default vaccine in parameters list
  p$vaccine_default    = o$vaccine_default
  p$vaccine_acceptance = o$vaccine_acceptance
  
  # Define which vaccine to be used for each group
  vaccine_type   = rep(o$vaccine_default, length(o$vaccine_priority))
  p$vaccine_type = setNames(vaccine_type, names(o$vaccine_priority))
  
  # Vaccine effect is currently coming from options - want to proper parameterise this in the future
  p$vaccine_effect   = o$vaccine_effect
  p$vaccine_criteria = o$vaccine_criteria
  
  # Similar for transmission blocking and disease prevention scalers
  if (p$vaccine_effect == "mixed") {
    p$transmission_blocking_scaler = o$transmission_blocking_scaler
    p$disease_prevention_scaler    = o$disease_prevention_scaler
    
  } else {  # If vaccine_effect is not mixed, trivialise these parameters
    p$transmission_blocking_scaler = 1
    p$disease_prevention_scaler    = 1
  }
  
  # Sanity check on the value of vaccine_effect
  if (!p$vaccine_effect %in% c("transmission", "severity"))
    stop("Vaccine effect '", p$vaccine_effect, "' not recognised")
  
  # Vaccine efficacy comes straight from pre-formated data
  p$vaccine_efficacy = o$vaccine_efficacy
  p$n_doses          = o$n_doses
  
  # TODO: Consider immunity decay over time here too...
  
  # TODO: Move these to vaccine properties file so they can be unique for each vaccine
  immunity_growth = 5 # Growth rate for immunity following vaccination
  immunity_mid = p$vaccine_immunity_days / 2 # Inflection point is half the period until full immunity
  
  # Loop though the vaccines we're modelling
  vaccine_profiles = list()
  for (vaccine_name in names(o$vaccine_dict)) {
    
    # Properties that influence immunity / severity
    this_efficacy = p$vaccine_efficacy[[vaccine_name]]
    
    # Generate a logistic function to define temporal dynamics
    temporal_profile = logistic(1 : o$n_dates, 
                                slope = immunity_growth, 
                                mid   = immunity_mid, 
                                upper = this_efficacy)
    
    # Store in a list that'll be converted to datatable
    vaccine_profiles[[vaccine_name]] = temporal_profile
  }
  
  # Store temporal dynamics of all vaccines in datatable
  p$vaccine_profile = as.data.table(vaccine_profiles)
  
  # ---- Non-pharmaceutical interventions ----
  
  # NOTE: These can be altered into the future using alter_model_params (see scenarios.R)
  
  # Project into the future by assuming most recent policy remains
  response_df = data.table(date = o$dates_all) %>%
    left_join(d$response, by = "date") %>%
    fill(osi_level)
  
  # Store extrapolated NPI effect in parameter list
  p$npi_level = pmin(1, response_df$osi_level * p$npi_scaler)
  
  # Allow potential to examine effect of NPI adherence
  p$npi_adherence = rep(1, o$n_dates)
  
  # Alter adherence at certain time points
  p$npi_adherence[o$dates_all >= "2020-12-22"] = p$holiday_adherence
  
  # ---- Capacity constraints ----
  
  # If/when this capacity data is temporal, use this...
  
  # # Project into the future by assuming current capacity stays constant
  # icu_capacity_df = data.table(date = o$dates_all) %>%
  #   left_join(d$capacity, by = "date") %>%
  #   fill(icu_capacity, .direction = "downup")
  # 
  # # Store ICU capacity for this canton
  # p$icu_capacity = icu_capacity_df$icu_capacity
  
  # ... unitl then, just use it directly ...
  p$icu_capacity = d$capacity
  
  # ---- Remove redundant fields ----
  
  # These have done their job - remove to avoid confusion when parameters are actually used in the model
  
  # Miscellaneous parameters
  p[c("import_date", "import_variants")] = NULL
  
  # Infectiousness-related parameters
  p[c("beta_reduction", "contacts_reduction", "susceptibility")] = NULL
  
  # Delay-related parameters
  p[c("diagnosis_delay_care", "diagnosis_delay_home", 
      "diagnosis_delay_mild", "death_reporting_delay")] = NULL
  
  # Disease, death, and variant related parameters
  p[c("improved_care_factor", "variant_infectivity_factor", "variant_severity_factor")] = NULL
  # p[c("severe_factor", "critical_factor", "seek_hospital", 
  #     "critical_death_icu", "critical_death_non_icu")] = NULL
  
  # Finally, order elements alphabetically to make it easy to search
  p = p[order(names(p))]
  
  # ---- Plotting ---- 
  
  # Produce any plots if desired
  if (do_plot == TRUE) {
    
    # Call plotting functions (see plotting.R)
    plot_viral_load(o, p)  # Viral load profile
    plot_durations(o, p)   # Duration distribution
    plot_seasonality(o, p, d)    # Seasonality effect
    plot_vaccine_efficacy(o, p)  # Vaccine efficacy profiles e
  }
  
  return(p)
}

# ---------------------------------------------------------
# Transform model parameters for use in calibration and back again
# ---------------------------------------------------------
transform_params = function(o, values, rev = FALSE, is_named = FALSE) { 
  
  # Sanity check on the length of input values
  if (!is_named & nrow(o$calibration_df) != length(values))
    stop("Inconsistent parameter vector length")
  
  # We can also apply on a canton subset, but the vector must be named 
  if (is_named) param_names = names(values)
  else param_names = o$calibration_df$names
  
  browser() # TODO: Update this function with simplier transformation approach
  
  # Initiate vector for transformed values
  transformed_params = as.numeric(values)
  
  # Loop through parameters
  for (i in 1 : length(param_names)) {
    param = o$calibration_df[names == param_names[i]]
    
    # We can either transform for calibration or back again to run model
    if (rev) trans_dist = param$retrans_dist
    else     trans_dist = param$trans_dist
    
    # The trivial case is defined by 'none'
    if (tolower(trans_dist) != "none") {
      
      # Get the transformation function
      trans_func = get(trans_dist)
      
      # Apply the tranformation to parameter best estimate
      transformed_params[i] = trans_func(as.numeric(values[i]))
    }
  }
  
  return(transformed_params)
}

