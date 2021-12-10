###########################################################
# PARSE INPUT 
#
# Load yaml input file of interest and parse all inputs such
# that they can be interpreted by the transmission model.
#
# Note that it can be expensive to parse large, high-dimensional
# arrays, so by default read_array is set to FALSE.
#
###########################################################

# ---------------------------------------------------------
# Read yaml input file and perform some additional parsing
# ---------------------------------------------------------
parse_yaml = function(o, scenario, read_array = FALSE) {
  
  # NOTE: The initial parsing work is done by read_yaml from yaml package
  
  # Paths to default and user-defined yaml input files
  yaml_file_default = file.path(o$pth$config, "default.yaml")
  yaml_file_user    = file.path(o$pth$input, paste0(o$analysis_name, ".yaml"))
  
  # Check that the user-defined file actually exists 
  if (!file.exists(yaml_file_user))
    stop("\n You are attempting to run analysis: ", o$analysis_name, 
         "\n  but no yaml file was found with this name", 
         "\n  (missing file: ", yaml_file_user, ")")
  
  # Load default model parameters
  y = read_yaml(yaml_file_default)
  
  # Load user-defined model parameters
  y_user = read_yaml(yaml_file_user)
  
  # Overwrite any parameter defaults for which we have user-defined values 
  #
  # NOTE: A few sanity checks on the user-defined inputs are perfomed here
  y = overwrite_defaults(y, y_user)
  
  # ---- Apply scenario of choice ----
  
  # Apply values provided within specified scenario block
  #
  # NOTE: Some general checks on the content of all scenarios are perfomed here
  y = parse_scenarios(o, y, scenario, read_array)
  
  # If we only need to read scenario names, return out here with named vector
  if (scenario == "*read*")
    return(y$scenario_names)
  
  # ---- Time ----
  
  # Easier to consider n_days_init as a positive integer
  n_days_init = -y$n_days_init
  if (n_days_init < 0)
    stop("It is an error to initialise the epidemic in the future, set n_days_init < 0")
  
  # Some parameters need to be length of initial plus projection period
  n_days_total = y$n_days + n_days_init
  
  # ---- Demographics ----
  
  # Hardcode maximum upper age to 90
  age_max = 90
  
  # We're not yet interpolating between ages - throw an error if insufficiant age bins
  if (length(y$demography) != age_max)
    stop("Single-year age bins must be provided from year 1 up to year ", age_max)
  
  # Multiple age pyramid by population size to model
  y$demography = unlist(y$demography) / sum(unlist(y$demography))
  y$demography = round(y$demography * y$population_size)
  
  # Correct any rounding errors by adding/subtracting from first age bin
  round_err = sum(y$demography) - y$population_size
  y$demography[1] = y$demography[1] - round_err
  
  # Age structure (in years) - defines upper bound of bin
  y$age = list(all = 1 : age_max - 1, groups = seq(10, age_max, by = 10))  
  
  # ---- Network ----
  
  # Sanity check on network structure
  if (!y$network_structure %in% qc(random, age, layers))
    stop("Network structure '", y$network_structure, "' not recognised")
  
  # Network layers are not quite ready yet - they will be soon!
  if (y$network_structure == "layers")
    stop("Network 'layers' are currently under development", 
         " - these will be available in version 3.0")
  
  # Convert comma-seperated string to vector of strings
  y = str2vec(y, "network_layers")
  y = str2vec(y, "contact_matrix_countries")
  
  # Sanity check on network layers if needed
  if (y$network_structure == "layers") {
    
    # Network layers recognised / unrecognised
    known_layer    = y$network_layers %in% qc(household, school, workplace)
    unknown_layers = paste(y$network_layers[!known_layer], collapse = ", ")
    
    # Throw an error if there are any unrecognised layers
    if (unknown_layers != character(1))
      stop("Unrecognised network layer(s): '", unknown_layers)
  }
  
  # Throw an error if country codes are in not in ISO-2 format
  if (any(sapply(y$contact_matrix_countries, nchar) != 2))
    stop("Use ISO Alpha-2 country codes for 'contact_matrix_countries' parameters")
  
  # ---- Risk groups ----
  
  # Convert nested list into datatable
  y$risk_groups = y$risk_groups %>% 
    lapply(as.data.table) %>% 
    rbindlist() %>%
    mutate(name = names(y$risk_groups), .before = 1)
  
  # Sanity check that upper age is no more than age_max
  if (any(y$risk_groups$age_upper > age_max))
    stop("Invalid risk group age range")
  
  # ---- Importation ----
  
  # Scale initial and constant imports per 100,000
  y$import_constant = y$import_constant / 1e5
  
  # ---- NPI effect ----
  
  # Most basic case: NPI effect is a single value moving forward
  if (is.numeric(y$npi_effect) && length(y$npi_effect) == 1) {
    
    # Convert to a vector - one value per time point
    y$npi_effect = rep(y$npi_effect, y$n_days)
    
    # Consider special case: load OCHI from OSI for given country
  } else if ("fn" %in% names(y$npi_effect)) {
    
    message("  > Pulling OSI data")
    
    # Dates to load data for (only go as far as yesterday - we'll fill the rest)
    date_from = format_date(y$npi_effect$start)
    date_to   = min(date_from + y$n_days - 1, format_date(today() - 1))
    
    # Construct call to API endpoint and convert to json format
    api_call = paste0(o$osi_api, date_from, "/", date_to)
    api_data = fromJSON(rawToChar(httr::GET(api_call)$content))
    
    # Data is horribly hidden, a few lapply's needed to extract what we need
    osi_fn = function(x) rbindlist(lapply(x, as.data.table), fill = TRUE)
    osi_df = rbindlist(lapply(api_data$data, osi_fn), fill = TRUE)
    
    # Filter for only this country
    osi_country = osi_df %>% 
      filter(country_code == y$npi_effect$country) %>% 
      arrange(date_value) %>%
      mutate(day = 1 : n(), 
             stringency = stringency / 100) %>%
      select(day, stringency)
    
    # Extract stringency and extend if necessary
    y$npi_effect = data.table(day = 1 : y$n_days) %>%
      left_join(osi_country, by = "day") %>%
      fill(stringency) %>%
      pull(stringency)
  }
  
  # ---- Seasonality ----
  
  # Parse seasonality evaluation points
  y$seasonality_fn$x = parse_fn(y$seasonality_fn$x)
  
  # Evaluate seasonality function
  seasonality_profile = parse_fn(y$seasonality_fn)
  
  # Check this generates a vector and has 365 values
  if (length(seasonality_profile) != 365)
    stop("Seasonality profile must be a vector of 365 values")
  
  # Inverse seasonality (so peak is during coldest period) and normalise
  seasonality_trans   = seasonality_profile - min(seasonality_profile)
  seasonality_inverse = 1 - seasonality_trans / max(seasonality_trans)
  
  # Scale seasonality effect based on inverse temprerature and seasonality_scaler
  seasonality_scale = 1 + y$seasonality_scaler * (seasonality_inverse - 1)
  
  # We may want to temporally 'shift' so we can initiate at a particular time of year
  seasonality_shift = c(seasonality_scale[(y$seasonality_shift + 1) : 365], 
                        seasonality_scale[seq_len(y$seasonality_shift)])
  
  # Extend (or clip) scaled seasonality profile to appropriate number of days
  y$seasonality = rep_len(seasonality_shift, y$n_days)
  
  # Also generate seasonality curve for initialisation period
  y$seasonal_init = rev(rep_len(rev(seasonality_shift), n_days_init))
  
  # Remove redundant items
  y[c("seasonality_fn", "seasonality_scaler", "seasonality_shift")] = NULL
  
  # ---- Prognosis probabilities ----
  
  # Parse function: age-related probabiliy of severe disease given person is symptomatic
  y$severe_symptom_age  = parse_fn(y$severe_symptom_age,  along = list(x = y$age$groups))
  
  # Parse function: gge-related probabiliy of critical disease given person is severe
  y$critical_severe_age = parse_fn(y$critical_severe_age, along = list(x = y$age$groups))
  
  # Parse function: age-related probabiliy of death given person is critical
  y$death_critical_age  = parse_fn(y$death_critical_age,  along = list(x = y$age$groups))
  
  # ---- Disease state durations ----
  
  # Define max possible viral shedding days to be period of infectiousness for average severe case
  n_shed_days = y$durations$presymptomatic$mean + y$durations$infectious_severe$mean
  
  # Most common infectious period - used to calculate effective reproduction number
  y$infectious_period = y$durations$presymptomatic$mean + y$durations$infectious_mild$mean
  
  # Initiate list to store parsed functions
  durations_fn = list()
  
  # Loop through durations we want to construct funcions for
  for (this_duration in names(y$durations)) {
    
    # Parse the function call, including default first argument (# of people)
    parsed_fn = parse_fn(fn_args  = y$durations[[this_duration]], 
                         along    = list(n = "length(id)"), 
                         evaluate = FALSE)
    
    # Bound below by one day and round to nearest integer
    bounded_fn = paste0("pmax(1, round(", parsed_fn, "))")
    
    # Create a string which - once evaluated - will assign the function to the durations_fn list
    assign_fn = paste0("durations_fn$", this_duration, " = function(id) {", bounded_fn, "}")
    
    # Evaluate the function assignment
    eval(parse(text = assign_fn))
  }
  
  # Overwrite list with corresponding functions
  y$durations = durations_fn
  
  # ---- Waning immunity ----
  
  # Initiate a temporal vector of acquired immunity effect over time
  a = y$acquired_immunity
  a$profile = rep(a$value, n_days_total)
  
  # If simulating long enough, also consider acquired immunity decay
  if (n_days_total > a$decay_init) {
    
    # Linear decay from initial immunity to zero
    decay_vec = seq(a$profile[a$decay_init], 0, length.out = a$decay_duration)
    
    # Apply this decay - clipping vector if necessary
    decay_idx = 1 : min(a$decay_duration, n_days_total - a$decay_init)
    a$profile[a$decay_init + decay_idx] = decay_vec[decay_idx]
    
    # After this point, assume acquired immunity is zero
    if (max(decay_idx) == a$decay_duration)
      a$profile[(a$decay_init + max(decay_idx) + 1) : n_days_total] = 0
  }
  
  # Overwrite acquired_immunity item with temporal profile
  y$acquired_immunity = a$profile
  
  # ---- Testing and diagnosis ----
  
  # Values to evaluate testing probabilities for - normalised ages
  along_vals = list(x = y$age$all / max(y$age$all))
  
  # Loop through the different types of testing
  for (test_type in names(y$testing)) {
    test_info = y$testing[[test_type]]
    
    # Evaluate the user-defined function 
    age_dist = parse_fn(fn_args = test_info$age_dist, along = along_vals)
    
    # Normalise to a mean of one then multiply through by probability AND sensitivity
    test_info$age_prob = (age_dist / mean(age_dist)) * 
      test_info$probability * test_info$sensitivity
    
    # Remove redundant items and store age-probability values
    y$testing[[test_type]] = list.remove(test_info, "age_dist")
  }
  
  # Parse 'when' mass testing events should happen
  y$testing$mass_testing$when = parse_fn(fn_args = y$testing$mass_testing$when)
  
  # ---- Viral load profile ----
  
  # Either a constant daily average (multiplied by beta in model)
  if (y$viral_load_infectivity == FALSE) {
    y$viral_load_profile = rep(1 / n_shed_days, n_shed_days)
    
  } else {  # Or a profile represnted by viral load...
    
    # Evaluate user-defined function for all days for which viral shedding is possible
    viral_load_shape = parse_fn(fn_args = y$viral_load_shape_fn, 
                                along   = list(x = paste("1", ":", n_shed_days)))
    
    # Normalise so we have the peak at 1 (multiplied by beta in model)
    y$viral_load_profile = viral_load_shape / max(viral_load_shape)
  }
  
  # Remove redundant items
  y[c("viral_load_infectivity", "viral_load_shape_fn")] = NULL
  
  # ---- Viral variants ----
  
  # Append trivial default values to primary variant - infectivity & severity always relevant to this
  variant_primary = append(y$variant_primary, list(import_day = 0, infectivity = 1, severity = 1))
  
  # Convert lists into single dataframe
  y$variants = list2dt(y$variants_novel) %>%
    filter(import_number > 0) %>%
    mutate(import_number = ceiling(import_number / 1e5 * y$population_size)) %>%
    bind_rows(variant_primary) %>%
    arrange(import_day)
  
  # Remove redundant items
  y[c("variant_primary", "variants_novel")] = NULL
  
  # ---- Vaccine properties ----
  
  # If using a custom vaccine, apply user-defined vaccine properties
  if (y$vaccine_type == "custom") {
    v = append(list(id = y$vaccine_type), y$vaccine_custom)
    
  } else {  # ... or extract a named vaccine
    v = append(list(id = y$vaccine_type), y$vaccine_defintion[[y$vaccine_type]])
  }
  
  # Transmission blocking effect cannot be more than total efficacy
  if (v$transmission_blocking >= v$efficacy)
    stop("Vaccine transmission blocking effect is greater than total efficacy")
  
  # The remaining effect comes from infections that do not result in severe disease
  v$disease_prevention = v$efficacy - v$transmission_blocking
  
  # Remove redundant items
  y[c("vaccine_type", "vaccine_defintion", "vaccine_custom")] = NULL
  
  # ---- Vaccine profile ----
  
  # Evaluate vaccine immunity growth function
  growth_shape = parse_fn(fn_args = v$growth_fn, 
                          along   = list(x = paste("1", ":", n_days_total)))
  
  # Multiply through by maximum efficacy for vaccine immunity scale up
  v$profile = growth_shape * v$efficacy
  
  # Preallocate for booster dose profile
  v$booster = rep(NA, n_days_total)
  
  # Leave NA's for growth period - need to calculate these as we go
  v$booster[1 : n_days_total] = v$booster_efficacy
  
  # If simulating long enough, also consider vaccine immunity decay
  if (n_days_total > v$decay_init) {
    
    # Linear decay from current immunity to zero
    decay_vaccine = seq(v$profile[v$decay_init], 0, length.out = v$decay_duration)
    decay_booster = seq(v$booster[v$decay_init], 0, length.out = v$decay_duration)
    
    # Apply this decay - clipping vector if necessary
    decay_idx = 1 : min(v$decay_duration, n_days_total - v$decay_init)
    v$profile[v$decay_init + decay_idx] = decay_vaccine[decay_idx]
    v$booster[v$decay_init + decay_idx] = decay_booster[decay_idx]
    
    # After this point, assume immunity is zero
    if (max(decay_idx) == v$decay_duration) {
      v$profile[(v$decay_init + max(decay_idx) + 1) : n_days_total] = 0
      v$booster[(v$decay_init + max(decay_idx) + 1) : n_days_total] = 0
    }
  }
  
  # Rename list item to make it clear that this is a temporal profile
  names(v)[names(v) == "booster"] = "booster_profile"
  
  # Remove redundant items
  v[c("growth_fn")] = NULL
  
  # ---- Vaccine rollout ----
  
  # Join vaccine history and rollout tables
  v$groups = list2dt(y$vaccine_history) %>% 
    full_join(list2dt(y$vaccine_rollout), by = "id") %>%
    mutate(priority = as.integer(priority)) %>%
    arrange(priority)
  
  # Check for any missing values - will occur if priority groups do not align
  if (any(is.na(v$groups)))
    stop("Non-conformable vaccine history and vaccine rollout defintions")
  
  # Sanity check: n_days_init should be greater than all vaccine init start dates
  if (any(v$groups$init_start < y$n_days_init))
    stop("It is an error to initiate vaccination before epidemic outbreak")
  
  # Sanity check: init dates should be neagtive
  if (any(v$groups[, .(init_start, init_end)] > 0))
    stop("Vaccine history 'init_start' and 'init_end' must be non-positive")
  
  # Sanity check: end date not before start date
  if (any(v$groups[, init_start > init_end]))
    stop("Vaccine history 'init_start' cannot be greater than 'init_end'")
  
  # Sanity check: 
  if (any(v$groups[init_start == init_end, init_coverage] > 1e-6))
    stop("Vaccine history 'init_coverage' must be zero if 'init_start' == 'init_end'")
  
  # Sanity check: mas scale-up dates should positive
  if (any(v$groups[, .(max_start, max_end)] < 0))
    stop("Vaccine rollout 'max_start' and 'max_end' must be non-negative")
  
  # Sanity check: end date not before start date
  if (any(v$groups[, max_start > max_end]))
    stop("Vaccine rollout 'max_start' cannot be greater than 'max_end'")
  
  # Sanity check: Max coverage should be at least init coverage
  if (any(v$groups[, init_coverage - max_coverage] > 1e-6))
    stop("Vaccine 'init_coverage' cannot be greater than 'max_coverage'")
  
  # Sanity check: We cannot scale-up further if there is no time to do so
  if (any(v$groups[max_start == max_end, max_coverage - init_coverage] > 1e-6))
    stop("Vaccine 'max_coverage' must be equal to 'init_coverage' if 'max_start' == 'max_end'")
  
  # Append number of vaccine groups and initial vaccine coverage
  v$n_groups = nrow(v$groups)
  
  # Ensure the priority ranking makes sense - require groups are defined 1 : n
  if (!identical(v$groups$priority, 1 : v$n_groups))
    stop("Priority of vaccine groups must be ranked 1 : n")
  
  # Vaccine booster details - totally fine to have NAs here
  v$booster_groups = v$groups[, .(id)] %>% 
    left_join(list2dt(y$booster_rollout)[probability > 0], by = "id") %>%
    rename(vaccine_group = id) %>%
    mutate(probability = pmin(pmax(probability, 0, na.rm = TRUE), 1), 
           start       = pmin(start, force_start, na.rm = TRUE))
  
  # Day last vaccination given for each group
  group_end = v$groups[, ifelse(max_end > 0, max_end, init_end)]
  
  # Check that all groups are finished vaccinating before boosters are forced
  force_check = group_end >= v$booster_groups$force_start
  force_fail  = v$groups$id[sapply(force_check, isTRUE)]
  
  # Sanity check that boosters are not forced before everyone in this group is vaccinated
  if (length(force_fail) > 0)
    stop("Vaccine group(s) are forced boosters before everyone is even vaccinated: ",
         paste(force_fail, collapse = ", "))
  
  # Preallocate temporal vaccine coverage matrix (one column per priority group)
  v$coverage = t(matrix(data = v$groups$init_coverage, 
                        ncol = y$n_days, 
                        nrow = v$n_groups, 
                        dimnames = list(v$groups$id)))
  
  # Loop through the priority groups
  for (i in seq_len(v$n_groups)) {
    this_group = v$groups[i, ]
    
    # Skip this process if no max coverage to reach
    if (this_group$max_start < y$n_days && 
        this_group$max_end - this_group$max_start > 0) {
      
      # Linear scale up from current to maximal coverage
      growth_vec = seq(this_group$init_coverage, this_group$max_coverage, 
                       length.out = this_group$max_end - this_group$max_start + 1)
      
      # Date indicies of growth scale up - may be truncated by n_days
      growth_idx = seq_along(growth_vec) + this_group$max_start
      growth_idx = growth_idx[growth_idx <= y$n_days]
      
      # Apply this coverage - clipping vector if necessary
      v$coverage[growth_idx, i] = growth_vec[-1][seq_along(growth_idx)]
      
      # Set remaining days to max coverage
      if (this_group$max_end < y$n_days)
        v$coverage[max(growth_idx) : y$n_days, i] = this_group$max_coverage
    }
  }
  
  # Sanity check that vaccine coverage does not contain NAs
  if (any(is.na(v$coverage)))
    stop("Vaccine coverage matrix contains NAs")
  
  # Remove redundant items
  y[c("vaccine_history", "vaccine_rollout")] = NULL
  
  # We can now append the vaccine list to output list
  y$vaccine = v
  
  # ---- Model states and model metrics ----
  
  # Disease, care, and prognosis states
  y = parse_model_states(o, y)
  
  # Model metrics to be tracked
  y = parse_model_metrics(o, y)
  
  # ---- Model output groupings ----
  
  # Sanity check that max number of infections allowed is positive
  if (y$max_infections < 1)
    stop("Parameter 'max_infections' must be a positive integer")
  
  # List of all possible metric groups, and their respective groups
  y$count = list(age           = y$age$all, 
                 infections    = 1 : y$max_infections, 
                 variant       = y$variants$id, 
                 vaccine_group = c(v$groups$id, "none"))
  
  # Dictionaries for variants and vaccine groups
  variant_dict = setNames(y$variants$name, y$variants$id)
  vaccine_dict = setNames(v$groups$name,   v$groups$id) %>% 
    c(setNames("Ineligible", "none"))
  
  # Append to dict list (current contains metric names and IDs)
  y$dict = list.append(y$dict, variant = variant_dict, vaccine_group = vaccine_dict)
  
  # Remove redundant items
  y[c("max_infection_count")] = NULL
  
  # ---- Final formatting ----
  
  # Order elements alphabetically to make it easy to search
  y = y[order(names(y))]
  
  # Provide both the raw yaml file and the parsed parameters used in the model
  yaml = list(raw = as.yaml(y_user), parsed = y)
  
  return(yaml)
}

# ---------------------------------------------------------
# Overwrite any parameter defaults for which we have user-defined values 
# ---------------------------------------------------------
overwrite_defaults = function(y, y_overwrite) {
  
  # Sanity checks on user-defined inputs
  do_checks(y1 = y, y2 = y_overwrite)
  
  # ---- Special case: named vaccines ----
  
  # Error message to help point user in the right direction
  use_custom = "use 'custom' to apply alternative properties"
  
  # Check a vaccine been defined (if not, a default is used)
  if (!is.null(y_overwrite$vaccine_type)) {
    
    # Allow only named vaccines or 'custom'
    if (!y_overwrite$vaccine_type %in% c("custom", names(y$vaccine_defintion)))
      stop("Vaccine '", y_overwrite$vaccine, "' is not defined - ", use_custom)
    
    # If 'custom' vaccine has not been selected, it doesn't make sence to define custom properties
    if (y_overwrite$vaccine_type != "custom" && !is.null(y_overwrite$vaccine_custom))
      stop("It is an error to select a named vaccine but also define custom vaccine properties")
  }
  
  # Do not allow the overwriting of vaccine definitions period - use 'custom' for this
  if (!is.null(y_overwrite$vaccine_defintion))
    stop("It is an error to redefine named vaccines - ", use_custom)
  
  # NOTE: If you really do want to redefine named vaccines, do this in default.yaml
  
  # ---- Special case: booster force start ----
  
  # Booster start date uses FALSE to denote off - convert to integer NA for class consistency
  na_fn = function(x) rapply(x, function(w) ifelse(isFALSE(w), as.integer(NA), w), how = "replace")
  
  # Apply to default yaml
  y$booster_rollout = na_fn(y$booster_rollout)
  
  # ... and also to user-defined yaml (if defined)
  if (!is.null(y_overwrite$booster_rollout))
    y_overwrite$booster_rollout = na_fn(y_overwrite$booster_rollout)
  
  # ---- Special case: demography ----
  
  # If defined by user, overwrite
  if (!is.null(y_overwrite$demography))
    y$demography = y_overwrite$demography
  
  # Then remove from user-defined list
  y_overwrite$demography = NULL
  
  # ---- Special case: scenarios ----
  
  # The one exception for overwriting is scenarios, which we'll concatenate
  y$scenarios = c(y$scenarios, y_overwrite$scenarios)
  
  # User defined scenarios can now be removed from this list
  y_overwrite$scenarios = NULL
  
  # ---- Name all unnamed lists with unique IDs ----
  
  # Start by naming scenarios with their unique IDs
  scenario_ids = unlist(lapply(y$scenarios, function(x) x$id))
  y$scenarios  = setNames(y$scenarios, scenario_ids)
  
  # Then for the lists within each scenario (as scenarios can have nested unnamed lists)
  for (scenario in names(y$scenarios))
    y$scenarios[[scenario]] = name_lists(y$scenarios[[scenario]])
  
  # Then do the same for all other unnames lists
  y = name_lists(y, y_overwrite = y_overwrite)
  
  # Also for what we still need to overwrite with
  y_overwrite = name_lists(y_overwrite)
  
  # ---- Overwrite defaults ----
  
  # All individual items to overwrite
  overwrite_items = names(unlist_format(y_overwrite))
  
  # Which of the these are functions - these are special cases
  fn_call_idx = grepl("\\$fn$", overwrite_items)
  fn_items = str_remove(overwrite_items[fn_call_idx], "\\$fn")
  
  # Skip this if no function items
  if (length(fn_items) > 0) {
    
    # Apply function using eval - this is the only way I can think of for n-nested lists
    #
    # NOTE: No checks needed as key-pairs can differ depending on function defined
    for (fn in fn_items)
      eval(parse(text = paste0("y$", fn, " = y_overwrite$", fn)))
    
    # All items that pertain to a function call
    fn_item_idx = paste0("^") %>%
      paste0(fn_items, collapse = "|") %>%
      str_replace_all("\\$", "\\\\$") %>%
      grepl(overwrite_items)
    
    # Drop these function items from recursive overwriting
    overwrite_items = overwrite_items[!fn_item_idx]
  }
  
  # Loop through whats left: values to overwrite with
  for (overwrite_item in overwrite_items) {
    
    # The value to overwrite with, and the original (to check for consistency)
    overwrite_value = eval(parse(text = paste0("y_overwrite$", overwrite_item)))
    default_value   = eval(parse(text = paste0("y$", overwrite_item)))
    
    # Normally we will be replacing some non-trivial value
    if (length(default_value) > 0) {
      
      # Check class consistency
      do_checks(y1 = setNames(default_value, overwrite_item), 
                y2 = setNames(overwrite_value, overwrite_item))
      
    } else {  # Otherwise value has been trivialised within name_lists
      
      # However class should still be consistent
      if (class(overwrite_value) != class(default_value))
        stop("Inconsistent data class in input yaml file: \n", 
             paste0(" ! ", overwrite_item, ": ", class(default_value), 
                    " -> ", class(overwrite_value), "\n"))
    }
    
    # If a string, wrap in quotes so eval knows it's not a variable
    if (is.character(overwrite_value))
      overwrite_value = paste0("'", overwrite_value, "'")
    
    # Ensure class is preserved by using as.xxx
    as_class = paste0("get('as.", class(default_value), "')(", overwrite_value, ")")
    
    # Apply the assignment using eval - this is the only way I can think of for n-nested lists
    eval(parse(text = paste0("y$", overwrite_item, " = ", as_class)))
  }
  
  return(y) 
}

# ---------------------------------------------------------
# Sanity checks on user-defined inputs
# ---------------------------------------------------------
do_checks = function(y1 = NULL, y2 = NULL) {
  
  # Names of paramters defined by user
  user_params = names(y2)
  
  # Check if any parameters have been repeated
  duplicate_items = user_params[duplicated(user_params)]
  
  # Throw an error if this is the case
  if (length(duplicate_items) > 0)
    stop(" ! Duplicated items in input yaml file: ", paste0(duplicate_items, collapse = ", "))
  
  # Check if any unknown items been defined
  unknown_items = setdiff(user_params, names(y1))
  
  # Throw an error if this is the case
  if (length(unknown_items) > 0)
    stop(" ! Unrecognised items in input yaml file: ", paste0(unknown_items, collapse = ", "))
  
  # Data class of all user-defined parameters
  class_user = unlist(lapply(y2, class))
  
  # Default data class of these same parameters
  class_def = unlist(lapply(y1, class))[user_params]
  
  # Parameters that will cause problems as they do not have the same data class
  errors = names(class_def)[!compare.list(class_def, class_user)]
  
  # This isn't a problem if these are 'function lists'
  for (err_param in errors) {
    
    # Check if either class is a list
    is_list = c(class_def[err_param], class_user[err_param]) == "list"
    
    # If yes, continue
    if (any(is_list)) {
      
      # Extract the elements from the list of interest
      list_items = get(paste0("y", which(is_list)))[[err_param]]
      
      # If this is a 'function list', remove this parameter from vector of errors
      if (names(list_items)[[1]] == "fn")
        errors = setdiff(errors, err_param)
    }
  }
  
  # Any inconsistencies remain?
  if (length(errors) > 0) {
    
    # Throw an error reporting which parameters, what it should be, and what it is
    stop("Inconsistent data class in input yaml file: \n", 
         paste0(" ! ", errors, ": ", class_def[errors], " -> ", class_user[errors], "\n"))
  }
}

# ---------------------------------------------------------
# Apply names to unanmed lists with unique IDs
# ---------------------------------------------------------
name_lists = function(y, y_overwrite = NULL) {
  
  # Items that are themselves lists
  y_lists = y[unlist(lapply(y, is.list))]
  
  # Of those lists, the ones which have NULL names (ie unnamed lists)
  y_null = y_lists[unlist(lapply(lapply(y_lists, names), is.null))]
  
  # Of those unnamed lists, the ones which have unique IDs
  y_id = y_null[unlist(lapply(y_null, function(x) names(x[[1]])[1] == "id"))]
  
  # Apply the IDs as the name of the previously unnamed list
  y_named = lapply(y_id, function(x) setNames(x, list2dt(x, fill = TRUE)$id))
  
  # Which of this now named lists will be overwritten
  for (param in intersect(names(y_overwrite), names(y_named))) {
    
    # IDs of this parameters before and after the overwrite
    param_ids_default = unlist(lapply(y[[param]], function(x) x$id))
    param_ids = unlist(lapply(y_overwrite[[param]], function(x) x$id))
    
    # We'll want to reset all values if not identical
    if (!identical(param_ids, param_ids_default)) {
      
      # Store the basic structure of each item
      param_format = y_named[[param]][1]
      
      # Loop through elements are extract class
      for (param_el in names(param_format[[1]])) {
        param_class = class(param_format[[1]][[param_el]])
        
        # Reset value by trivialising, but retraining class
        eval(parse(text = paste0("param_format[[1]]$", param_el, " = ", param_class, "(0)")))
      }
      
      # Set up the necessary number of items and apply the basic structure
      y_named[[param]] = NULL
      y_named[[param]][param_ids] = param_format
    }
  }
  
  # Store these newly named lists
  y[names(y_named)] = y_named
  
  return(y)
}

# ---------------------------------------------------------
# Construct function call as a string to be evaluated
# ---------------------------------------------------------
parse_fn = function(fn_args, along = NULL, evaluate = TRUE) {
  
  # Extract function name from input list
  fn = fn_args$fn
  
  # Remove this function field to leave only function arguments
  fn_args$fn = NULL
  
  # We may want to append the values to evaluate the function 'along'
  if (!is.null(along))
    fn_args = append(along, fn_args)
  
  # Collapse all key-value arguments into a comma-seperated string
  args_string = paste0(names(fn_args), " = ", fn_args, collapse = ", ")
  
  # Concatenate function name with arguments
  fn_call = paste0(fn, "(", args_string, ")")
  
  # Either return that string evaluated...
  if (evaluate == TRUE) {
    fn_eval = eval(parse(text = fn_call))
    
    return(fn_eval)
    
    # ... or simply the function string
  } else {
    return(fn_call)
  }
}

# ---------------------------------------------------------
# Parse user-defined scenarios - in a seperate function for readability
# ---------------------------------------------------------
parse_scenarios = function(o, y, scenario, read_array) {
  
  # Append scenario ID to main y list, simply for reference
  y$.id = scenario
  
  # Trivial process for the baseline scenario - just store 'name'
  if (scenario == "baseline") {
    y$.name = y$scenarios$baseline$name
    
  } else {  # Otherwise, work to do...
    
    # First flatten any array scenarios
    if (read_array == TRUE)
      y = parse_arrays(o, y, scenario)
    
    # Extract short and long names of all scenarios we've defined
    scenarios_id   = names(y$scenarios)
    scenarios_name = unlist(lapply(y$scenarios, function(x) x$name))
    
    # ---- Sanity checks on scenario definitions ----
    
    # Any scenarios for which we do not have a unique 'name' item
    unset_names = scenarios_id[!scenarios_id %in% names(scenarios_name)]
    duplicate_names = scenarios_id[duplicated(scenarios_name)]
    
    # Throw an error if any names are missing
    if (length(unset_names) > 0)
      stop("! All scenarios must have a 'name' item: ", paste0(unset_names, collapse = ", "))
    
    # Throw an error if any names are duplicated
    if (length(duplicate_names) > 0)
      stop("! All scenarios must have a unique 'name': ", paste0(duplicate_names, collapse = ", "))
    
    # If we only want to read scenario names, return out here
    if (scenario == "*read*")
      return(list(scenario_names = scenarios_name))
    
    # Check the selected scenario actually exists
    if (!scenario %in% scenarios_id) {
      
      # If it's a parent of an array scenario, provide a useful error message
      # if (scenario %in% get_array_parents(names(y$scenarios)))
      #   stop("Scenario '", scenario, "' is an array parent - you must define a child to simulate")
      
      # Otherwise just report that such a scenario isn't defined
      stop("Scenario '", scenario, "' not recognised")
    }
    
    # ---- Apply scenario values ----
    
    # Trivial process for the baseline scenario
    if (scenario != "baseline") {
      
      # Index of the alternative scenario we want to apply
      scenario_idx = which(scenarios_id == scenario)
      
      # Remove id and name items to leave only parameters to overwrite with
      y_scenario = list.remove(y$scenarios[[scenario_idx]], c("id", "name"))
      
      # Use overwrite function to check and apply item values
      y = overwrite_defaults(y, y_scenario)
    }
    
    # Keep reference to this scenario 'name'
    y$.name = scenarios_name[[scenario]]
  }
  
  # Can now safely remove scenarios field
  y$scenarios = NULL
  
  return(y)
}

# ---------------------------------------------------------
# Flatten and parse parent array scenarios into their children
# ---------------------------------------------------------
parse_arrays = function(o, y, scenario) {
  
  # Full names of all scenarios
  scenarios_name = unlist(lapply(y$scenarios, function(x) x$name))
  
  # ---- Scenario arrays ----
  
  # Indentify where (if anywhere) arrays are defined
  array_idx   = grepl("*\\.array\\.*", names(unlist(y$scenarios)))
  array_items = names(unlist(y$scenarios))[array_idx]
  
  # Names of parent scenarios that have array items
  parents = get_array_parents(array_items)
  
  # Loop through all parent array scenarios
  for (parent in parents) {
    
    # Array items associated with this parent
    array_scen_idx   = grepl(paste0("^", parent, "\\."), array_items)
    array_scen_items = str_remove(array_items[array_scen_idx], paste0(parent, "\\."))
    
    # Strip away to leave unique variables we have arrays for
    array_vars = unique(str_remove(array_scen_items, "\\.array\\..*"))
    array_vars = str_replace_all(array_vars, "\\.", "$")
    
    # Preallocate lists to store array IDs, names, and values
    array_vals = array_ids = array_names = list()
    
    # Preallocate lists to store array variable IDs and names
    var_ids = var_names = character()
    
    # Loop through variables that have arrays
    for (array_var in array_vars) {
      
      # String which - when evaluated - indexes the deep nested array details for this variable 
      array_str = paste("y$scenarios", parent, array_var, "array", sep = "$")
      
      # Evaluate the string to obtain a list of inputs to R's seq function
      array_eval = eval(parse(text = array_str))
      
      # Remove the ID and name items, and set the function name to seq
      array_fn = c(fn = "seq", list.remove(array_eval, c("id", "name")))
      
      # Evaluate the seq function
      all_vals = parse_fn(array_fn)
      
      # Use child scenario indices for IDs
      all_ids = paste0(array_eval$id, 1 : length(all_vals))
      
      # Store these IDs and values in seperate lists
      array_ids[[array_var]]  = all_ids
      array_vals[[array_var]] = all_vals
      
      # Also store child scenario names as a named vector
      array_names[[array_var]] = setNames(paste0(array_eval$name, ": ", all_vals), all_ids)
      
      # Also store full name of variable as defined by user
      var_ids   = c(var_ids,   array_eval$id)
      var_names = c(var_names, array_eval$name)
    }
    
    # Now take the full factorial of IDs and associated names
    array_scen_ids   = unite(expand.grid(array_ids),   "x", sep = ".")$x
    array_scen_names = unite(expand.grid(array_names), "x", sep = ", ")$x
    
    # Do the same for the values, and concatenate it all together
    array_scen_df = expand.grid(array_vals) %>%
      mutate(id   = paste0(parent, ".", array_scen_ids),
             name = paste0(scenarios_name[[parent]], 
                           " (", array_scen_names, ")"), 
             parent = parent)
    
    # ---- Construct child scenarios ---- 
    
    # Store details of the parent array scenario
    array_details = y$scenarios[[parent]]
    
    # Loop through each individual child
    for (i in 1 : nrow(array_scen_df)) {
      this_scen = array_scen_df[i, ]
      
      # Start with the basic structure of the parent
      scen_details = list_modify(array_details, 
                                 id   = this_scen$id, 
                                 name = this_scen$name)
      
      # Reduce array_df down to just values we want to use for this child scenario
      val_df = select(this_scen, -id, -name, -parent)
      
      # Loop through the values
      for (val_name in names(val_df)) {
        
        # String to be evaluated to set the individual value
        val_str = paste0("scen_details$", val_name, " = ", val_df[[val_name]])
        
        # Evaluate the string to set the value
        eval(parse(text = val_str))
      }
      
      # Append this child scenario
      y$scenarios[[this_scen$id]] = scen_details
    }
    
    # Remove the parent scenario
    y$scenarios[[parent]] = NULL
    
    # ---- Store array details ----
    
    # Only store array details if not reading scenario names
    if (scenario != "*read*") {
      
      # Initiate array info list with basic scenario ID dataframe
      array_info = list(meta = select(array_scen_df, parent, scenario = id, name))
      
      # Append info about the variables - IDs and names
      array_info$vars = data.table(variable_id   = var_ids, 
                                   variable_name = var_names, 
                                   variable_call = array_vars)
      
      # Store all the values to be simulated
      array_info$values = array_scen_df %>%
        select(scenario = id, all_of(array_vars)) %>%
        rename(setNames(array_vars, var_ids)) %>%
        pivot_longer(cols = -scenario, names_to = "variable_id") %>%
        as.data.table()
      
      # Save this info for use when collating all results
      saveRDS(array_info, file = paste0(o$pth$arrays, parent, ".rds"))
    }
  }
  
  return(y)
}

# ---------------------------------------------------------
# Parse disease, care, and prognosis states
# ---------------------------------------------------------
parse_model_states = function(o, y) {
  
  # Load excel file
  y$model_flows = read_excel(o$pth$states, 1) %>%
    select(-notes) %>%
    as.data.table()
  
  # Distinct disease, care, and prognosis states
  y$model_states$disease   = unique(na.omit(y$model_flows$disease))
  y$model_states$care      = unique(na.omit(y$model_flows$care))
  y$model_states$prognosis = unique(na.omit(y$model_flows$prognosis))
  
  # Prognosis states relating to severe disease - used when considering variant severity factor
  all_prognosis    = setdiff(y$model_states$prognosis, "none")
  severe_prognosis = setdiff(all_prognosis, c("asym", "mild"))
  
  # Store the indicies of severe and non-severe disease for indexing prognosis matrix
  y$prognosis_idx = list(severe     = all_prognosis %in% severe_prognosis, 
                         non_severe = !all_prognosis %in% severe_prognosis)
  
  # Once diagnosed, only those with certain prognoses will go into isolation
  y$model_states$iso = unique(filter(y$model_flows, care_state == "iso")$prognosis_state)
  
  # Concatenate all states for easy referencing - used to plot states over time in model()
  y$model_states$all = c(y$model_states$disease, y$model_states$care)
  
  return(y)
}

# ---------------------------------------------------------
# Parse model metrics to collate and plot
# ---------------------------------------------------------
parse_model_metrics = function(o, y) {
  
  # Apply formating, convert to dt, append descriptions, and filter
  metrics_df = rbindlist(y$model_metrics, fill = TRUE) %>% 
    mutate(metric = names(y$model_metrics), .before = 1) %>% 
    rename(grouping = by) %>% 
    right_join(read_excel(o$pth$metrics, 1), by = "metric") %>%
    mutate(aggregate  = as.logical(aggregate), 
           temporal   = as.logical(temporal), 
           cumulative = as.logical(cumulative), 
           scaled     = as.logical(scaled), 
           grouping   = ifelse(metric == "n_infections", "infections", grouping),  # Special case
           grouping   = ifelse(metric == "variant_prevalence", "variant", grouping),  # Special case
           grouping   = ifelse(is.na(grouping), "na", grouping)) %>%
    filter(report == TRUE) %>%
    select(-report, -notes)
  
  # All the various types of groupings requested by the user
  groupings = str_split(str_remove_all(metrics_df$grouping, " "), ",")
  
  # Metrics with grouping (repeat elements if we have multiple groupings)
  grouping_metrics = rep(metrics_df$metric, unlist(lapply(groupings, length)))
  
  # Convert to named vector - used in update_output in model.R
  all_groupings = setNames(unlist(groupings), grouping_metrics)
  
  # Groupings known and understood by the model
  known_groupings = c("age", "variant", "vaccine_group", "infections", "none", "na")
  
  # Are any unknown
  unknown_groupings = setdiff(unique(all_groupings), known_groupings)
  if (length(unknown_groupings > 0))
    stop("Unknown metric grouping(s): ", paste(unknown_groupings, collapse = ", "))
  
  # ---- Parse the input ----
  
  # Store the key details - metric groupings and whether temporal and/or cumulative
  y$metrics$df = select(metrics_df, -description, -cumulative_description)
  y$metrics$groupings = all_groupings
  
  # Store metric dictionary ('description' column)
  y$dict$metric = setNames(metrics_df$description, metrics_df$metric)
  
  # Also store cumulative metric dictionary ('cumulative_description' column)
  cumulative_df = filter(metrics_df, cumulative == TRUE)
  y$dict$cumulative = setNames(cumulative_df$cumulative_description, cumulative_df$metric)
  
  # Remove redundant items
  y[c("model_metrics")] = NULL
  
  return(y)
}

# ---------------------------------------------------------
# Extract parent scenario names from a vector of child names
# ---------------------------------------------------------
get_array_parents = function(scenarios, multidim = FALSE) {
  
  if (!is.null(names(scenarios)))
    stop("Should this be allowed?")
  
  # Array scenario IDs contain periods
  parents_df = scenarios %>%
    str_split_fixed( "\\.", n = 2) %>%
    as.data.table() %>%
    setnames(qc(parent, sub)) %>%
    filter(sub != "")
  
  # Extract unique parents
  parents = unique(parents_df$parent)
  
  # We may also want to determine if arrays are multi-dimensional  
  if (multidim == TRUE) {
    
    # Loop through array parents
    for (this_parent in parents) {
      children = subset(parents_df, parent == this_parent)$sub
      
      # Multi-dimensional arrays identified by period symbols in child names
      is_multidim = all(grepl("\\.", children))
      
      # Remove this parent if not a multi-dimensional array
      if (!is_multidim)
        parents = setdiff(parents, this_parent)
    }
  }
  
  return(parents)
}

# ---------------------------------------------------------
# Extract children scenario names given a parent
# ---------------------------------------------------------
get_array_children = function(scenarios, parents) {
  
  # Initiate child character vector
  children = character()
  
  # Iteraate through any parents provided (often only one)
  for (parent in parents) {
    
    # Expression to match to identify children
    parent_exp = paste0("^", parent, "\\.*")
    
    # Indices of such scenarios
    children_idx = grepl(parent_exp, scenarios)
    
    # Extract scenario names
    children = c(children, scenarios[children_idx])
  }
  
  return(children)
}

