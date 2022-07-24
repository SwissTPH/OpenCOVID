###########################################################
# LOAD INPUT
#
# Load and format input files for use in various stages of 
# the model pipeline.
#
###########################################################

# ---------------------------------------------------------
# Load and format disease, care, and prognosis states
# ---------------------------------------------------------
format_model_states = function(o) {
  
  # Load excel file
  o$model_flows = read_excel(o$pth$states, 1) %>%
    select(-notes) %>%
    as.data.table()
  
  # Distinct disease, care, and prognosis states
  o$disease_states   = unique(na.omit(o$model_flows$disease))
  o$care_states      = unique(na.omit(o$model_flows$care))
  o$prognosis_states = unique(na.omit(o$model_flows$prognosis))
  
  # Concatenate all states for easy referencing - used to plot states over time in model()
  o$all_states = c(o$disease_states, o$care_states)  # Also prognosis_states? Depends if want to plot them
  
  return(o)
}

# ---------------------------------------------------------
# Load calibration csv and extend for canton-specific parameters
# ---------------------------------------------------------
format_calibration_parameters = function(o, test_param) {
  
  # Load calibration details from file
  o$parameter_table = read_excel(o$pth$calib_file, 1) %>%
    mutate(fixed = 1 - fitted) %>%
    select(names, values = prior_mean, fixed, global, 
           prior_sd, prior_weight, lower_bound, upper_bound) %>% 
    as.data.frame()
  
  # ---- Force parameters off for alternative functionality ----
  
  # Alternative functionality for testing a single paramter (see unit_tests.R)
  if (!is.null(test_param)) {
    
    # Index of parameter we want to test
    test_idx = o$parameter_table$names == test_param
    
    # Fix all parameters aside from the one we want to test
    o$parameter_table$fixed[] = 1
    o$parameter_table$fixed[test_idx] = 0
  }
  
  # ---- Override canton names if necessary ----
  
  # Check for my_options file - override cantons option only 
  o = override_options(o, subset = "cantons", quiet = TRUE)
  
  # Overwrite if we want to run all cantons
  if ("all" %in% tolower(o$cantons)) {
    
    # Check that other cantons have not been defined - that'd just be confusing
    if (length(o$cantons) > 1)
      stop("If running 'all' cantons you should not also select individual cantons - see o$cantons")
    
    # Overwrite the value of o$cantons
    o$cantons = o$all_cantons
  }
  
  # Logical vector of valid canton codes defined by user
  valid_cantons = o$cantons %in% o$all_switzerland
  
  # Throw an error if any of these are not valid
  if (!all(valid_cantons))
    stop("Canton codes not recognised: ", paste0(o$cantons[!valid_cantons], collapse = ", "))
  
  # ---- Seperate fixed and calibrated (global and local) parameters ----
  
  # Seperate calibration parameters from fixed parameters
  param_df = o$parameter_table %>% filter(fixed == 0)
  fixed_df = o$parameter_table %>% filter(fixed == 1) %>% 
    select(names, values, lower_bound, upper_bound, prior_sd)
  
  # Throw an error if all parameters are fixed - this doesn't seem like a valid use case
  if (nrow(param_df) == 0)
    stop("Zero parameters are being calibrated (see ", basename(o$pth$calib_file), ")")
  
  # Seperate global parameters and define scope
  global_df = param_df %>% filter(global == 1) %>% mutate(scope = "global")
  
  # Dataframe for canton-specifc params
  local_df = param_df %>% filter(global == 0)
  
  # If no non-global parameters, define trivial 'scope' variable
  if (nrow(local_df) == 0) {
    
    # A slightly hacky way of ensuring we maintain the correct column order
    var_names = setdiff(names(local_df), "names")
    local_df  = local_df %>% mutate(scope = NA) %>% select(names, scope, one_of(var_names))
    
  } else { # Otherwise, define the scope of each local parameter
    
    # Full factorial expansion for canton-specific parameters
    local_df = base::merge(data.frame(scope = o$cantons), local_df)
    local_df = unite(local_df, "names", "names", "scope", sep = ".", remove = FALSE)
  }
  
  # ---- Prepare output ----
  
  # Concatenate these dataframes together
  calibration_df = rbind(global_df[names(local_df)], local_df)
  
  # Using data tables to improve indexing speed
  o$fixed_df       = as.data.table(fixed_df)
  o$calibration_df = as.data.table(calibration_df)
  
  return(o)
}

# ---------------------------------------------------------
# Load calibration multiplers csv and apply defaults where necessary
# ---------------------------------------------------------
format_calibration_multipliers = function(o) {
  
  # Increase importance of fitting different metrics 
  calib_mult = read.csv(o$pth$calib_mult)  # File with canton-specific metric multipliers
  
  # ---- Sanity checks on inputs provided ----
  
  # Reference to file and dir so we can direct the user there
  file_dir = substring(str_remove(o$pth$calib_mult, o$pth$code), 2)
  file_ref = paste0(" - see '", file_dir, "'")
  
  # Ensure a canton is referenced at the most once
  if (any(duplicated(calib_mult$canton)))
    stop("Canton references must be unique", file_ref)
  
  # ---- Use default weightings unless otherwise specified ----
  
  # INTERPRETATION: 1 means normal weight, > 1 means more weight, < 1 means less weight
  
  # Expand out to consider all cantons being modelled and set default values were necessary
  calib_mult_all = data.frame(canton = o$cantons) %>% 
    left_join(calib_mult, by = "canton")
  
  # Inidices of cantons that are defined in this table
  defined_idx = calib_mult_all$canton %in% calib_mult$canton
  
  # Index of default settings
  default_idx = calib_mult$canton == "default"
  
  # Apply default settings to any canton not defined
  calib_mult_all[!defined_idx, -1] = calib_mult[default_idx, -1]
  
  # Store this table in o
  o$calibration_multiplier = calib_mult_all
  
  return(o)
}

# ---------------------------------------------------------
# Load and format SARS-CoV-2 variant properties
# ---------------------------------------------------------
format_variant_properties = function(o) {
  
  # Load variant properties file
  variants_file = read.csv(o$pth$variants) %>%
    mutate(primary    = as.logical(primary), 
           calibrated = as.logical(calibrated))
  
  # ---- Seperate out primary variant ----
  
  # Details of primary variant (dominant variant in first wave)
  variant_primary = variants_file %>%
    filter(primary == TRUE)
  
  # Details of later variants
  variant_mutation = variants_file %>%
    filter(primary == FALSE) %>%
    mutate(import_date   = format_date(import_date), 
           import_number = as.numeric(import_number), 
           infectivity_factor = as.numeric(infectivity_factor), 
           severity_factor    = as.numeric(severity_factor)) %>%
    select(-primary)
  
  # Vector of names for all variants - primary first
  #
  # NOTE: Reconstructing name vector rather than using variants_file$name to ensure primary first
  v_names = c(variant_primary$name, variant_mutation$name)
  
  # Factors of all variants - including trivial factor for primary variant
  infectivity = c(1, variant_mutation$infectivity_factor)
  severity    = c(1, variant_mutation$severity_factor)
  
  # Logical vector of which variant is flagged for potential calibration
  calibrated = c(FALSE, variant_mutation$calibrated)
  
  # ---- Sanity checks on inputs provided ----
  
  # Reference to file and dir so we can direct the user there
  file_dir = substring(str_remove(o$pth$variants, o$pth$code), 2)
  file_ref = paste0(" - see '", file_dir, "'")
  
  # We must have one variant selected as the primary variant
  if (nrow(variant_primary) == 0)
    stop("No variant has been selected as 'primary'", file_ref)
  
  # ... but no more than one
  if (nrow(variant_primary) > 1)
    stop("Multiple variants selected as 'primary'", file_ref)
  
  # Ensure vectors are numeric
  if (any(is.na(infectivity))) stop("Infectivity values not numeric", file_ref)
  if (any(is.na(severity)))    stop("Severity values not numeric", file_ref)
  
  # We want one variant to be identified for potential calibration, but it should not be the primary
  if (variant_primary$calibrated == TRUE)
    stop("It is an error to attempt to calibrate the 'primary' reference variant", file_ref)
  
  # ... it must be one of the later mutations
  if (sum(variant_mutation$calibrated) == 0)
    stop("Identify one non-primary variant for potential calibration", file_ref)
  
  # ... but we are not currently supporting the calibration of multiple variants at once
  if (sum(variant_mutation$calibrated) > 1)
    stop("Calibrating multiple variants at once is currently not supported", file_ref)
  
  # ---- Further formatting ----
  
  # Initiate a list for variant information
  o$variants = list(primary  = variant_primary$name, 
                    mutation = variant_mutation$name, 
                    import_date   = variant_mutation$import_date, 
                    import_number = variant_mutation$import_number,
                    calibrated    = calibrated, 
                    infectivity_factor = setNames(infectivity, v_names), 
                    severity_factor    = setNames(severity,    v_names))
  
  # Construct a named vector of AKAs to operate as a dictionary
  o$variants_dict = setNames(c(variant_primary$aka,  variant_mutation$aka), 
                             c(variant_primary$name, variant_mutation$name))
  
  # Easy access variant names and the number of variants (including primary)
  o$variant_names = names(o$variants_dict)
  o$n_variants    = length(o$variants_dict)
  
  # Assocaite each variant to a number - this will be used in model.R for speed
  o$variants_idx = setNames(1 : o$n_variants,  names(o$variants_dict))
  
  # TODO: Throw an error if non-unique dates for variant import dates, it adds a lot of avoidable
  #       code to deal with the edge cases.
  
  return(o)
}

# ---------------------------------------------------------
# Load and format vaccine properties
# ---------------------------------------------------------
format_vaccine_properties = function(o) {
  
  # Load vaccine properties file
  vaccines = read.csv(o$pth$vaccines)
  
  # Construct a named vector of AKAs to operate as a dictionary
  o$vaccine_dict = setNames(vaccines$aka, vaccines$name)
  
  # Provide vaccine efficacy as a list
  vaccine_efficacy   = split(vaccines$efficacy, seq(nrow(vaccines)))
  o$vaccine_efficacy = setNames(vaccine_efficacy, vaccines$name)
  
  # Number of doses as a named vector
  o$n_doses = setNames(vaccines$n_doses, vaccines$name)
  
  # Name of the vaccine top be used as a default
  o$vaccine_default = vaccines$name[vaccines$default]
  
  # Assocaite each vacine priority group to a number - this will be used in model.R for speed
  o$vaccine_priority_idx = setNames(1 : length(o$vaccine_priority), 
                                    names(o$vaccine_priority))
  
  # ---- Sanity checks on inputs provided ----
  
  # Reference to file and dir so we can direct the user there
  file_dir = substring(str_remove(o$pth$vaccines, o$pth$code), 2)
  file_ref = paste0(" - see '", file_dir, "'")
  
  # Check that exactly one vaccine is defined as the default
  if (length(o$vaccine_default) != 1)
    stop("One (and only one) vaccine should be defined as the default", file_ref)
  
  return(o)
}

# ---------------------------------------------------------
# Load model metrics to collate and plot
# ---------------------------------------------------------
format_model_metrics = function(o) {
  
  # Load excel file and remove metrics we're not modelling
  metrics_df = read_excel(o$pth$metrics, 1) %>%
    mutate(capture    = as.logical(capture), 
           aggregate  = as.logical(aggregate), 
           temporal   = as.logical(temporal), 
           cumulative = as.logical(cumulative), 
           scaled     = as.logical(scaled), 
           grouping   = ifelse(tolower(grouping) == "na", "na", grouping), 
           calibrate_grouping = ifelse(tolower(calibrate_grouping) == "na", "na", calibrate_grouping), 
           delay      = ifelse(tolower(delay) == "none", NA, delay)) %>%
    filter(capture == TRUE) %>%
    select(-capture, -notes) %>%
    as.data.table()
  
  # All the various types of groupings requested by the user
  groupings = str_split(str_remove_all(metrics_df$grouping, " "), ",")

  # Metrics with grouping (repeat elements if we have multiple groupings)
  grouping_metrics = rep(metrics_df$metric, unlist(lapply(groupings, length)))

  # Convert to named vector - used in update_output in model.R
  all_groupings = setNames(unlist(groupings), grouping_metrics)
  
  # ---- Sanity checks on inputs provided ----
  
  # Reference to file and dir so we can direct the user there
  file_dir = substring(str_remove(o$pth$metrics, o$pth$code), 2)
  file_ref = paste0(" - see '", file_dir, "'")
  
  # Grouping known and understood by the model
  known_groupings = c("age", "variant", "vaccine_priority", "infections", "none", "na")
  
  # Are any unknown
  unknown_groupings = setdiff(unique(all_groupings), known_groupings)
  if (length(unknown_groupings > 0))
    stop("Unknown grouping(s): ", paste(unknown_groupings, collapse = ", "), file_ref)
  
  # ---- Parse the input ----
  
  # Store the key details - metric groupings and whether temporal and/or cumulative
  o$all_metrics = select(metrics_df, -description, -cumulative_description)
  o$metric_groupings = all_groupings
  
  # Store metric dictionary ('description' column)
  o$metric_dict = setNames(metrics_df$description, metrics_df$metric)
  
  # Also store cumulative metric dictionary ('cumulative_description' column)
  cumulative_df = filter(metrics_df, cumulative == TRUE)
  o$cumulative_dict = setNames(cumulative_df$cumulative_description, cumulative_df$metric)
  
  return(o)
}

