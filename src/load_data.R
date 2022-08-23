###########################################################
# LOAD DATA
#
# Load epi data that may - or may not - be required for
# calibration. Also able to generate synthetic data to test
# fitting algorithm.
#
###########################################################

# ---------------------------------------------------------
# Load emperical data from source (or generate synthetic data)
# ---------------------------------------------------------
load_data = function(o, p, synthetic = NULL) {
  
  # Shorthand for calibration options
  opts = p$calibration_options
  
  # If desired, generate synthetic data instead of actual data
  #
  # NOTE: This is used to test/debug calibration process
  if (opts$data_source == "synthetic") {
    fit_data = load_synthetic(o, synthetic)
    
    # We're done, return out early
    return(fit_data)
  }
  
  # ---- Otherwise load empirical data ----
  
  message(" - Loading empirical data: ", opts$data_source)    
  
  # Pull data from ECDC
  if (opts$data_source == "ecdc") {
    
    # TODO: Load actual data
    fit_data = data.table(date          = 1 : 90, 
                          confirmed     = 35,
                          hospital_beds = 8, 
                          deaths        = 0)
  }
  
  # Pull data from xxxx
  if (opts$data_source == "xxxx") {
    
    # TODO: Load actual data
    fit_data = data.table(date      = 1 : 90, 
                          confirmed = 35)
  }
  
  # Pull data from covidbaseau
  if (opts$data_source == "covidbaseau") {
    
    # Load from url
    raw_data = read.csv(o$covidbaseau_api) 
    names(raw_data) = raw_data[1, ]
    
    # A touch of cleaning
    source_data = raw_data %>%
      select(1, 2, 4) %>%
      slice(2 : n()) %>%
      mutate(date          = format_date(Date), 
             hospital_beds = as.numeric(Hospital), 
             icu_beds      = as.numeric(ICU)) %>%
      select(date, hospital_beds, icu_beds) %>%
      mutate(confirmed = hospital_beds * 10) %>%  # Fake cases data used for testing
      as.data.table()
    
    # Apply burn in period and format for use in fitting process
    fit_data = source_data %>%
      filter(date <= opts$data_end) %>%
      slice_tail(n = opts$data_days + opts$data_burn_in) %>%
      mutate(date = 1 : n()) %>%
      pivot_longer(cols     = -date, 
                   names_to = "metric") %>%
      mutate(value = ifelse(date <= opts$data_burn_in, NA, value)) %>%
      pivot_wider(id_cols    = date, 
                  names_from = "metric") %>%
      as.data.table()
  }
  
  return(fit_data)
}

# ---------------------------------------------------------
# Load synthetic data (and generate it if it doesn't already exist)
# ---------------------------------------------------------
load_synthetic = function(o, synthetic) {
  
  # File name to save data to / load data from
  save_file = paste0(o$pth$fitting, "synthetic_data.rds")
  
  # If file exists we may want to just reload this
  if (file.exists(save_file) && !o$force_regenerate_synthetic) {
    
    message(" - Loading synthetic data")
    
    # Simply load previous saved data
    fit_data = readRDS(save_file)
    
  } else { # Otherwise run the model to generate fake data 
    
    message(" - Generating synthetic data")
    
    # Create new synthetic data by running the model
    fit_data = generate_synthetic(o, synthetic)
    
    # Save this file for potential later re-use
    saveRDS(fit_data, file = save_file)
  }
  
  return(fit_data)
}

# ---------------------------------------------------------
# Generate synthetic data by running the model
# ---------------------------------------------------------
generate_synthetic = function(o, synthetic) {
  
  # Check input is provided and is in list format
  if (!is.list(synthetic) || length(synthetic) == 0)
    stop("Synthetic parameter values must be provided as a non-trivial list")
  
  # Append flag for performing fit
  fit_list = c(synthetic, .perform_fit = TRUE)
  
  # Run model with these pre-defined parameters (see model.R)
  result = model(o, "baseline", seed = 1, fit = fit_list, verbose = "bar")
  
  # TODO: Offset this with a bit of random noise...
  
  # Extract model outcomes, this is what we'll fit to
  fit_data = result$output %>%
    select(date, metric, value) %>%
    pivot_wider(names_from = "metric", 
                values_from = value) %>%
    as.data.table()
  
  return(fit_data)
}

# ---------------------------------------------------------
# Pull OSI data from API endpoint and save in data cache
# ---------------------------------------------------------
pull_osi = function(o, y) {
  
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
  
  # Save data in cache (needed as cluster nodes don't have internet access)
  saveRDS(osi_df, file = paste0(o$pth$cache, "osi.rds"))  
}

