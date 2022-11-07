###########################################################
# LOAD DATA
#
# Load epi data that may - or may not - be required for
# calibration. Also able to generate synthetic data to test
# fitting algorithm.
#
###########################################################

# ---------------------------------------------------------
# Parent function for extracting relevant fitting data
# ---------------------------------------------------------
load_data = function(o, fit, synthetic = NULL) {
  
  # Check flag for whether we need to load data
  if (fit$.use_data == TRUE) {
    
    # Extract subset of inputs relevant for loading appropriate data
    opts = c(type = fit$input$calibration_type, fit$input$calibration_options)
    
    # Append country codes (converted from ISO2) using countrypackage package
    #
    # NOTE: countrycode is installed as a dependency of socialmixr
    opts$country_iso3 = countrycode::countrycode(opts$country, "iso2c", "iso3c")
    opts$country_name = countrycode::countrycode(opts$country, "iso2c", "country.name")
    
    # Identify data source regardless of capitalisation
    opts$data_source = lapply(opts$data_source, toupper)
    
    # Load epi data regardless of what we're fititng to
    fit = load_epi(o, opts, fit, synthetic)
    
    # Load Re estimates if necessary (may or may not use epi data)
    if (fit$.use_Re == TRUE)
      fit = load_Re(o, opts, fit)
  }
  
  # If fitting to user-defined Re, take this directly from yaml file
  if (fit$.use_data == FALSE)
    fit$target = fit$input$Re
  
  # Set trivial value if data/target irrelevant
  if (is.null(fit$target)) fit$target = NA
  if (is.null(fit$data))   fit$data   = NA
  
  # Display target Re if appropriate
  if (!is.na(fit$target))
    message("  > Re target: ", round(fit$target, digits = 2))
  
  # Quick plot of fitting data
  # g = ggplot(fit$data, aes(x = date, y = value, colour = metric)) +
  #   geom_line(size = 3, alpha = 0.5) +
  #   facet_grid(metric ~ type, scales = "free_y")
  
  return(fit)
}

# ---------------------------------------------------------
# Load emperical epidemiological data from source (or generate synthetic data)
# ---------------------------------------------------------
load_epi = function(o, opts, fit, synthetic) {
  
  # If desired, generate synthetic data instead of actual data
  #
  # NOTE: This is used to test/debug calibration process
  if (opts$data_synthetic == TRUE) {
    fit$data = load_synthetic(o, opts, synthetic)
    
    # We're done, return out early
    return(fit)
  }
  
  # Otherwise load empirical epidemiological data...
  message(" - Loading epi data: ", opts$data_source$epi)
  
  # ---- Source: ECDC ----
  
  # All dates we're interested in
  dates_df = get_data_dates(opts)
  
  # Pull data from ECDC
  if (opts$data_source$epi == "ECDC") {
    
    # Load raw case and death data from ECDC
    raw_data = read.csv(o$ecdc_api$cases, fileEncoding = "UTF-8-BOM")
    
    # Select only columns of interest and convert dates to R-interpretable
    data_cases = raw_data %>%
      select(date      = dateRep,
             confirmed = cases,
             deaths    = deaths,
             country   = geoId) %>%
      filter(country == opts$country) %>%
      select(-country) %>%
      mutate(date = format_date(date)) %>%
      # Keep only dates of interest and melt down...
      right_join(y  = dates_df,
                 by = "date") %>%
      pivot_longer(cols = c(-date, -day), 
                   names_to = "metric") %>%
      arrange(metric, date) %>%
      # Summarise time period if desired...
      change_time_period(dates_df, opts$data_period) %>%
      setDT()
    
    # Throw error if no case data for this country
    if (nrow(data_cases) == 0)
      stop("No ECDC case data found for country '", opts$country)
    
    # Load raw hospital occupancy data from ECDC
    raw_data = read.csv(o$ecdc_api$hosp, fileEncoding = "UTF-8-BOM")
    
    # We'll also scale 'per 100k' values to pop size to ensure consistency
    scale_pop = opts$country_pop / 1e5
    
    # Select only columns of interest and convert dates to R-interpretable
    data_hospital = raw_data %>%
      filter(country == opts$country_name) %>%
      mutate(date = format_date(date)) %>%
      select(date, indicator, value) %>%
      pivot_wider(id_cols    = date, 
                  names_from = indicator) %>%
      right_join(y  = dates_df,
                 by = "date") %>%
      arrange(date) %>%
      # Scale 'weekly' and 'per 100k' appropriately...
      mutate(across(contains("Weekly"),   ~ . / 7), 
             across(contains("per 100k"), ~ . * scale_pop)) %>% 
      rename(any_of(o$data_dict$ecdc)) %>%
      # Back fill weekly data to represent that value over the whole week...
      fill(any_of(names(o$data_dict$ecdc)), .direction = "up") %>% 
      pivot_longer(cols = c(-date, -day), 
                   names_to = "metric") %>%
      filter(!is.na(value)) %>%
      arrange(metric, date) %>% 
      # Summarise time period if desired...
      change_time_period(dates_df, opts$data_period) %>%
      setDT()
    
    # Throw an error if any unknown data metrics are present 
    unknown_metrics = setdiff(unique(data_hospital$metric), names(o$data_dict$ecdc))
    if (length(unknown_metrics) > 0)
      stop("Unknown data variables for country ", opts$country, 
           ": ", paste(unknown_metrics, collapse = ", "))
    
    # Throw warning if no hospital data for this country
    if (nrow(data_hospital) == 0)
      warning("No ECDC hospital data found for country ", opts$country)
    
    # Filter dates of interest for fitting and plotting
    fit$data = rbind(data_cases, data_hospital)
    
    # Check no data is negative
    if (any(fit$data$value < 0))
      stop("Negative data values identified")
  }
  
  # Apply pop scaler to every epi metric (-> per 100,000 people)
  fit$data[, value := value / scale_pop]
  
  return(fit)
}

# ---------------------------------------------------------
# Load or calculate Re estimates
# ---------------------------------------------------------
load_Re = function(o, opts, fit) {
  
  # Synthetic data has this taken care of
  if (fit$.use_synthetic == TRUE) {
    Re_values = fit$data[metric == "Re" & !is.na(value), value]
    
    # Throw an error if no Re found - likely need to regenerate
    #
    # NOTE: This happens if synthetic was previously generated for epi_data
    if (length(Re_values) == 0)
      stop("No Re identified in synthetic data - try regenerating (force_regenerate_synthetic)")
    
  } else {  # Re is to be calculated or loaded
    
    # ---- Source: EpiEstim using pre-loaded epi data ----
    
    # Sources for which we must calculate Re using EpiEstim
    calculate_Re_sources = c("ECDC")
    
    # Check whether user wishes to use such a source for Re
    if (opts$data_source$Re %in% calculate_Re_sources) {
      
      message(" - Calculating Re from epi data")
      
      # Extract incidence data we'll use for fitting
      incidence = fit$data[metric == "confirmed", value]
      
      # Ensure confirmed cases is being reported in the data
      if (length(incidence) == 0)
        stop("No incidence data identified (required when fitting 'Re_data')")
      
      # TODO: Remove hardcoding
      serial_int = list(mean = 8.5, std = 2.5)
      
      # Calculate Re, returning various details
      Re_info = calculate_Re(incidence, serial_int, fit$input$Re_window)
      
      # Extract vector of Re values
      Re_values = Re_info$Re_df %>% 
        filter(!is.na(value)) %>% 
        pull(value)
    }
    
    # ---- Source: ETH ----
    
    # Pull data from ETH
    if (opts$data_source$Re == "ETH") {
      
      message(" - Loading pre-calculated Re: ", opts$data_source$Re)
      
      # Insert the ISO3 code into url
      country_api = str_replace(o$eth_api, "<country>", opts$country_iso3)
      
      # Load raw data
      raw_data = read.csv(country_api)
      
      # Take the mean Re across all metric sources
      Re_values = raw_data %>%
        filter(estimate_type == "Cori_slidingWindow") %>%
        select(date, data_type, median_R_mean) %>%
        mutate(date = format_date(date)) %>%
        group_by(date) %>%
        summarise(value = mean(median_R_mean, na.rm = TRUE)) %>%
        ungroup() %>%
        # Take only dates of interest...
        right_join(y  = get_data_dates(opts), 
                   by = "date") %>%
        filter(value > 0) %>%
        pull(value)
    }
  }
  
  # ---- Extract Re fitting target ----
  
  # Values we wish to take the mean over (consistent with number of fit_days)
  Re_idx = length(Re_values) - fit$input$Re_fit$n_days : 0
  
  # Take the average over these days 
  fit$target = mean(Re_values[Re_idx])
  
  # Remove data we may have used to generate this Re
  fit$data = NULL
  
  return(fit)
}

# ---------------------------------------------------------
# Load synthetic data (and generate it if it doesn't already exist)
# ---------------------------------------------------------
load_synthetic = function(o, opts, synthetic) {
  
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
    fit_data = generate_synthetic(o, opts, synthetic)
    
    # Save this file for potential later re-use
    saveRDS(fit_data, file = save_file)
  }
  
  return(fit_data)
}

# ---------------------------------------------------------
# Generate synthetic data by running the model
# ---------------------------------------------------------
generate_synthetic = function(o, opts, synthetic) {
  
  # Check input is provided and is in list format
  if (!is.list(synthetic) || length(synthetic) == 0)
    stop("Synthetic parameter values must be provided as a non-trivial list")
  
  # Append flag for performing fit
  fit_list = c(synthetic, .perform_fit = TRUE)
  
  # Run model with these pre-defined parameters (see model.R)
  result = model(o, "baseline", seed = 1, fit = fit_list, verbose = "bar")
  
  # Extract model outcomes & offset with some noise - this is what we'll fit to
  fit_data = result$output %>%
    select(date, metric, value) %>%
    mutate(value = value * rnorm(n(), mean = 1, sd = opts$offset_synthetic)) %>%
    as.data.table()
  
  return(fit_data)
}

# ---------------------------------------------------------
# Filter to only dates of interest
# ---------------------------------------------------------
get_data_dates = function(opts) {
  
  # # Throw an error if we have insufficient data
  # if (max(all_data$date) < opts$data_end)
  #   stop("Data only available up to ", max(all_data$date) - 1, 
  #        " (requested ", opts$data_end, ")")
  
  # Calibrating to epi data: go back to start of burn in
  if (opts$type == "epi_data")
    days = list(fit = opts$data_days, burn = opts$data_burn_in)
  
  # Calibrating to Re: only a small window of data needed
  if (opts$type == "Re_data")
    days = list(fit = opts$data_days, burn = 0)
  
  # Dates of data we'll fit to
  date_end   = format_date(opts$data_end)
  date_start = date_end - days$fit - days$burn + 1
  
  # Sequence of dates for fitting and plotting
  all_dates  = seq(date_start,  date_end,  by = "day")
  
  # Datatable of fitting dates and day indices
  dates_df = data.table(day  = 1 : length(all_dates), 
                        date = all_dates) %>%
    filter(day > days$burn)
  
  return(dates_df)
}

# ---------------------------------------------------------
# Change time period from day to week, month, quater, or year
# ---------------------------------------------------------
change_time_period = function(data, dates_df, data_period) {
  
  # A trivial process is looking at daily data
  if (data_period != "day") {
    
    # Changing daily data to different time period 
    data %<>% 
      mutate(date = ceiling_date(date, data_period)) %>%
      group_by(metric, date) %>%
      summarise(value = sum(value)) %>% 
      ungroup() %>%
      inner_join(dates_df, by = "date") %>%
      filter(!is.na(value)) %>%
      select(date, day, metric, value) %>%
      setDT()
  }
  
  return(data)
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

