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
  
  # Initiate trivial output
  Re_df = NULL
  
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
      list[fit, Re_df] = load_Re(o, opts, fit)
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
  
  return(list(fit, Re_df))
}

# ---------------------------------------------------------
# Load emperical epidemiological data from source (or generate synthetic data)
# ---------------------------------------------------------
load_epi = function(o, opts, fit, synthetic) {
  
  # If desired, generate synthetic data instead of actual data
  #
  # NOTE: This is used to test/debug calibration process
  if (opts$data_synthetic == TRUE) {
    fit$data = load_synthetic(o, synthetic)
    
    # We're done, return out early
    return(fit)
  }
  
  # Otherwise load empirical epidemiological data...
  message(" - Loading epi data: ", opts$data_source$epi)
  
  # ---- Source: ECDC ----
  
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
      mutate(date = format_date(date)) %>%
      select(-country) %>%
      as.data.table()
    
    # Throw error if no case data for this country
    if (nrow(data_cases) == 0)
      stop("No ECDC case data found for country '", opts$country)
    
    # Load raw hospital occupancy data from ECDC
    raw_data = read.csv(o$ecdc_api$hosp, fileEncoding = "UTF-8-BOM")
    
    # Select only columns of interest and convert dates to R-interpretable
    data_hospital = raw_data %>%
      filter(country == opts$country_name) %>%
      mutate(date = format_date(date)) %>%
      select(date, indicator, value) %>%
      pivot_wider(id_cols    = date, 
                  names_from = indicator) %>%
      select(date, icu_beds = "Daily ICU occupancy", 
             hospital_beds  = "Daily hospital occupancy") %>%
      as.data.table()
    
    # Throw warning if no hospital data for this country
    if (nrow(data_hospital) == 0)
      warning("No ECDC hospital data found for country '", opts$country)
    
    # Filter dates of interest for fitting and plotting
    fit$data = rbind(filter_data_dates(opts, data_cases), 
                     filter_data_dates(opts, data_hospital)) %>%
      arrange(type, metric, date) 
    
    # Check no data is negative
    if (any(fit$data$value < 0))
      stop("Negative data values identified")
  }
  
  # Apply pop scaler to every epi metric (-> per 100,000 people)
  fit$data[, value := value / fit$pop_scaler]
  
  return(fit)
}

# ---------------------------------------------------------
# Load or calculate Re estimates
# ---------------------------------------------------------
load_Re = function(o, opts, fit) {
  
  # Sythentic data has this taken care of
  if (fit$.use_synthetic == TRUE) {
    Re_values = fit$data[!is.na(value), value]
    
    # Plotting data is irrelevant here
    fit$data = fit$data[type != "plot", ]
    
    # As is a dataframe of all Re values
    Re_df = NULL
    
  } else {  # Re is to be calculated or loaded
    
    # ---- Source: EpiEstim using pre-loaded epi data ----
    
    # Sources for which we must calculate Re using EpiEstim
    calculate_Re_sources = c("ECDC")
    
    # Check whether user wishes to use such a source for Re
    if (opts$data_source$Re %in% calculate_Re_sources) {
      
      message(" - Calculating Re from epi data")
      
      # Extract incidence data we'll use for fitting
      incidence = fit$data %>%
        filter(type   == "fit", 
               metric == "confirmed") %>%
        pull(value)
      
      # Ensure confirmed cases is being reported in the data
      if (length(incidence) == 0)
        stop("No incidence data identified (required when fitting 'Re_data')")
      
      # TODO: Remove hardcoding
      serial_int = list(mean = 8.5, std = 2.5)
      
      # Calculate Re, returning various details
      Re_info = calculate_Re(incidence, serial_int, fit$input$Re_window)
      
      # Extract vector of Re values (ignoring NAs)
      Re_df     = Re_info$Re_df %>% select(date, value)
      Re_values = as.numeric(na.omit(Re_df$value))
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
      clean_data = raw_data %>%
        filter(estimate_type == "Cori_slidingWindow") %>%
        select(date, data_type, median_R_mean) %>%
        mutate(date = format_date(date)) %>%
        group_by(date) %>%
        summarise(Re = mean(median_R_mean, na.rm = TRUE)) %>%
        mutate(Re = ifelse(is.nan(Re), NA, Re)) %>%
        as.data.table()
      
      # Filter dates of interest for fitting and plotting
      Re_data = filter_data_dates(opts, clean_data)
      
      # Extract vector of Re values
      Re_df     = Re_data[type == "plot", ] %>% select(date, value)
      Re_values = Re_data[type == "fit", value]
      
      # Only 'data' we take forward is plotting points (include Re 'data')
      fit$data = rbind(fit$data, Re_data)
    }
  }
  
  # ---- Extract Re fitting target ----
  
  # Values we wish to take the mean over (consistent with number of fit_days)
  Re_idx = length(Re_values) - fit$input$Re_fit$n_days : 0
  
  # Take the average over these days 
  fit$target = mean(Re_values[Re_idx])
  
  # Only 'data' we take forward is plotting points
  fit$data = fit$data[type == "plot", ]
  
  return(list(fit, Re_df))
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
  model_output = result$output %>%
    select(date, metric, value) %>%
    as.data.table()
  
  # Replicate data for fitting and plotting
  fit_data = rbind(mutate(model_output, type = "fit"), 
                   mutate(model_output, type = "plot"))
  
  return(fit_data)
}

# ---------------------------------------------------------
# Filter to only dates of interest
# ---------------------------------------------------------
filter_data_dates = function(opts, all_data) {
  
  # Throw an error if we have insufficient data
  if (max(all_data$date) < opts$data_end)
    stop("Data only available up to ", max(all_data$date) - 1, 
         " (requested ", opts$data_end, ")")
  
  # Calibrating to epi data: go back to start of burn in
  if (opts$type == "epi_data")
    days = list(fit = opts$data_days, burn = opts$data_burn_in)
  
  # Calibrating to Re: only a small window of data needed
  if (opts$type == "Re_data")
    days = list(fit = opts$data_days, burn = 0)
  
  # Dates of data we'll fit to
  fit_end   = format_date(opts$data_end)
  fit_start = fit_end - days$fit - days$burn + 1
  
  # Calibrating to epi data: plot this very same data
  #
  # TODO: Should we record all future data here? Then differentiate when actually plotting?
  if (opts$type == "epi_data") {
    plot_end   = fit_end # max(all_data$date)
    plot_start = fit_start
  }
  
  # Calibrating to Re: plot everything into the 'future' (after data_end)
  if (opts$type == "Re_data") {
    plot_end   = max(all_data$date)
    plot_start = fit_end + 1
  }
  
  # Sequence of dates for fitting and plotting
  fit_dates  = seq(fit_start,  fit_end,  by = "day")
  plot_dates = seq(plot_start, plot_end, by = "day")
  
  # Datatable of fitting dates and day indices
  fit_df = data.table(day  = 1 : length(fit_dates), 
                      date = fit_dates, 
                      type = "fit") %>%
    filter(day > days$burn)
  
  # Datatable of plotting dates and day indices
  plot_df = data.table(day  = 1 : length(plot_dates), 
                       date = plot_dates, 
                       type = "plot") %>%
    filter(day > days$burn)
  
  # Join the data and melt to long format
  fit_data = rbind(fit_df, plot_df) %>%
    left_join(all_data, by = "date") %>%
    pivot_longer(cols = names(all_data)[-1], 
                 names_to = "metric") %>%
    filter(!is.na(value)) %>%
    select(date = day, metric, value, type) %>%
    arrange(type, metric, date) %>%
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

