###########################################################
# LOAD DATA
#
# Load demographic, epidemiological, and response data.
#
###########################################################

# ---------------------------------------------------------
# Extract epi and response data for all cantons
# ---------------------------------------------------------
load_data = function(o, from_cache = FALSE, plot_data_sources = FALSE, quiet = FALSE) {
  
  # ---- Load from cache if desired ----
  
  # If data has been cached we can simply load the RDS file
  if (from_cache == TRUE) {
    
    if (!quiet) message(" - Loading data from cache")
    
    # Load the cached data - throw an error if file does not exist
    err_msg = "Data has not yet been cached - either load directly or ensure prior caching"
    d = try_load(o$pth$data, "data_cache", msg = err_msg, sep = TRUE) 
    
    return(d)
  } 
  
  if (!quiet) message(" - Loading data")
  
  # ---- Pull FOPH data (once per day) ----
  
  # Assume we'll need to (re)download FOPH data
  do_download = TRUE
  
  # Skip all conditions if we're forcing a new download
  if (o$forceload_foph == FALSE) {
    
    # We may not need to if it already exists
    if (file.exists(o$foph_zip)) {
      
      # Cancel the download if it's not a recent file
      if (file.info(o$foph_zip)$ctime > now() - hours(24))
        do_download = FALSE
    }
  }
  
  # Check for FOPH zip in /temp location
  if (do_download == TRUE) {
    
    if (!quiet) message("  > Pulling latest data from FOPH")
    
    # Pull the latest data and save in temp location
    quiet(download.file(o$foph_data_url, o$foph_zip))
  }
  
  # ---- Load all data for all cantons ----
  
  # Easy access temporary bool - do we want national level data
  o$national_data = "CH" %in% o$cantons
  
  # Load all data from the appropriate sources (for epi data see o$data_source)
  epi_data        = load_data_epi(o, plot_data = plot_data_sources)  # Sources: FOPH, OpenZH, and ETH github
  pop_data        = load_data_pop(o)       # Source: FOPH
  vaccine_data    = load_data_vaccine(o)   # Source: FOPH
  response_data   = load_data_response(o)  # Source: SwissTPH github
  capacity_data   = load_data_capacity(o)  # Source: icumonitoring.ch via polybox (credit Thom Van Boeckel)
  weather_daily   = load_weather_daily(o, pop_data)    # Source: Meteo Schweiz via opendata.swiss
  weather_monthly = load_weather_monthly(o, pop_data)  # Source: Meteo Schweiz
  
  # Initiate data list
  d = list()
  
  # ---- Main canton loop ----
  
  # Initiate plotting dataframes if required
  if (plot_data_sources == TRUE)
    diagnoses_plot_df = weather_plot_df = NULL
  
  # Loop through modelled cantons
  for (canton_name in o$cantons) {
    
    # ---- Demographic data ----
    
    # Demographic breakdown by age (needs to be consistent with o$age_groups)
    pop_canton = pop_data %>% 
      filter(canton == canton_name) %>% 
      select(-canton)
    
    # Convert to named vector
    #
    # NOTE: Adding 'years' to make it more obvious names are strings and must be indexed as so
    ages_canton = setNames(pop_canton$value, paste0(pop_canton$age, "years"))
    
    # Store in data list
    d[[canton_name]]$demog = ages_canton
    
    # ---- Epidemiological data ----
    
    # Epi data for this canton
    epi_canton = epi_data %>% 
      filter(canton == canton_name) %>% 
      select(-canton)
    
    # Store epi data for this canton
    d[[canton_name]]$epi = epi_canton
    
    # ---- Outbreak start date ----
    
    # Indices of dates with more than one confirmed case
    cases_pos = filter(epi_canton, metric == "confirmed", value > 1)
    cases_idx = which(o$dates_data %in% cases_pos$date)
    
    # We'll assume epidemic kicks off when we have 3 consecutive days of this
    continued_cases_idx = cases_idx[cumsum(diff(cases_idx) == 1) >= 3]
    
    # We'll consider the first time this happens as epidemic kick off
    outbreak_start = continued_cases_idx[1] - 2
    
    # Alternatively we could start more recently if manually defined
    manual_start = which(o$dates_data %in% o$manual_start_date)
    
    # We'll use the greater of these two to initiate proceedings
    d[[canton_name]]$outbreak_start = max(outbreak_start, manual_start)
    
    # ---- Diagnosis data ----
    
    # Confirmed case data for this canton
    diagnoses_canton = epi_data %>%
      filter(canton == canton_name, 
             metric == "confirmed") %>%
      select(date, n_diagnosed = value)
    
    # Amount we'll shift the smoothing to achieve centering
    n_shift = ceiling(o$diagnosis_ma_days / 2)
    
    # Fit a simple moving average model to the data
    sma_model = smooth::sma(diagnoses_canton$n_diagnosed,
                            order = o$diagnosis_ma_days,
                            h = n_shift,
                            interval = "none")
    
    # Shift the vector to recenter the moving average
    sma_centre = c(sma_model$fitted[-(1 : n_shift)], sma_model$forecast)
    
    # Store these fitted and centred moving averages
    diagnoses_canton = diagnoses_canton %>% 
      mutate(sma_fitted = round(as.vector(sma_model$fitted)), 
             sma_centre = round(sma_centre)) %>%
      right_join(data.table(date   = o$dates_data, 
                            canton = canton_name), by = "date") %>%
      arrange(date)
    
    # If we want to plot, concatenate testing outcomes
    if (plot_data_sources == TRUE)
      diagnoses_plot_df = rbind(diagnoses_plot_df, diagnoses_canton)
    
    # Simplify what we provide to the model
    d[[canton_name]]$diagnoses = diagnoses_canton %>%
      select(date, number_diagnosed = sma_centre) %>%
      replace_na(list(number_diagnosed = 0))
    
    # ---- Vaccine rollout data ----
    
    # vaccine rollout data for this canton
    vaccine_canton = vaccine_data %>% 
      filter(canton == canton_name) %>% 
      select(date, number_vaccines, number_doses)
    
    # Store vaccine data for this canton
    d[[canton_name]]$vaccine = vaccine_canton
    
    # ---- Response measure data ----
    
    # Response data for this canton
    response_canton = response_data %>% 
      filter(canton == canton_name) %>% 
      select(-canton)
    
    # Store response data for this canton
    d[[canton_name]]$response = response_canton
    
    # ---- ICU capacity data ----
    
    # Capacity data for this canton
    capacity_canton = capacity_data %>% 
      filter(canton == canton_name)
    
    # Store capacity data for this canton
    d[[canton_name]]$capacity = capacity_canton$icu_capacity
    
    # ---- Weather data ----
    
    # Temperature data for this canton
    weather_canton = weather_daily %>% 
      filter(canton == canton_name) %>% 
      select(-canton)
    
    # Amount we'll shift the smoothing to achieve centering
    n_shift = ceiling(o$weather_ma_days / 2)
    
    # Fit a simple moving average model to the data
    sma_model = smooth::sma(weather_canton$temperature, 
                            order = o$weather_ma_days, 
                            h = n_shift, 
                            interval = "none")
    
    # Shift the vector to recenter the moving average
    sma_centre = c(sma_model$fitted[-(1 : n_shift)], sma_model$forecast)
    
    # Format into a datatable of all data time points
    weather_canton = data.table(date = o$dates_data, 
                                canton = canton_name, 
                                sma_fitted = as.vector(sma_model$fitted), 
                                sma_centre = sma_centre) %>%
      left_join(weather_canton, by = "date")
    
    # If we want to plot, concatenate testing outcomes
    if (plot_data_sources == TRUE)
      weather_plot_df = rbind(weather_plot_df, weather_canton)
    
    # Simplify what we provide to the model
    d[[canton_name]]$weather = weather_canton %>%
      select(date, temperature = sma_centre)
    
    # Also append monthly averages for this canton - used for future projections
    d[[canton_name]]$weather_extrap = weather_monthly %>%
      filter(canton == canton_name) %>% 
      select(-canton)
    
    # Also store all cantonal weather if modelling at the national level
    if (o$national_data == TRUE)
      d[[canton_name]]$national_weather = weather_daily
  }
  
  # Perform some basic sanity checks on the data we've just compiled
  data_checks(o, d)
  
  # Save the data to file - this can then be used by cluster jobs
  saveRDS(d, file = file.path(o$pth$data, "data_cache.rds"))
  
  # Finally, plot testing and temperature data if required (see plotting.R)
  #
  # NOTE: We do this here rather than when we load the raw data as we want to
  #       examine the smoothed, interpolated, and extrapolated equivalents
  if (plot_data_sources == TRUE) {
    plot_smoothed_data(o, weather_plot_df,   "temperature", "Temperature")
    plot_smoothed_data(o, diagnoses_plot_df, "n_diagnosed", "Number of diagnosis")
  }
  
  return(d)
}

# ---------------------------------------------------------
# Load epi data - anything that can aligned to when calibrating
# ---------------------------------------------------------
load_data_epi = function(o, plot_data = FALSE, by_age = FALSE) {
  
  # Initiate data template of all date-canton combinations
  data_df = base::merge(data.frame(date   = o$dates_data), 
                        data.frame(canton = o$all_switzerland)) # all_cantons
  
  if (by_age == TRUE)
    data_df = base::merge(data_df, data.frame(age = o$age_groups))
  
  # Use this to initiate full datatables for the different sources
  foph_data = openzh_data = as.data.table(data_df)
  
  # ---- FOPH data ----
  
  # Loop through metrics in FOPH data dictionary
  for (load_metric in names(o$foph_dict)) {
    
    if (by_age == TRUE) {
      
      # File name of interest
      foph_file = paste0("data/COVID19", o$foph_dict[[load_metric]], "_geoRegion_AKL10_w.csv")
      
      # Load the data and remove what we don't need
      this_data = data.table(read.csv(unz(o$foph_zip, foph_file))) %>%
        select(datum, geoRegion, altersklasse_covid19, entries) %>%
        setNames(c("date", "canton", "age", load_metric)) %>%
        mutate(date = as.character(date),
               year = format_date(paste0(substr(date, 1, 4), "-01-01")), 
               week = as.numeric(substr(date, 5, 6)), 
               week = as.character(lubridate::weeks(week)), 
               date = year + period(week, units = "week"),
               date = round_date(date, "week")) %>%
        select(-year, -week) %>%
        mutate(date   = format_date(date), 
               canton = as.character(canton), 
               age    = str_remove_all(age, " "))
      
      # Join with full dates, cantons, and age groups
      foph_data = foph_data %>%
        left_join(this_data, by = c("date", "canton", "age"))
      
    } else {
      
      # File name of interest
      foph_file = paste0("data/COVID19", o$foph_dict[[load_metric]], "_geoRegion.csv")
      
      # Load the data and remove what we don't need
      this_data = data.table(read.csv(unz(o$foph_zip, foph_file))) %>%
        select(datum, geoRegion, entries) %>%
        setNames(c("date", "canton", load_metric)) %>%
        mutate(date   = format_date(date), 
               canton = as.character(canton))
      
      # Join with full dates and cantons
      foph_data = foph_data %>%
        left_join(this_data, by = c("date", "canton"))
    }
    
    # Close connection to the file
    close(unz(o$foph_zip, foph_file))
  }
  
  if (by_age == TRUE)
    return(foph_data)
  
  # ---- OpenZH data ----
  
  # Loop through metrics in OpenZH data dictionary
  for (load_metric in names(o$openzh_dict)) {
    
    # File name of interest
    openzh_file = paste0("covid19_", o$openzh_dict[[load_metric]], "_switzerland_openzh.csv")
    
    # Load data for this metric and select cantons
    #
    # NOTE: Apply date filter here as occasionally get negative (ie non-monotonic) values for todays data
    this_data = read.csv(paste0(o$openzh_path, openzh_file)) %>%
      mutate(date = format_date(Date)) %>%
      filter(date <= max(o$dates_data)) %>%
      select(date, one_of(o$all_switzerland))
    
    # If cumulative, convert to temporal
    if (grepl("^cumulative", load_metric))
      this_data = mutate_at(.tbl = this_data, 
                            .vars = vars(o$all_switzerland), 
                            .funs = cum_to_temporal, o)
    
    # Melt to long format
    this_data = this_data %>%
      pivot_longer(cols      = -date, 
                   names_to  = "canton", 
                   values_to = str_remove(load_metric, "cumulative_"))
    
    # Join with full dates and cantons
    openzh_data = openzh_data %>%
      left_join(this_data, by = c("date", "canton"))
  }
  
  # ---- Combine sources ----
  
  # Data variables that are common across data sets (ie data for which we choose the source)
  common_vars = intersect(names(o$data_source), intersect(names(foph_data), names(openzh_data)))
  unique_vars = setdiff(names(o$data_source), common_vars)
  
  # Sources chosen for these data variables
  chosen_vars = paste0(common_vars, ".", o$data_source[common_vars])
  
  # Select variables from sources as requested and melt to long format
  epi_data = full_join(foph_data, openzh_data,
                       by = c("date", "canton"), 
                       suffix = c(".foph", ".openzh")) %>%
    select(date, canton, one_of(c(unique_vars, chosen_vars))) %>%
    setNames(c("date", "canton", unique_vars, common_vars)) %>%
    select(one_of(c("date", "canton", names(o$data_source)))) %>%
    pivot_longer(cols = -c("date", "canton"), 
                 names_to = "metric") %>%
    mutate(grouping = "none", group = NA) %>%
    select(date, canton, metric, grouping, group, value) %>%
    as.data.table()
  
  # ---- Additional data / estimates ----
  
  # Load variant data - assume representative for all cantons
  variant_data = fread(o$variant_url, showProgress = FALSE) %>% 
    rename(B117 = b117, S501Y_V2 = s501yv2) %>% 
    replace_na(list(B117 = 0, S501Y_V2 = 0)) %>% 
    mutate(D614G = pmax(n - B117 - S501Y_V2, 0), 
           year = format_date(paste0(year, "-01-01")), 
           week = lubridate::weeks(week - 1), 
           date = round_date(year + week, "week")) %>%
    select(-n, -year, -week) %>%
    pivot_longer(cols = o$variant_names, 
                 names_to = "group") %>%
    group_by(date, group) %>%
    summarise(cases = sum(value)) %>%
    group_by(date) %>%
    mutate(value = 100 * cases / sum(cases), 
           metric = "variant_prevalence", 
           grouping = "variant") %>%
    select(-cases) %>%
    expand_grid(canton = o$cantons) %>%
    select(names(epi_data))
  
  # Append to epi datatable
  epi_data = rbind(epi_data, variant_data)
  
  # Load seroprevalence estimates from manually downloaded csv file
  seroprev_data = read.csv(o$pth$seroprev) %>% 
    filter(canton %in% o$cantons) %>% 
    mutate(date = format_date(date)) %>% 
    pivot_longer(cols = seroprevalence, 
                 names_to = "metric") %>% 
    mutate(grouping = "na", group = NA)
  
  # Append to epi datatable
  epi_data = rbind(epi_data, seroprev_data)
  
  # ---- Final formatting touches ----
  
  # Remove trivial entries and sort
  epi_data = epi_data %>%
    filter(!is.na(value)) %>%
    arrange(canton, metric, grouping, group, date)
  
  # Examine the two main data sources - won't want to do this often
  if (plot_data == TRUE) {
    plot_epi_data(o, foph_data, openzh_data)  # See plotting.R
    
    # Also plot testing patterns over time
    plot_diagnosis_ratio(o, epi_data)
  }
  
  return(epi_data)
}

# ---------------------------------------------------------
# Load population data
# ---------------------------------------------------------
load_data_pop = function(o, force_national = FALSE) {
  
  # TODO: Plot the demography pyramid when plot_data_sources is on
  
  # Load data, limit max age, and summarise over gender
  this_data = rio::import(o$demog_url, which = "Population nach Alter, sex, KTN") %>%
    setNames(qc(canton, gender, age_num, value)) %>%
    filter(canton %in% o$all_cantons) %>%
    mutate(age = pmin(age_num, max(o$all_ages))) %>%
    group_by(canton, age) %>%
    summarise(value = sum(value))
  
  # A bit more to do for national level
  if (isTRUE(o$national_data) || isTRUE(force_national)) {
    
    # Summaries over all Swiss cantons
    swiss_data = this_data %>% 
      group_by(age) %>%
      summarise(value = sum(value)) %>%
      mutate(canton = "CH")
    
    # Also report back pop data for each canton as we need to pop weight other data types
    if (isTRUE(o$national_data))
      this_data = bind_rows(swiss_data, this_data)
    
    # TODO: This is a temporary fix, can be removed when hardcoding of risk groups in parameters.R has been removed
    
    # TEMP: If forcing, take only national level data
    if (isTRUE(force_national))
      this_data = swiss_data
  }
  
  # Final formatting touches
  pop_data = this_data %>%
    select(canton, age, value) %>%
    as.data.table()
  
  return(pop_data)
}

# ---------------------------------------------------------
# Load vaccination rollout data
# ---------------------------------------------------------
load_data_vaccine = function(o, by_age = FALSE) {
  
  # Default behaviour: not by age
  if (by_age == FALSE) {
    
    # File names of interest
    vaccine_file = "data/COVID19VaccPersons_v2.csv"
    n_doses_file = "data/COVID19VaccDosesAdministered.csv"
    
    # Load the data and remove what we don't need
    all_vaccinated = data.table(read.csv(unz(o$foph_zip, vaccine_file))) %>%
      filter(age_group == "total_population", 
             type == "COVID19AtLeastOneDosePersons") %>%
      select(date, canton = geoRegion, 
             number_vaccines = entries) %>%
      filter(canton %in% o$cantons) %>%
      mutate(date = format_date(date))
    
    # Load the data and remove what we don't need
    all_doses = data.table(read.csv(unz(o$foph_zip, n_doses_file))) %>%
      select(date, canton = geoRegion, 
             number_doses = entries) %>%
      filter(canton %in% o$cantons) %>%
      mutate(date = format_date(date))
    
    # Close connection to the files
    close(unz(o$foph_zip, vaccine_file))
    close(unz(o$foph_zip, n_doses_file))
    
    # Join with full combination of dates and cantons
    vaccine_data = data.frame(date = o$dates_data) %>%
      merge(data.frame(canton = o$cantons)) %>%
      left_join(all_vaccinated, by = c("date", "canton")) %>%
      left_join(all_doses, by = c("date", "canton")) %>%
      replace_na(list(number_vaccines = 0, 
                      number_doses = 0)) %>%
      as.data.table() 
  }
  
  # Alternative behaviour: by age
  if (by_age == TRUE) {
    
    # File name of interest
    age_file = "data/COVID19VaccPersons_AKL10_w_v2.csv"
    
    # Load the data and remove what we don't need
    vaccine_age = data.table(read.csv(unz(o$foph_zip, age_file))) %>%
      filter(age_group_type == "age_group_AKL10", 
             type == "COVID19AtLeastOneDosePersons") %>%
      select(date, canton = geoRegion, 
             age = altersklasse_covid19, 
             number_vaccines = entries) %>%
      filter(canton %in% o$cantons) %>%
      mutate(date = as.character(date),
             year = format_date(paste0(substr(date, 1, 4), "-01-01")),
             week = as.numeric(substr(date, 5, 6)),
             week = as.character(lubridate::weeks(week)),
             date = year + period(week, units = "week"),
             date = round_date(date, "week")) %>%
      select(-year, -week) %>%
      mutate(date   = format_date(date),
             canton = as.character(canton),
             age    = str_remove_all(age, " "))
    
    # Join with full combination of dates and cantons
    vaccine_data = expand_grid(date = o$dates_data, 
                               canton = o$cantons, 
                               age = unique(vaccine_age$age)) %>%
      left_join(vaccine_age, by = c("date", "canton", "age")) %>%
      as.data.table() 
    
    # Close connection to the file
    close(unz(o$foph_zip, age_file))
  }
  
  return(vaccine_data)
}

# ---------------------------------------------------------
# Load response data
# ---------------------------------------------------------
load_data_response = function(o) {
  
  # URL to pull from depends on whether we want national or cantonal data
  if (o$national_data) url = o$github_national_url else url = o$github_data_url
  
  # Simply load the csv from GitHub using datatables fread function
  this_data = fread(url, showProgress = FALSE) %>% 
    mutate(date = format_date(date), 
           osi_level = oxford_containment_health / 100) %>%
    select(date, canton, osi_level) %>%
    arrange(canton, date)
  
  # Join with full combination of dates and cantons
  response_data = data.frame(date = o$dates_data) %>%
    merge(data.frame(canton = o$cantons)) %>%
    left_join(this_data, by = c("date", "canton")) %>%
    as.data.table() 
  
  return(response_data)
}

# ---------------------------------------------------------
# Load capacity data
# ---------------------------------------------------------
load_data_capacity = function(o) {
  
  # Load capacity from manually downloaded csv file
  capacity_data = read.csv(o$pth$capacity) %>% 
    filter(canton %in% o$cantons) %>% 
    as.data.table()
  
  return(capacity_data)
}

# ---------------------------------------------------------
# Load daily weather data - used for dates_data
# ---------------------------------------------------------
load_weather_daily = function(o, pop_data) {
  
  # Some extra work needed for national case - use all data available
  if (o$national_data) cantons = o$all_cantons else cantons = o$cantons
  
  # Load station details (includes URLs to data)
  #
  # NOTE: Canton names already linked here so no need for o$pth$stations
  stations = read.delim(o$stations_url, sep = ";", encoding = "latin1") %>%
    select(station = Station, 
           canton  = Canton, 
           url_current = URL.Current.year,
           url_past    = URL.Previous.years..verified.data.) %>%
    filter(canton %in% cantons)
  
  # Throw an error if no data at all
  if (nrow(stations) == 0)
    stop("No weather data identified for cantons: ", paste0(cantons, collapse = ", "))
  
  # Loop over stations to extract data for
  station_data = list()
  for (i in 1 : nrow(stations)) {
    this_station = stations[i, ]
    
    # Load current and past data (need past data for 2020 weather)
    current_data = fread(this_station$url_current, showProgress = FALSE)
    past_data    = fread(this_station$url_past,    showProgress = FALSE)
    
    # Load and format data from this station
    this_data = rbindlist(list(past_data, current_data)) %>% 
      mutate(date    = format_date(date), 
             station = this_station$station, 
             canton  = this_station$canton) %>%
      select(date, station, canton, temperature = tre200dx) %>%
      filter(date >= min(o$dates_data), 
             date <= max(o$dates_data)) %>%
      mutate(temperature = as.numeric(temperature))
    
    # Row bind as we iterate
    station_data[[i]] = this_data
  }
  
  # Some cantons have multiple stations - deal with this by taking the mean across dates
  canton_data = rbindlist(station_data) %>% 
    group_by(date, canton) %>% 
    summarise(temperature = mean(temperature))
  
  # Join with full combination of dates and cantons
  daily_data = data.frame(date = o$dates_data) %>%
    merge(data.frame(canton = cantons)) %>%
    left_join(canton_data, by = c("date", "canton")) %>%
    as.data.table()
  
  # Some extra work needed for national case
  if (o$national_data) {
    
    # Population size of all cantons with weather data
    canton_pop = pop_data %>% 
      filter(canton %in% unique(canton_data$canton)) %>%
      group_by(canton) %>%
      summarize(pop_size = sum(value)) %>%
      as.data.table()
    
    # Total population of all cantons with weather data
    total_pop = sum(canton_pop$pop_size)
    
    # Take the population weighted temperature across all cantons with data
    swiss_data = daily_data %>%
      left_join(canton_pop, by = "canton") %>%
      mutate(weighted = temperature * pop_size / total_pop) %>%
      group_by(date) %>% 
      summarise(swiss_temperature = sum(weighted, na.rm = TRUE)) %>%
      mutate(canton = "CH") %>%
      select(date, canton, temperature = swiss_temperature) %>%
      as.data.table() 
    
    # Append national level weather data
    daily_data = rbind(daily_data, swiss_data)
  }
  
  # Throw an error if we are modelling cantons with no weather data
  missing_cantons = setdiff(o$cantons, unique(daily_data$canton))
  if (length(missing_cantons) > 0)
    stop("No weather data identified for cantons: ", paste0(missing_cantons, collapse = ", "))
  
  return(daily_data)
}

# ---------------------------------------------------------
# Load monthly weather data - used for future projections
# ---------------------------------------------------------
load_weather_monthly = function(o, pop_data) {
  
  # Some extra work needed for national case - use all data available
  if (o$national_data) cantons = o$all_cantons else cantons = o$cantons
  
  # Download zip container if required
  if (!file.exists(o$climate_zip))
    quiet(download.file(o$climate_url, o$climate_zip))
  
  # Link each weather station to a canton
  stations_dict = read.csv(o$pth$stations, sep = ";")
  
  # Abbreviated month names in German - tried to find a nicer solution, gave up
  months_de = setNames(1 : 12, qc(Jan, Feb, Mar, Apr, Mai, Jun, Jul, Aug, Sep, Okt, Nov, Dez))
  
  # Load station data, link to cantons, and take the mean for each canton
  monthly_data = read.delim(unz(o$climate_zip, o$climate_file), 
                            encoding = "latin1", skip = 8) %>%
    mutate(station = str_replace_all(Station, " ", "_")) %>%
    left_join(stations_dict, by = "station") %>%
    filter(canton %in% cantons) %>%
    pivot_longer(cols = names(months_de), 
                 names_to = "month") %>%
    group_by(canton, month) %>%
    summarise(temperature = mean(value)) %>%
    mutate(month = months_de[month]) %>%
    arrange(canton, month) %>%
    as.data.table()
  
  # Close the file connection
  close(unz(o$climate_zip, o$climate_file))
  
  # Some extra work needed for national case
  if (o$national_data) {
    
    # Population size of all cantons with weather data
    canton_pop = pop_data %>%
      filter(canton %in% unique(monthly_data$canton)) %>%
      group_by(canton) %>%
      summarize(pop_size = sum(value)) %>%
      as.data.table()
    
    # Total population of all cantons with weather data
    total_pop = sum(canton_pop$pop_size)
    
    # Take the population weighted temperature across all cantons with data
    swiss_data = monthly_data %>%
      left_join(canton_pop, by = "canton") %>%
      mutate(weighted = temperature * pop_size / total_pop) %>%
      group_by(month) %>% 
      summarise(swiss_temperature = sum(weighted, na.rm = TRUE)) %>%
      mutate(canton = "CH") %>%
      select(canton, month, temperature = swiss_temperature) %>%
      as.data.table() 
    
    # plot_df = rbind(mutate(monthly_data, national = FALSE),
    #                 mutate(swiss_data,   national = TRUE)) %>%
    #   mutate(canton = factor(canton, levels = c(o$all_cantons, "CH")))
    # 
    # g = ggplot(plot_df, aes(x = month, y = temperature, group = canton)) +
    #   geom_line(aes(colour = national), size = 1.5) +
    #   scale_colour_manual(values = c("grey", "blue"))
    
    # Rename the datatable and we're done
    monthly_data = swiss_data
  }
  
  # Throw an error if we are modelling cantons with no monthly data
  #
  # NOTE: We have more stations here, so likely a subset of the missing cantons previously identified
  missing_cantons = setdiff(o$cantons, unique(monthly_data$canton))
  if (length(missing_cantons) > 0)
    stop("No monthly climate data identified for cantons: ", paste0(missing_cantons, collapse = ", "))
  
  return(monthly_data)
}

# ---------------------------------------------------------
# Convert data from cumulative to temporal
# ---------------------------------------------------------
cum_to_temporal = function(x, o) {
  
  # Throw an error if input is not a vector
  if (!is.vector(x))
    stop("Input must be a vector")
  
  # Shift all points to the right and front pad with zero
  shifted = c(0, x[1 : (length(x) - 1)])
  
  # We then want the difference between these vectors
  x_temporal = x - shifted
  
  # Check for any negative values here
  negative_idx = x_temporal < 0
  
  # Check flag for removing negative data points 
  if (any(negative_idx, na.rm = TRUE) & !is.null(o$set_negative_data)) {
    warning("Removing negative data point(s)")
    
    # Set this point to some default (NA or zero)
    x_temporal[x_temporal < 0] = NA # o$set_negative_data
  }
  
  return(x_temporal)
}

# ---------------------------------------------------------
# Perform some basic sanity checks on the compiled data
# ---------------------------------------------------------
data_checks = function(o, d) {
  
  # TODO: Ideally add some checks on size and value of 'epi' data here
  
  # Loop through the cantons
  for (canton in o$cantons) {
    data = d[[canton]]
    
    # Easy reference canton name if we need to through an error
    canton_ref = paste0(" (canton '", canton, "')")
    
    # Check we have non-trivial demography data
    if (length(data$demog) != length(o$all_ages))
      stop("Inconsistent age pyramid identified in demography data", canton_ref)
    
    # Check we have a non-trivial epidemic start date
    if (is.na(data$outbreak_start))
      stop("No outbreak start date identified", canton_ref)
    
    # Check we have a non-trivial epidemic start date
    if (length(data$capacity) == 0)
      stop("No ICU capacity data identified", canton_ref)
    
    # Check monthly weather data for one row a month and non-trivial values
    if (nrow(data$weather_extrap) != 12 || any(is.na(data$weather_extrap)))
      stop("Monthly weather data not well defined", canton_ref)
    
    # Loop through data sets compiled as dataframes
    for (data_type in qc(diagnoses, response, weather)) { # capacity
      data_df = data[[data_type]]
      
      # Dataframe should have exactly the same number of rows as data days (o$dates_data)
      if (nrow(data_df) != length(o$dates_data))
        stop("Data type '", data_type, "' has inconsistent size", canton_ref)
      
      # If all entries (excluding dates) are NA then we probably have a problem
      if (all(is.na(select(data_df, -date))))
        stop("Data type '", data_type, "' has only trivial entires", canton_ref)
    }
  }
  
  # All checks passed, we're good to continue
}

# ---------------------------------------------------------
# Run the model with default parameters and fake the data accordingly
# ---------------------------------------------------------
create_synthetic_data = function(o, d) {
  
  message(" - Generating synthetic data")
  
  stop("This function needs to be updated before it is used again")
  
  # Set the random number seed generator if desired
  if (o$emulator_reproducible == TRUE)
    set.seed(1)
  
  # Each calibration parameter needs to be defined in known_params
  if (!identical(o$calibration_df$names, names(o$known_params)))
    stop("Input 'known_params' (a named vector) is not well defined")
  
  # Loop through cantons to overwrite data
  for (canton in o$cantons) {
    
    # Initiate a dataframe for the fake data we'll create - this will overwrite d$epi
    synthetic_data = data.frame(date = o$dates_all, canton = canton)
    
    # Assign these model parameters and run model
    p = get_parameters(o, d[[canton]], o$known_params)  # See parameters.R
    m = model(o, p)  # See model.R
    
    browser() # TODO: Model outcomes now long in format...
    
    # First draw *actual* count then sum those over time - think the below is correct
    
    # Create fake data with Poisson draws from simulation output
    for (metric in d$metrics)
      synthetic_data[[metric]] = rpois(nrow(m), m[[metric]])
    
    browser() # TODO: Data now long in format...
    
    # Overwrite epi data with this synthetic data
    d[[canton]]$epi = synthetic_data
  }
  
  # Store the known parameters in d
  d$known_params = o$known_params
  
  # Overwrite the cached file for use on the cluster
  saveRDS(d, file = file.path(o$pth$data, "data_cache.rds"))
  
  return(d)
}

# ---------------------------------------------------------
# Shortcut to visualise cantonal data
# ---------------------------------------------------------
run_data_exploration = function(o) {
  
  # Only continue if specified by do_step
  if (!is.element(0, o$do_step)) return()
  
  message("* Performing inital data exploration")
  
  # Call main function to produce plot
  load_data(o, plot_data_sources = TRUE)
}

