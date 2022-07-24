###########################################################
# OPTIONS
#
# Set key options for all things model related.
#
# Any numerical, logical, or character value or vector defined
# in the main set_options() function can be overridden through
# the non-mandatory \config\my_options.csv file. Note that
# this my_options file is git ignored (indeed, that is the very
# point of such a file), so each user will need to create one.
#
###########################################################

# ---------------------------------------------------------
# Set model options and assumptions
# ---------------------------------------------------------
set_options = function(do_step = NA, test_param = NULL, quiet = FALSE) {
  
  if (!quiet) message("* Setting options")
  
  # Reset R's most annoying default options
  options(stringsAsFactors = FALSE, scipen = 999, dplyr.summarise.inform = FALSE)
  
  # Initiate options list
  o = list(do_step = do_step)
  
  # Detect user - required for checking cluster jobs
  o$user = Sys.info()[["user"]]
  
  # ---- Analysis settings ---- 
  
  # Define a analysis name to differentiate previous analyses / results
  o$analysis_name = "paper_march2021"
  
  # Define a calibration file to be used
  o$calibration_name = "fit_20210308"
  
  # Swiss cantons to simulate ('all' for all cantons together, 'CH' for national level)
  o$cantons = "CH"
  
  # NOTE: Avoid using period symbols (.) in scenario and strategy names, and ensure all names are unique...
  
  # Define alternative scenarios (use "none" for baseline only, see scenarios.R)
  o$scenarios_past   = "none"  # eg: "past_evaluation"
  o$scenarios_future = "none"  # eg: c("relax_npi", "remove_npi")
  
  # Define intervention strategies (use "none" to turn off, see strategies.R)
  o$strategies = "none"  # eg: "tf8_astrazeneca"
  
  # Set analysis name and create output directory system
  o = set_dirs(o)  # See directories.R
  
  # ---- Time parameters ----
  
  # Manual start date (overwritten if before epidemic outbreak)
  o$manual_start_date = format_date("2020-01-01")
  
  # Use data up to and including this date
  o$data_end_date = format_date("2021-03-05")
  
  # Number of days to project into the future
  o$n_future_days = 180
  
  # ---- Data sources ----
  
  # Define data source for available metrics
  o$data_source = list(confirmed           = "openzh",  # "foph" or "openzh" (data looks similar)
                       deaths              = "openzh",  # "foph" or "openzh" (data looks similar)
                       hospital_admissions = "foph",    # "foph" or "openzh" (OpenZH has only a few cantons)
                       hospital_beds       = "openzh",  # "openzh" only
                       icu_beds            = "openzh")  # "openzh" only
  
  # TODO: Could also use cumulative_hospital_released from OpenZH - would need to output this from the model
  
  # URL for all epi, response, and some testing data - put together by RD, available from github
  o$github_data_url = "https://raw.githubusercontent.com/SwissTPH/COVID_measures_by_canton/master/COVID_measures_KT_epi_openZH.csv"
  
  # We can also use the Swiss TPH github for national level response measures
  o$github_national_url = "https://raw.githubusercontent.com/SwissTPH/COVID_measures_by_canton/master/COVID_measures_CH_only.csv"
  
  # URL for more detailed (and complete) testing data from FOPH
  #
  # TODO: Need a better solution for this - the link changes each day
  o$foph_data_url = "https://www.covid19.admin.ch/api/data/20211027-weppzp5n/downloads/sources-csv.zip"
  
  # Dictonary for metrics and their respective csv files
  o$foph_dict = qc(confirmed = "Cases", 	
                   deaths    = "Death", 
                   hospital_admissions = "Hosp")
  
  # Force download new FOPH data, otherwise this is done every 24 hours
  o$forceload_foph = FALSE  # Useful if just updated foph_data_url
  
  # Path stem to pull OpenZH data directly
  o$openzh_path = "https://raw.githubusercontent.com/daenuprobst/covid19-cases-switzerland/master/"
  
  # URLs for OpenZH ICU and hospitalisation data
  # 
  # TODO: Variable 'released' could also be used for calibrating to
  o$openzh_dict = qc(cumulative_confirmed = "cases", 
                     cumulative_deaths    = "fatalities", 
                     hospital_beds = "hospitalized", # TODO: Investigate difference with 'hospitalized_total'
                     icu_beds      = "icu")
  
  # URL to variant prevalence data
  o$variant_url = "https://raw.githubusercontent.com/covid-19-Re/variantPlot/master/data/data.csv"
  
  # URL to canton demographic data
  o$demog_url = paste0("https://www.bag.admin.ch/dam/bag/en/dokumente/mt/k-und-i/aktuelle-ausbrueche-pandemien/", 
                       "2019-nCoV/covid-19-basisdaten-bevoelkerungszahlen.xlsx.download.xlsx/Population_Size_BFS.xlsx")
  
  # URL for weather station information
  o$stations_url = "https://data.geo.admin.ch/ch.meteoschweiz.klima/nbcn-tageswerte/liste-download-nbcn-d.csv"
  
  # URL for monthly climate data
  o$climate_url  = "https://data.geo.admin.ch/ch.meteoschweiz.klima/normwerte/normwerte.zip"
  o$climate_file = "np8110/nvrep_np8110_tre2dymx_d.txt"
  
  # Teporary location to store zipped folders
  o$foph_zip    = file.path("", "tmp", "foph.zip")
  o$climate_zip = file.path("", "tmp", "meteo.zip")
  
  # Deal with negative data by overwriting with some value (NA or 0 would be best)
  #
  # NOTE: This may happen as we need to convert from cumulative, so may find some issues
  o$set_negative_data = NA  # Set to NULL to turn off (not recommended)
  
  # Number of days for moving average smoothing
  o$diagnosis_ma_days = 7   # Number of confirmed cases
  o$weather_ma_days   = 14  # Temperature data
  
  # ---- Model structure ----
  
  # Load disease, care, and prognosis states (see load_input.R)
  o = format_model_states(o)
  
  # TODO: As we're now improving age-related processes, o$age_groups can likely go, but o$all_ages will need to remain
  
  # Age structure (in years)
  o$age_groups = seq(10, 90, by = 10)  # Defines upper bound of bin
  o$all_ages   = 0 : (max(o$age_groups) - 1)
  
  # Define type of contact network to create
  #
  # OPTIONS:
  #  "age"    := Age structured network informed by socialmixr contact matrices
  #  "random" := Simple random network with no particular structure
  o$network_structure = "age"  # This is case sensitive
  
  # Which countries from the POLYMOD study should be used to create the contact matrix, as a character vector
  o$contact_matrix_countries = c("Germany", "Italy", "France")
  
  # Maximum number of individuals to model - if required a scaler is applied to satisfy this constraint
  o$max_population_size = 9e5  # 300k-400k is good for cantons, ~1m for national level
  
  # Replace any people that have died with new susceptibles
  o$replace_dead = FALSE  # TODO: Implement this in model.R
  
  # Use viral load to quantify infectiousness (rather than a constant)
  o$viral_load_infectivity = TRUE
  
  # Load and format SARS-COV-2 variant properties
  o = format_variant_properties(o)
  
  # ---- Vaccine settings ----
  
  # Assumed effect of vaccination
  o$vaccine_effect = "transmission"  # OPTIONS: "transmission", "severity", or "mixed"
  
  # If 'mixed', define scalers for effect on transmission blocking and disease prevention
  #
  # NOTE: If vaccine_effect != "mixed", these factors are trivial (see parameters.R)
  o$transmission_blocking_scaler = 0.8
  o$disease_prevention_scaler    = 0.75
  
  # Define vaccine priority groups - not sure this really belongs here
  #
  # NOTES: 
  #  1) Must be valid syntax as this is directly parsed and evaluated in model.R
  #  2) Groups should be mutually exclusive to avoid overwriting
  #  3) There is no need to include everyone in these groups (eg <18 may not be vaccinated)
  o$vaccine_priority = c(group1 = 'age >= 75',  # 756,400 people
                         group2 = '(age >= 65 & age < 75) | healthcare_worker == TRUE | comorbidities == TRUE',  # 850,000 + 621,600 + 560,000 people
                         group3 = 'age >= 51 & age < 64',  # 1,243,000 people
                         group4 = 'age >= 50 & age < 51',  # Those in 'communal facilities' (100,000 people)
                         group5 = 'age >= 18 & age < 50')  # All others (4,414,600 people)
  
  # Vaccine acceptance rates within each risk group (complement will not recieve a vaccine)
  o$vaccine_acceptance = c(group1 = 0.75, 
                           group2 = 0.75, 
                           group3 = 0.75,
                           group4 = 0.75,
                           group5 = 0.75)
  
  # Specify conditions for which people CANNOT receive a vaccine
  o$vaccine_criteria = "vaccine_priority != 'none' & is.na(days_vaccinated)"
  
  # Start vaccination campaign (needed as data is NOT YET backdated)
  o$vaccine_start = format_date("2021-01-01")
  
  # Load and format properties of the various vaccines available
  o = format_vaccine_properties(o)
  
  # ---- Other intervention settings ----
  
  # Future assumption of proportion of infected cases getting diagnosed
  #
  # OPTIONS:
  #  "linear_model" := Linear projection over time
  #  "constant"     := Project a constant which is the mean of the last 14 days
  o$future_diagnoses = "constant"
  
  # Date at which improved care was introduced in Switzerland
  #
  # TODO: This should be uploaded as 'data' and not as an 'option'
  o$improved_care_date = format_date("2020-06-01")
  
  # ---- Emulator settings ----
  
  # Set random seed generator prior to training emulator
  o$emulator_reproducible = TRUE
  
  # Number of parameter sets to initially sample
  o$n_init_samples = 4000
  
  # GPs are not scalable with large parameter sets - limit what we train with
  #
  # NOTE: This only comes into play if n_init_samples > max_train_samples
  o$max_train_samples = 1000
  
  # Number of seeds to run for each parameter set 
  #
  # NOTE: We'll use the same value when performing adaptive sampling
  o$n_seeds_emulator = 4
  
  # Proportion of points to hold back for training
  o$train_test_split = 0.15
  
  # Select GP kernel function
  o$gp_kernel = "Matern5_2" # "Gaussian", "Matern5_2", or "Matern3_2"
  
  # Maximum number of iterations of GP algorithm
  o$gp_max_iter = 10000
  
  # ---- Adaptive sampling settings ----
  
  # Maximum number of adaptive sampling iterations to perform
  o$resample_iters = 2  # Set to zero to turn off
  
  # Acquisition function used to evaluate points for resampling
  #
  # OPTIONS: "expected_improvement" or xxx
  o$acquisition_funcion = "expected_improvement" # TODO: Just using EI for now
  
  # Number of randomly sampled points to evaluate for resampling
  o$acquisition_points = 1e5
  
  # Limit the number of points resampled
  o$min_value_accepted   = 0 # 1e3  # Reject anything below this acquisition function value
  o$max_samples_accepted = 100  # Upper limit for points to accept on each adaptive sampling iteration
  
  # Limit how close these points can be
  o$distance_method = "manhattan"  # Manhattan better than Euclidean for 'high' dimensions
  o$distance_tolerance = 1e-3  # Function of number of calibrated parameters (x*10^sqrt(n))
  
  # ---- Calibration settings ----
  
  # Load parameter file and extend by canton (see load_input.R)
  o = format_calibration_parameters(o, test_param)
  
  # Set random seed generator prior to calibration
  o$calibration_reproducible = TRUE
  
  # Number of MCMC chains to generate
  o$mcmc_n_chains = 10
  
  # Flag to use the cluster when mcmc_n_chains > 1 (ignored if mcmc_n_chains = 1)
  o$mcmc_cluster = FALSE
  
  # Number of MCMC iterations (adaptive period = mcmc_opt_stages x mcmc_opt_freq)
  o$mcmc_iters      = 10000
  o$mcmc_opt_freq   = 500
  o$mcmc_opt_stages = 10
  
  # Initial step size
  o$mcmc_init_steps = 1
  
  # Define values of parameters we want to recover when fitting to synthetic data
  #
  # NOTES:
  #  1) Use for calibration testing purposes - set to NULL to turn off
  #  2) Each "unfixed" parameter in o$calibration_df should have a value defined here
  #  3) Non-global parameters must be defined as canton-specific (eg xxx_ZH)
  #  4) Not a good option to define in 'my_options' file (as it is a list) - define it here directly
  o$known_params = NULL  # c("contacts.BE" = 15, "npi_scaler.BE" = 0.8)
  
  # ---- Likelihood settings ----
  
  # Flag for recalculating likelihood for existing simulations
  #
  # NOTE: This can be used to test different weighting settings for existing simulations
  o$recalculate_fit = FALSE
  
  # Proportion of calibration jobs we can accept failing
  o$failure_tolerance = 0.05
  
  # Weight more heavily to recent data points - value between 0 and 1
  o$recent_data_weight = 0.5  # Applies a linear weighting from first point to last
  
  # Ignore the first n days of data in calibration (from first confirmed case)
  o$ignore_data_days = 0
  
  # Load calibration multiplers (see load_input.R)
  o = format_calibration_multipliers(o)
  
  # ---- Uncertainty settings ----
  
  # Which parameter set to use for 'best estimate' projection
  #
  # OPTIONS:
  #  "simulated" := The best of all actually simulated parameter sets
  #  "emulated"  := The best parameter set identified by MCMC of emulated space
  o$best_parameter_set = "simulated"
  
  # How should the 'best' simulated parameter set be quantified
  #
  # OPTIONS:
  #  "mean"    := The best mean likelihood across all seeds - any sets with NAs are discarded
  #  "mean_na" := As above but NAs are just ignored and can still be selected
  #  "max"     := The very best seed simulation (not so robust but occasionally produces better results)
  o$best_parameter_by = "mean_na"
  
  # Which parameter set to use for 'best estimate' projection
  #
  # OPTIONS:
  #  "median"    := Median of uncertainty simulations (stochastic and parameter uncertainty)
  #  "mean"      := Mean of uncertainty simulations (stochastic and parameter uncertainty)
  #  "simulated" := Directly use best fitting simulation as best estimate (not ideal)
  o$best_estimate_simulation = "mean"  # TODO: Implement 'simulated' if needed
  
  # Number of parameters sets to sample from MCMC posteriors
  o$n_parameter_samples = 0  # Not including best parameter set
  
  # Number of seeds to run for each analysis scenario and/or strategy (including baseline)
  o$n_seeds_analysis = 10  # Not to be confused with 'n_seeds_emulator'
  
  # Quantiles for credibility intervals
  o$quantiles = c(0.1, 0.9)  # c(0.025, 0.975)
  
  # ---- Cluster settings ----
  
  # Flag for overwriting any existing calibration simulations
  o$overwrite_fit = FALSE
  
  # Flag for overwriting any existing analysis simulations
  o$overwrite_analysis = TRUE
  
  # Choose cluster partition for parallel jobs
  o$cluster_partition = "covid19" # OPTIONS: "covid19" (best option) or "scicore"
  
  # Set thresholds (in terms of max_population_size) for when to use 'med' and 'big' job options for the cluster
  #
  # NOTE: Any value of max_population_size ABOVE these will use 'job size' options
  o$job_size_thresholds = c("small_job" = 0, "med_job" = 5e5, "big_job" = 9e5) 
  
  # If running on scicore partition, you also need to define a job time 
  # 
  # NOTES:
  #  1) This should be an OVER-ESTIMATE of time needed per job - this is important as jobs over this will fail
  #  2) For more details, see: https://wiki.biozentrum.unibas.ch/display/scicore/4.+Queues+and+partitions
  o$job_time  = qc(small_job = "00:30:00", med_job = "02:00:00", big_job = "02:00:00")  # Use "HH:MM:SS" or "D-HH:MM:SS" format
  o$job_queue = qc(small_job = "30min",    med_job = "6hours",   big_job = "6hours")  # Queue to use (check out scicore wiki if unfamiliar)
  
  # Memory to allocate for each cluster job
  #
  # NOTES: 
  #  1) General principle - the higher this is set, the less resources you get
  #  2) Anything upto (and incuding) 3.7GB will give maximum resources (~500 cores at a time)
  #  3) For populations sizes of over 1m, you'll likely want something around 8GB
  o$job_memory = qc(small_job = "3700M", med_job = "3700M", big_job = "8GB")
  
  # Some jobs require some extra juice, force a higher memory requirement
  #
  # NOTE: Currently used for summarising analysis outputs - may be a large number of simulations
  o$force_memory = "16GB"
  
  # Set an upper limit for jobs that can be run at any one time
  o$job_limit = 600
  
  # Define name of log and error files
  o$log_file = "scicore_log.txt"
  o$err_file = "scicore_error.txt"
  
  # Action to take if user is already running cluster jobs
  o$cluster_conflict_action = "error"
  
  # ---- Output settings ----
  
  # Flag for producing csv outputs (see results.R for what outputs are to be produced)
  o$output_csv = FALSE
  
  # Lower bound of age groups for plotting - bounded above by max(o$all_ages)
  o$plot_ages = c(0, 18, 50, 65, 75, 85)
  
  # Count up to and including n infections per person
  o$max_infection_count = 3
  
  # ---- Plotting settings ----
  
  # TODO: Decide which metrics details are inherited...
  
  # Load model metrics to collate and plot (see load_input.R)
  o = format_model_metrics(o)
  
  # Default metrics to plot if otherwise undefined	
  o$default_plot_metrics = c("confirmed", "deaths", "hospital_beds", "icu_beds")
  
  # Colour packages and palettes for scenarios, metrics, and cantons (see colour_scheme in myRfunctions.R)
  o$palette_scenario = "brewer::dark2"
  o$palette_metric   = "base::rainbow"
  o$palette_canton   = "viridis::viridis"
  
  # Colour packages and palettes for groupings (see colour_scheme in myRfunctions.R)
  o$palette_age     = "brewer::set2"
  o$palette_variant = "brewer::set1"
  o$palette_vaccine_priority = "brewer::set1"
  
  # Define some nice properties for baseline metric plots
  o$baseline_name   = "Baseline scenario"
  o$baseline_colour = "grey50" # ""#A4A4A4"  # Light grey
  
  # Grey colour for current date dashed line
  o$data_colour = "#555555"  # Dark grey
  o$dash_colour = "#808080"  # Even darker grey
  
  # Parameter colour scheme: prior, posterior, and global
  o$prior_colour     = "darkred"
  o$posterior_colour = "blue1"
  o$global_colour    = "blue4"
  
  # Space between x tick marks
  o$x_tick_dates = "1 month" # "2 weeks"
  o$x_tick_format = "%b-%y" # Month-year: "%b %y"; Day-month: "%d %b"
  
  # Saved figure size
  o$save_width  = 14
  o$save_height = 10
  
  # ---- Plotting flags ----
  
  # Turn figures on or off
  o$plot_baseline     = FALSE  # Standard baseline calibration plots
  o$plot_future       = FALSE  # Future projection plots with any defined future scenarios
  o$plot_past         = FALSE  # Past evaluation plots with any defined past scenarios
  o$plot_strategies   = FALSE  # Plot results of any defined strategies
  o$plot_elements     = FALSE  # Plot full factorial of strategy elements (requires plot_strategies)
  o$plot_cumulative   = FALSE  # Cumulative outcome figures
  o$plot_data_sources = FALSE  # Visualise the different data sources
  o$plot_emulator     = FALSE  # Plot emulator performance
  o$plot_parameters   = FALSE  # Parameter posterior plots
  
  # Manuscript plotting flags 
  o$plot_fig2 = TRUE
  o$plot_fig3 = TRUE
  o$plot_fig4 = TRUE
  o$plot_fig5 = TRUE
  o$plot_fig6 = TRUE
  o$plot_fig7 = TRUE
  
  # Supplement plotting flags 
  o$plot_sup1 = TRUE
  
  # ---- Advanced functionality ----
  
  # Override options set in my_options file
  #
  # NOTE: Cantons already dealt with in load_input.R
  o = override_options(o, ignore = "cantons", quiet = quiet)
  
  # Convert intermediary options to functional options
  o = functional_options(o)
  
  # Sanity checks for particular options
  options_checks(o)
  
  # ---- Display analysis details ----
  
  if (!quiet) {
    
    # Display both analysis and calibration names
    message(" - Calibration name: ", o$calibration_name)
    message(" - Analysis name: ",    o$analysis_name)
    
    # Display list of cantons we're modelling
    message(" - Canton(s): ", paste(o$cantons, collapse = ', '))
  }
  
  return(o)
}

# ---------------------------------------------------------
# Override options set in my_options file
# ---------------------------------------------------------
override_options = function(o, subset = NULL, ignore = NULL, quiet = FALSE) {
  
  # If user has a 'my options' file, load it
  # 
  # NOTE: Only the first sheet (in terms of position) is considered - sheet names are not important
  if (file.exists(o$pth$my_options)) {
    my_options = read.csv(o$pth$my_options)
    
    # Subset options to override if necessary
    if (!is.null(subset))
      my_options = filter(my_options, option %in% subset)
    
    # Remove any options we do not want to override
    if (!is.null(ignore))
      my_options = filter(my_options, !option %in% ignore)
    
    # Continue if we have options to override
    if (nrow(my_options) > 0) {
      
      if (!quiet) message(" - Overriding options using config/my_options.csv file")
      
      # Throw an error if there are entries that are not well defined options
      unrecognised = my_options$option[!my_options$option %in% names(o)]
      if (length(unrecognised) > 0)
        stop("Unrecognised entries in 'my options' file: ", paste(unrecognised, collapse = ", "))
      
      # Throw an error if there are multiple entries for any one option
      duplicates = my_options$option[duplicated(my_options$option)]
      if (length(duplicates) > 0)
        stop("Duplicate entries in 'my options' file: ", paste(duplicates, collapse = ", "))
      
      # Variable class of the options we wish to overwrite
      class_conversion = paste0("as.", lapply(o[my_options$option], class))
      
      # Iterate through the options and overwrite with value of correct class
      for (i in seq_len(nrow(my_options))) {
        
        # Format entry - primary job is to seperate out multiple values
        format_value = base::trimws(str_split(my_options$value[i], ",")[[1]])
        
        # Overwrite the converted option as defined in 'my options' file
        o[[my_options$option[i]]] = get(class_conversion[i])(format_value)
      }
    }
  }
  
  return(o)
}

# ---------------------------------------------------------
# Convert intermediary options to functional options
# ---------------------------------------------------------
functional_options = function(o) {
  
  # This is done here as user may have overwritten key options through 'my options' file
  
  # ---- Dates ----
  
  # Set initial date and extend (any data prior to this will be ignored)
  dates_data   = seq(format_date("2020-01-01"), format_date(o$data_end_date), by = "day")
  dates_future = seq(max(dates_data) + 1, max(dates_data) + o$n_future_days, by = "day")
  
  # Concatenate all dates
  o$dates_data = dates_data
  o$dates_all  = c(dates_data, dates_future)
  o$n_dates    = length(o$dates_all)
  
  # Remove intermediary options
  o[c("data_end_date", "n_future_days")] = NULL
  
  # ---- Cantons ----
  
  # Load Swiss canton dictionary from file
  canton_dict_df = read.csv(o$pth$cantons) %>%
    filter(code %in% o$cantons)
  
  # Format into named vector for accessibility
  o$canton_dict = setNames(canton_dict_df$name, canton_dict_df$code)
  
  # ---- Scenarios and strategies ----
  
  # Covert scenario vectors to list and remove trivial
  scenarios   = list(past = o$scenarios_past, future = o$scenarios_future)
  o$scenarios = lapply(scenarios, function(x) x[x != "none"])
  
  # Remove intermediary options
  o[c("scenarios_past", "scenarios_future")] = NULL
  
  # ---- Model outputs ----
  
  # We'll count all ages - we can then aggregate however we like
  o$count_age = o$all_ages
  
  # Check we have a valid value of max_infection_count
  if (!isTRUE(o$max_infection_count >= 1))
    stop("Option 'max_infection_count' must be non-negative")
  
  # We'll stop the count when a person has this many infections
  o$count_infections = 1 : round(o$max_infection_count)
  
  # Variants is more straightforward - count them all
  o$count_variant = names(o$variants_dict)
  o$count_vaccine_priority = c("none", names(o$vaccine_priority))
  
  # Remove intermediary options
  o[c("max_infection_count")] = NULL
  
  # ---- Cluster specifications ----
  
  # Size of job we are running (small, medium, or big)
  job_size_idx = max(which(o$max_population_size > o$job_size_thresholds))
  job_size = names(o$job_size_thresholds[job_size_idx])
  
  # Use options defined for this job size
  o$job_memory = o$job_memory[[job_size]]
  o$job_time   = o$job_time[[job_size]]
  o$job_queue  = o$job_queue[[job_size]]
  
  # Overwrite queue if using covid19 partition
  if (o$cluster_partition == "covid19")
    o$job_queue = "covid19"
  
  return(o)
}

# ---------------------------------------------------------
# Sanity checks for particular options
# ---------------------------------------------------------
options_checks = function(o) {
  
  # If running national level, check individual cantons have not been defined
  if ("CH" %in% o$cantons && length(o$cantons) > 1)
    stop("If running 'CH' you cannot select other cantons - see o$cantons")
  
  # If running national level, check individual cantons have not been defined
  if (!identical(o$cantons, unique(o$cantons)))
    stop("Not unique cantons defined - see o$cantons")
  
  # TODO: Add other checks here for things such as dates, probabilities, etc
}

# ---------------------------------------------------------
# Store a subset of options after calibrating & simulating
# ---------------------------------------------------------
store_options = function(o, step) {
  
  # Options we always want to store - all related to model structure
  always_store = qc(model_flows, 
                    disease_states, 
                    care_states, 
                    prognosis_states, 
                    all_states, 
                    age_groups, 
                    all_ages,
                    network_structure, 
                    contact_matrix_countries, 
                    max_population_size, 
                    replace_dead, 
                    viral_load_infectivity, 
                    variants, 
                    vaccine_efficacy, 
                    vaccine_default, 
                    vaccine_priority, 
                    vaccine_criteria)
  
  # TODO: Steps "emulator" and "calibration" should be seperated...
  
  # Other options of interest when storing calibration files
  if (step == "calibration")
    also_store = qc(calibration_name, 
                    dates_data, 
                    data_source, 
                    set_negative_data, 
                    diagnosis_ma_days, 
                    weather_ma_days, 
                    calibration_df, 
                    fixed_df, 
                    calibration_multiplier, 
                    emulator_reproducible, 
                    n_init_samples, 
                    n_seeds_emulator, 
                    train_test_split, 
                    gp_kernel, 
                    gp_max_iter, 
                    resample_iters, 
                    acquisition_funcion, 
                    acquisition_points, 
                    min_value_accepted, 
                    max_samples_accepted, 
                    distance_method, 
                    distance_tolerance, 
                    calibration_reproducible, 
                    mcmc_n_chains, 
                    mcmc_iters, 
                    mcmc_opt_freq, 
                    mcmc_opt_stages, 
                    mcmc_init_steps, 
                    # known_params, 
                    recent_data_weight, 
                    prior2data_weight, 
                    hyperprior_weight, 
                    ignore_data_days)
  
  # Other options of interest when storing analyses files
  if (step == "analysis")
    also_store = qc(analysis_name, 
                    dates_all, 
                    n_dates, 
                    best_parameter_set, 
                    best_parameter_by, 
                    n_parameter_samples, 
                    n_seeds_analysis, 
                    quantiles)
  
  # Subset the o list which is to be stored
  o_subset = o[c("user", also_store, always_store)]
  
  return(o_subset) 
}

