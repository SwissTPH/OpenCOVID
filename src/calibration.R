###########################################################
# CALIBRATION
#
# Fit parameters to achieve a given R_eff at the 'initial'
# point for forward projections.
#
# This is a basic implementation and will be expanded in 
# version 2.1.
#
###########################################################

# ---------------------------------------------------------
# Parent function to fit contacts to desired R_eff
# ---------------------------------------------------------
run_calibration = function(o) {
  
  # Only continue if specified by do_step
  if (!is.element(1, o$do_step)) return()
  
  message("* Calibrating model")
  
  # TODO: Do this for different seeds and different starting points
  
  # ---- Parameters to be fitted ----
  
  # Load parsed parameters from yaml file
  fit_input = parse_yaml(o, scenario = "baseline")$parsed
  
  # Extract parameters to calibrate
  fit_df = list2dt(fit_input$calibration) %>%
    mutate(lower = as.numeric(lower), 
           upper = as.numeric(upper))
  
  # Take initial points from yaml defaults
  x0 = unlist(fit_input[fit_df$param])
  
  # ---- Optimise parametrs for R_eff ----
  
  # Provide full o list as additonal argument
  obj_args = append(o, list(params_to_fit = fit_df$param, model_seed = 1))
  
  # Run ASD algorithm to determine optimal contacts
  fitted_model = asd(objective_fn, x0, 
                     args = obj_args, 
                     lb = fit_df$lower, 
                     ub = fit_df$upper,
                     max_iters  = o$fit_iters_max,
                     plot_iters = o$fit_iters_plot)
  
  # Convert optimal result to list
  fit_output = as.list(fitted_model$x)
  
  # Compile fitting input and output
  fit_result = list(fit_input  = fit_input, 
                    fit_output = fit_output)
  
  # Save fitted input and output to file
  saveRDS(fit_result, file = paste0(o$pth$fitting, "fit_result.rds"))
}

# ---------------------------------------------------------
# Objective function to minimise - modelled and desired R_eff difference
# ---------------------------------------------------------
objective_fn = function(x, args) {
  
  # Model parameters to alter when running fitting simulation
  fit_params = as.list(setNames(x, args$params_to_fit))
  fit_list = c(fit_params, list(n_days = max(args$fit_days)))
  
  # Simulate model
  result = model(args, "baseline", 
                 seed = args$model_seed, 
                 fit  = fit_list, 
                 verbose = "none")
  
  # Days to index R_eff (considers initial phase for which it cannot be calculated)
  r_eff_idx = args$fit_days + result$input$infectious_period
  
  # Extract mean R_eff over these time points
  r_eff_df = filter(result$output, date %in% r_eff_idx)
  r_eff    = mean(r_eff_df$value)
  
  # Squared difference with desired R_eff
  objective_val = (r_eff - result$input$r_eff) ^ 2
  
  # Output in list format
  return(list(y = objective_val))
}

# ---------------------------------------------------------
# Alter key parameters in yaml file for fitting simulation
# ---------------------------------------------------------
fit_yaml = function(o, yaml, fit_list) {
  
  # Only needed if fitting OR using fitted parameters
  if (!is.null(fit_list)) {
    
    # Check if we are fitting here - number of days will be specified
    if (!is.null(fit_list$n_days)) {
      
      # Throw error if R_eff is not being reported - we need this to calibrate to
      if (!"R_effective" %in% yaml$parsed$metrics$df$metric)
        stop("You must turn on reporting of 'R_effective' for calibration")
      
      # Only interested in collecting one metric when fitting: R_effective
      yaml$parsed$metrics$df = filter(yaml$parsed$metrics$df, metric == "R_effective")
      yaml$parsed$metrics$groupings = yaml$parsed$metrics$groupings["R_effective"]
      
      # To do this we'll need to simulate n_days past infectious_period
      fit_list$n_days = fit_list$n_days + yaml$parsed$infectious_period
      
      # Sanity check that n_days from yaml file is at least max of fit_days
      if (yaml$parsed$n_days < fit_list$n_days)
        stop("Model calibration requires n_days >= ", fit_list$n_days)
    }
    
    # Loop through parameters to alter
    for (param in names(fit_list)) {
      
      # Check parameter is valid - throw error if not
      if (!param %in% names(yaml$parsed))
        stop("Attempting to alter invalid parameter '", param, "' during model fitting")
      
      # Re-define paramter value
      yaml$parsed[[param]] = fit_list[[param]]
    }
  }
  
  return(yaml)
}

