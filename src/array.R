###########################################################
# ARRAY
#
# Everything related to array gird and LHC scenarios in one place. 
# Functions for identifying, parsing, fitting, and summarising.
#
###########################################################

# ---------------------------------------------------------
# Parent function for operating on array scenarios
# ---------------------------------------------------------
run_arrays = function(o) {
  
  # Only continue if specified by do_step (or forced)
  if (!is.element(3, o$do_step)) return()
  
  message("* Operating on array scenarios")
  
  # ---- Identify array grid and LHC parents ----
  
  # All individual scenarios
  all_scenarios = names(parse_yaml(o, "*read*", read_array = TRUE))
  
  # Indentify any parents of these individual scenarios
  all_parents = get_array_parents(all_scenarios)
  
  # Report if no array scenarios identified
  if (length(all_parents) == 0)
    message(" - No array scenarios identified")
  
  # Determine if any parents are LHC scenarisos
  lhc_vec = is_parent_lhc(o, all_parents)
  
  # Throw an error if type of array could not be established
  unknown_type = all_parents[is.na(lhc_vec)]
  if (length(unknown_type) > 0)
    stop("Array type could not be established: ", paste(unknown_type, collapse = ", "))
  
  # Remaining parents are non-LHC array scenarios
  parents_lhc  = all_parents[lhc_vec]
  parents_grid = all_parents %>% setdiff(parents_lhc)
  
  # ---- Operate on these parents ----
  
  # Train a ML model on parameter inputs and chosen outputs
  predict_lhc_scenarios(o, all_scenarios, parents_lhc)
  
  # Squish outputs of all array children into one datatable
  #
  # NOTE: This makes it much quicker to plot heatmaps
  summarise_grid_scenarios(o, all_scenarios, parents_grid)
}

# ---------------------------------------------------------
# Train a ML model for all user-defined endpoints
# ---------------------------------------------------------
predict_lhc_scenarios = function(o, scenarios, parents) {
  
  # TODO: Do not repeat for endpoints that have previously been done (unless forced)
  
  # Parse user-defined endpoint and emulator options
  input = parse_yaml(o, "baseline")$parsed
  
  # Extract list of endpoint and associated details
  endpoints = input$scenario_lhc_endpoints
  
  # Loop through all distinct parent array scenarios
  for (parent in parents) {
    
    message(" - Array LHC scenario: ", parent)
    
    # ---- Load inputs ----
    
    message("  > Loading inputs")
    
    # Load LHC array scenario info 
    array_info = try_load(o$pth$array_info, parent)
    
    # Extract parameter sets simulated
    input_df = array_info$values %>%
      pivot_wider(id_cols    = scenario, 
                  names_from = variable_id) %>%
      setDT()
    
    # ---- Load outputs ----
    
    message("  > Loading outputs")
    
    # All children scenarios associated with this parent
    array_children = get_array_children(scenarios, parent)
    
    # Full path to all raw output files
    output_files = paste0(o$pth$scenarios, array_children, "_raw.rds")
    
    # Load all files into one long datatable
    output_df = rbindlist(lapply(output_files, readRDS))
    
    # ---- Training input -> output ML model ----
    
    message("  > Training predictor")
    
    # Loop through user-defined endpoints
    for (i in seq_along(endpoints)) {
      e = endpoints[[i]]
      
      # Extract endpoint of interest and join with parameter sets
      endpoint_df = output_df %>%
        filter(metric == !!e$metric) %>%
        group_by(scenario, seed) %>%
        summarise(value = get(e$summarise)(value)) %>%
        ungroup() %>%
        full_join(input_df, by = "scenario") %>%
        select(scenario, all_of(array_info$vars$variable_id), value) %>%
        setDT()
      
      # Train input -> endpoint predictor
      predictor = train_predictor(o, endpoint_df, array_info, input$emulator)
      
      # Save predictor for this endpoint to file
      saveRDS(predictor, paste0(o$pth$endpoints, parent, ".", e$id, ".rds"))
    }
    
    # ---- Visualise the predictions ----
    
    message("  > Plotting diagnostics")
    
    # Plot LHC scenario diagnostics
    plot_lhc_diagnostics(o, "LHC diagnostics", parent)
  }
}

# ---------------------------------------------------------
# Squish outputs of all array children into one datatable
# ---------------------------------------------------------
summarise_grid_scenarios = function(o, scenarios, parents) {
  
  # Loop through all distinct parent array scenarios
  for (parent in parents) {
    
    message(" - Array grid scenario: ", parent)
    
    # All children scenarios associated with this parent
    array_children = get_array_children(scenarios, parent)
    
    # TODO: Load all grid array children and combine into single df
    #       The idea being to speed up heatmap plotting in step 4
    
    # Save summarised grid in: o$pth$parents
  }
  
}

# ---------------------------------------------------------
# Flatten and parse parent array scenarios into their children
# ---------------------------------------------------------
parse_array_grid = function(o, y, scenario) {
  
  # Full names of all scenarios
  scenarios_name = unlist(lapply(y$scenarios, function(x) x$name))
  
  # ---- Identify array scenarios ----
  
  # Indentify where (if anywhere) arrays are defined
  array_idx   = grepl("*\\.array_grid\\.*", names(unlist(y$scenarios)))
  array_items = names(unlist(y$scenarios))[array_idx]
  
  # Names of parent scenarios that have array items
  parents = get_array_parents(array_items)
  
  # Loop through all parent array scenarios
  for (parent in parents) {
    
    # Array variables associated with this parent
    array_vars_idx = grepl(paste0("^", parent, "\\."), array_items)
    array_vars = array_items[array_vars_idx] %>%
      str_remove(paste0(parent, "\\.")) %>%
      str_remove("\\.array_grid\\..*") %>%
      str_replace_all("\\.", "$") %>%
      unique()
    
    # Loop through these variables
    vars_list = list()
    for (var in array_vars) {
      
      # String which - when evaluated - indexes the deep nested array details for this variable 
      array_str = paste("y$scenarios", parent, var, "array_grid", sep = "$")
      
      # Evaluate the string to obtain a list of inputs to R's seq function
      array_eval = eval_str(array_str)
      
      # Remove the ID and name items, and set the function name to seq
      array_fn = c(fn = "seq", list.remove(array_eval, c("id", "name")))
      
      # Evaluate the seq function
      all_vals = parse_fn(array_fn)
      
      # Use child scenario indices for IDs
      all_ids = paste0(array_eval$id, 1 : length(all_vals))
      
      # Store datatable of all key details
      var_df = data.table(dim  = array_eval$id, 
                          id   = all_ids, 
                          var  = var, 
                          val  = all_vals, 
                          desc = array_eval$name)
      
      # Store naming convention (if it has been provided)
      if (!is.null(array_eval$name))
        var_df$name = paste0(array_eval$name, ": ", all_vals)
      
      # Store datatable
      vars_list[[var]] = var_df
    }
    
    # Datatable of all values of each variable
    vars_df = rbindlist(vars_list, fill = TRUE) %>%
      mutate(across(c(dim, id), fct_inorder))
    
    # Expand grid to all individual array elements
    array_df = vars_list %>%
      lapply(function(x) x$id) %>%
      unique() %>%
      expand.grid() %>%
      unite("x", sep = ".", remove = FALSE) %>%
      pivot_longer(-x, names_to = "var_ref", 
                   values_to    = "id") %>%
      left_join(vars_df, by = "id") %>%
      pivot_wider(id_cols = x, 
                  names_from  = var, 
                  values_from = val) %>%
      mutate(id = paste0(parent, ".", x)) %>%
      select(id, all_of(array_vars)) %>%
      setDT()
    
    # Extract full names of all array scenarios
    array_names = vars_df %>%
      group_by(id) %>%
      slice_head(n = 1) %>%
      split(., f = .$dim) %>%
      lapply(function(x) x$name) %>%
      expand.grid() %>%
      unite("names", sep = ", ") %>%
      pull(names)
    
    # Ensure we have equal numbers
    #
    # NOTEL: This error will also be thrown if dim parts are of different sizes
    if (nrow(array_df) != length(array_names))
      stop("Inconsistent number of array scenario elements")
    
    # Apply names of elements to array datatable
    array_df = array_df %>%
      mutate(name = paste0(scenarios_name[[parent]], 
                           " (", array_names, ")"), 
             parent = parent) %>%
      select(id, parent, all_of(array_vars), name)
    
    # Construct individual child scenarios and remove parent
    y = create_array_children(y, parent, array_df)
    
    # ---- Store array details ----
    
    # Only store array details if not reading scenario names
    if (scenario == "*create*") {
      
      # Initiate array info list with basic scenario ID dataframe
      array_info = list(type = "grid", 
                        meta = array_df[, .(parent, scenario = id, name)])
      
      # Append info about the variables - IDs and names
      array_info$vars = vars_df[, .(dim, desc)] %>%
        rename(variable_id   = dim, 
               variable_name = desc) %>%
        group_by(variable_id) %>%
        slice_head(n = 1) %>%
        ungroup() %>%
        setDT()
      
      # Dimension details
      dim_df = vars_df[!duplicated(vars_df$dim), .(dim, var)]
      
      # Store all the values to be simulated
      array_info$values = array_df %>%
        rename(setNames(dim_df$var, dim_df$dim)) %>%
        select(scenario = id, all_of(dim_df$dim)) %>%
        pivot_longer(cols = -scenario, 
                     names_to = "variable_id") %>%
        setDT()
      
      # Save this info for use when collating all results
      saveRDS(array_info, file = paste0(o$pth$array_info, parent, ".rds"))
    }
  }
  
  return(y)
}

# ---------------------------------------------------------
# Generate LHC samples from full set of parameter bounds
# ---------------------------------------------------------
parse_array_lhc = function(o, y, scenario) {
  
  # Full names of all scenarios
  scenarios_name = unlist(lapply(y$scenarios, function(x) x$name))
  
  # ---- Identify LHC scenarios ----
  
  # Indentify where (if anywhere) arrays are defined
  array_idx   = grepl("*\\.array_lhc\\.*", names(unlist(y$scenarios)))
  array_items = names(unlist(y$scenarios))[array_idx]
  
  # Names of parent scenarios that have array items
  parents = get_array_parents(array_items)
  
  # Loop through all parent array scenarios
  for (parent in parents) {
    
    # Array variables associated with this parent
    array_vars_idx = grepl(paste0("^", parent, "\\."), array_items)
    array_vars = array_items[array_vars_idx] %>%
      str_remove(paste0(parent, "\\.")) %>%
      str_remove("\\.array_lhc\\..*") %>%
      str_replace_all("\\.", "$") %>%
      unique()
    
    # Number of array variables
    n_vars = length(array_vars)
    
    # Preallocate matrix for parameter bounds (required by LHC sampler)
    lhc_bounds = matrix(nrow = n_vars, ncol = 2)
    
    # Loop through these variables
    vars_list = list()
    for (i in 1 : n_vars) {
      
      # String which - when evaluated - indexes the deep nested array details for this variable 
      array_str = paste("y$scenarios", parent, array_vars[i], "array_lhc", sep = "$")
      
      # Evaluate the string to obtain bounds to be used in LHC sampler
      array_eval = eval_str(array_str)
      
      # Most important, extract the bounds for this parameter
      lhc_bounds[i, ] = c(array_eval$lb, array_eval$ub)
      
      # Also store ID and name of array element
      vars_list[[i]] = as.data.table(array_eval[qc(id, name)])
    }
    
    # ---- Create or load LHC samples ----
    
    # Shorthand for number of LHC samples to generate
    n = y$scenario_lhc_samples
    
    # Check whether we need to create the samples
    if (scenario == "*create*") {
      
      # If we want this to be reproducible, start a new seed (unique for each parent)
      if (o$scenario_reproducible == TRUE)
        set.seed(letterToNumber(gsub("_", "", parent)))
      
      # Apply LHC sample to generate points across the space
      lhc_samples = tgp::lhs(n, lhc_bounds)
    }
    
    # If only reading scenario names, we can set trivial samples
    if (scenario == "*read*") 
      lhc_samples = matrix(nrow = n, ncol = n_vars)
    
    # Otherwise we'll need to load previously created samples
    if (!scenario %in% c("*read*", "*create*")) {
      
      # Load samples and convert to wide matrix
      lhc_samples = try_load(o$pth$array_info, parent)$values %>%
        pivot_wider(id_cols    = "scenario", 
                    names_from = "variable_id") %>%
        select(-scenario) %>%
        as.matrix()
    }
    
    # Construct scenario IDs and names (simplier than arrays)
    lhc_id    = paste(parent, seq_len(n), sep = ".")
    lhc_names = paste(scenarios_name[parent], seq_len(n))
    
    # Construct datatable with created or loaded samples 
    lhc_df = data.table(parent = parent, 
                        id     = lhc_id, 
                        name   = lhc_names) %>%
      cbind(as_named_dt(lhc_samples, array_vars))
    
    # Construct individual child scenarios and remove parent
    y = create_array_children(y, parent, lhc_df)
    
    # ---- Store LHC details ----
    
    # Only store array details if not reading scenario names
    if (scenario == "*create*") {
      
      # Initiate array info list with basic scenario ID dataframe
      array_info = list(type = "lhc", 
                        meta = lhc_df[, .(parent, scenario = id, name)])
      
      # Append info about the variables - IDs and names
      array_info$vars = rbindlist(vars_list) %>%
        rename(variable_id   = id, 
               variable_name = name)
      
      # Set row names of LHC parameter bounds matrix
      rownames(lhc_bounds) = array_info$vars$variable_id
      
      # Append this updated bounds matrix
      array_info$bounds = lhc_bounds
      
      # Store all the values to be simulated
      array_info$values = lhc_df %>%
        select(scenario = id, all_of(array_vars)) %>%
        rename(setNames(array_vars, array_info$vars$variable_id)) %>%
        pivot_longer(cols = -scenario, 
                     names_to = "variable_id") %>%
        setDT()
      
      # Save this info for use when collating all results
      saveRDS(array_info, file = paste0(o$pth$array_info, parent, ".rds"))
    }
  }
  
  return(y)
}

# ---------------------------------------------------------
# Construct child scenarios from array details
# ---------------------------------------------------------
create_array_children = function(y, parent, array_df) {
  
  # Store details of the parent array scenario
  array_details = y$scenarios[[parent]]
  
  # Loop through each individual child
  for (i in 1 : nrow(array_df)) {
    this_scen = array_df[i, ]
    
    # Start with the basic structure of the parent
    scen_details = list_modify(array_details, 
                               id   = this_scen$id, 
                               name = this_scen$name)
    
    # Reduce array_df down to just values we want to use for this child scenario
    val_df = this_scen %>% select(-id, -name, -parent)
    
    # Loop through the values
    for (val_name in names(val_df)) {
      this_val = val_df[[val_name]]
      
      # Class of parameter in baseline - could be numeric or integer
      param_class = class(eval_str("y$", val_name))
      
      # If integer, ensure consistent class
      if (param_class == "integer")
        this_val = paste0("as.integer(", this_val, ")")
      
      # String to be evaluated to set the individual value
      eval_str("scen_details$", val_name, " = ", this_val)
    }
    
    # Append this child scenario
    y$scenarios[[this_scen$id]] = scen_details
  }
  
  # Remove the parent scenario
  y$scenarios[[parent]] = NULL
  
  return(y)
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
    parent_exp = paste0("^", parent, "\\..*")
    
    # Indices of such scenarios
    children_idx = grepl(parent_exp, scenarios)
    
    # Extract scenario names
    children = c(children, scenarios[children_idx])
  }
  
  return(children)
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
# Test whether array parent is a LHC scenario
# ---------------------------------------------------------
is_parent_lhc = function(o, parents) {
  
  # Preallocate vector output
  parent_lhc = rep(NA, length(parents))
  
  # Loop through array parents
  for (i in seq_along(parents)) {
    
    # Construct useful error message in case file cannot be found
    err_msg = paste0("Scenario '", parents[i], "' does not appear to be an array parent")
    
    # Load array info - error thrown if this doesn't (yet) exist
    array_info = try_load(o$pth$array_info, parents[i], msg = err_msg)
    
    # Set true for LHC, false for non-LHC array, and leave as NA otherwise
    if (array_info$type == "lhc")  parent_lhc[i] = TRUE
    if (array_info$type == "grid") parent_lhc[i] = FALSE
  }
  
  return(parent_lhc)
}

