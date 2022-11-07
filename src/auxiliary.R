###########################################################
# AUXILIARY FUNCTIONS
#
# A series of helpful R functions.
#
# Written by A.J.Shattock
###########################################################

# ---------------------------------------------------------
# Set as datatable and rename columns in one line
# ---------------------------------------------------------
as_named_dt = function(x, new_names) {
  
  # Convert to datatable
  dt = as.data.table(x)
  
  # Check new names are correct length
  old_names = names(dt)
  if (length(old_names) != length(new_names))
    stop("Inconsistent number of column names provided")
  
  # Set new column names
  named_dt = setnames(dt, old_names, new_names)
  
  return(named_dt)
}

# ---------------------------------------------------------
# Check if user is currently running any cluster jobs
# ---------------------------------------------------------
check_cluster_jobs = function(o, action = "error") {
  
  # Check number of running and pending jobs
  n_jobs = n_slurm_jobs(user = o$user)
  
  # If this is non-zero then action needed
  if (sum(unlist(n_jobs)) > 0) {
    
    # TODO: Extend this by offering additional options (eg warning or prompt)
    
    # Throw an error
    if (action != "none")
      stop("You currently have jobs submitted to the cluster, this may lead to unexpected results.")
  }
}

# ---------------------------------------------------------
# Clear the console
# ---------------------------------------------------------
clc = function() cat("\014")

# ---------------------------------------------------------
# Clear all figures
# ---------------------------------------------------------
clf = function() graphics.off()

# ---------------------------------------------------------
# Create colour scheme
# ---------------------------------------------------------
colour_scheme = function(map, pal = NULL, n = 1, ...) {
  
  # Has colour palette been defined
  if (is.null(pal)) {
    
    # That's ok as long as it's defined within the map argument
    if (!grepl("::", map))
      stop("Palette not defined - Use 'pal = my_pal' or 'map = my_map::my_pal'")
    
    # Seperate out the map and the palette
    pal = str_remove(map, ".*\\::")
    map = str_remove(map, "\\::.*")
  }
  
  # Initiate colours variable
  colours = NULL
  
  # Built in colour schemes
  if (map == "base")
    colours = get(pal)(n, ...)	
  
  # A load of colour maps from the pals package
  #
  # See: https://www.rdocumentation.org/packages/pals/versions/1.6
  if (map == "pals")
    colours = get(pal)(n, ...)
  
  # Stylish HCL-based colour maps
  #
  # See: https://colorspace.r-forge.r-project.org/articles/hcl_palettes.html
  if (grepl("_hcl$", map))
    colours = get(map)(palette = pal, n = n, ...)
  
  # Colour Brewer colour schemes
  if (map == "brewer")
    colours = brewer_pal(palette = first_cap(pal), ...)(n)
  
  # Viridis colour schemes
  if (map == "viridis")
    colours = viridis_pal(option = pal, ...)(n)
  
  # Throw an error if colours not yetr defined
  if (is.null(colours))
    stop("Colour map '", map, "' not recognised (supported: base, pals, hcl, brewer, viridis)")
  
  return(colours)
}

# ---------------------------------------------------------
# Create a log file (see function wait_for_jobs)
# ---------------------------------------------------------
create_bash_log = function(pth, log = NULL, err = NULL) {
  for (this_name in c(log, err)) {
    this_file = paste0(pth, this_name)
    if (file.exists(this_file)) file.remove(this_file)
    Sys.sleep(0.1)
    file.create(this_file)
    Sys.sleep(0.1)
  }
  return(paste0(pth, log))
}

# ---------------------------------------------------------
# Evaluate a string (in calling function environment) using eval
# ---------------------------------------------------------
eval_str = function(...)
  eval(parse(text = paste0(...)), envir = parent.frame(n = 1))

# ---------------------------------------------------------
# Biphasic exponential function
# ---------------------------------------------------------
exp_biphasic = function(x, peak, p, d1, d2, vmax, alpha, beta) {
  
  # Both exponential 'phases': short and long
  exp1 = exp(-x * log(2)/d1) * p
  exp2 = exp(-x * log(2)/d2) * (1-p)
  
  # Combine the phases
  bi_exp = peak * (exp1 + exp2)
  
  # Bound above and below
  bi_exp_scaled = vmax * (1 - 1 / (1 + (bi_exp / beta)^alpha))
  
  return(bi_exp_scaled)
}

# ---------------------------------------------------------
# Double exponential function
# ---------------------------------------------------------
exp_double = function(x, b, g, h, d) {
  y = (h * b * exp(-d * x)) / ((h - b) * exp(-g * x) + b)
  return(y)
}

# ---------------------------------------------------------
# Extract facets rows and columns from a ggplot object
# ---------------------------------------------------------
facet_dims = function(g) {
  g_layout = ggplot_build(g)$layout$layout
  n_rows = length(unique(g_layout$ROW))
  n_cols = length(unique(g_layout$COL))
  return(c(n_rows, n_cols))
}

# ---------------------------------------------------------
# A wrapper for tagger::tag_facets with a few extras
# ---------------------------------------------------------
facet_labels = function(g, ...) {
  
  # Default arguments
  args = list(tag = "panel", 
              tag_levels = c("A", "1"), 
              tag_suffix = "")
  
  # Overwrite these if desired
  args = list_modify(args, !!!list(...))
  
  # Call tag_facets function with these arguments
  g = g + do.call(tag_facets, args) 
  
  return(g)
}

# ---------------------------------------------------------
# Platform specific file separator - for readability
# ---------------------------------------------------------
file_sep = function() {
  platform_file_sep = .Platform$file.sep
  return(platform_file_sep)
}

# ---------------------------------------------------------
# Capitalise first letter of a string (or vector of strings)
# ---------------------------------------------------------
first_cap = function(string) {
  string_cap = paste0(toupper(substring(string, 1, 1)), substring(string, 2))
  return(string_cap)
}

# ---------------------------------------------------------
# Format heterogeneous styles of dates
# ---------------------------------------------------------
format_date = function(dates, convert = "ymd") {
  styles = c("dmy", "dmY", "ymd", "Ymd")
  dates = parse_date_time(dates, styles)
  dates = get(convert)(dates)
  return(dates)
}

# ---------------------------------------------------------
# Convert list to datatable
# ---------------------------------------------------------
list2dt = function(x, ...) {
  dt = rbindlist(lapply(x, as.data.table), ...)
  return(dt)
}

# ---------------------------------------------------------
# Logistic function
# ---------------------------------------------------------
logistic = function(x, slope, mid, lower = 0, upper = 1) {
  y = upper + (lower - upper) / (1 + (x / mid) ^ slope)
  return(y)
}

# ---------------------------------------------------------
# Inverse logistic function
# ---------------------------------------------------------
logistic_inv = function(x, slope, mid, lb, ub) {
  y = 1 - logistic(x, slope, mid, lower = 1 - ub, upper = 1 - lb)
  return(y)
}

# ---------------------------------------------------------
# Parameter transformation: put a probability on the real number line
# ---------------------------------------------------------
logit = function(p) {
  z = log(p / (1 - p))
  return(z)
}

# ---------------------------------------------------------
# Parameter transformation: inverse of the above
# ---------------------------------------------------------
logit_inv = function(p_logit) {
  z = exp(p_logit) / (exp(p_logit) + 1)
  return(z)
}

# ---------------------------------------------------------
# Number of running and pending jobs on the cluster
# ---------------------------------------------------------
n_slurm_jobs = function(user) {
  
  # Base sq command for user
  sq = paste("squeue -u", user)
  
  # Concatenate full commands
  slurm_running  = paste(sq, "-t running | wc -l")
  slurm_pending  = paste(sq, "-t pending | wc -l")
  slurm_ondemand = paste(sq, "-q interactive | wc -l")  # Interactive jobs
  
  # Function to get number of jobs (minus 1 to remove header row)
  get_jobs_fn = function(x) as.numeric(system(x, intern = TRUE)) - 1
  
  # System call to determine number of slurm processes
  n_running  = get_jobs_fn(slurm_running)
  n_pending  = get_jobs_fn(slurm_pending)
  n_ondemand = get_jobs_fn(slurm_ondemand)  # Interactive jobs
  
  # Compile into list
  #
  # NOTE: ondemand jobs are considered 'running', so discount these
  n_jobs = list(running = n_running - n_ondemand, pending = n_pending)
  
  return(n_jobs)
}

# ---------------------------------------------------------
# Normalise a vector of values to between 0 and 1
# ---------------------------------------------------------
normalise_0to1 = function(x, x_min = NULL, x_max = NULL, direction = "forward") {
  
  if (!tolower(direction) %in% c("forward", "backward"))
    stop("Normalisation direction must be either 'forward' or 'backward'")
  
  # Forward normalisation
  if (tolower(direction) == "forward") {
    
    # Take bounds from data unless given
    if (is.null(x_min)) x_min = min(x)
    if (is.null(x_max)) x_max = max(x)
      
    # Normalisation equation
    y = (x - x_min) / (x_max - x_min)
    
    # Append original min and max values, needed to backtransform
    attributes(y)$x_min = x_min
    attributes(y)$x_max = x_max
  }
  
  # Backward normalisation
  if (tolower(direction) == "backward") {
    
    # Take bounds from attitubutes of pre-normalised data unless given
    if (is.null(x_min)) x_min = attributes(x)$x_min
    if (is.null(x_max)) x_max = attributes(x)$x_max
    
    # Rearrange equation to solve for x
    #
    # NOTE: as.vector removes all attributes
    y = as.vector(x * (x_max - x_min) + x_min)
  }
  
  return(y)
}

# ---------------------------------------------------------
# Sub a directory name within a file path string
# ---------------------------------------------------------
pth_replace = function(pth, dir, dir_replace, sep = file_sep()) {
  new_pth = gsub(paste0(sep, dir, sep), 
                 paste0(sep, dir_replace, sep), pth)
  return(new_pth)
}

# ---------------------------------------------------------
# Suppress output from a function call
# ---------------------------------------------------------
quiet = function(x) { 
  sink_con = file("sink.txt")
  sink(sink_con, type = "output")
  sink(sink_con, type = "message")
  on.exit(sink(type   = "output"))
  on.exit(sink(type   = "message"), add = TRUE)
  on.exit(file.remove("sink.txt"),  add = TRUE)
  invisible(force(x)) 
}

# ---------------------------------------------------------
# Wrapper for consistent behaviour of base::sample when length(x) is one
# ---------------------------------------------------------
sample_vec = function(x, ...) x[sample.int(length(x), ...)]

# ---------------------------------------------------------
# Allow scale_fill_fermenter to accept custom palettes
# ---------------------------------------------------------
scale_fill_fermenter_custom = function(cols, guide = "coloursteps", na.value = "grey50", ...) {
  palette = ggplot2:::binned_pal(scales::manual_pal(cols))
  g = binned_scale("fill", "fermenter", palette, guide = guide, na.value = na.value, ...)
  return(g)
}

# ---------------------------------------------------------
# Initiate progress bar with normal-use options
# ---------------------------------------------------------
start_progress_bar = function(n_tasks) {
  pb = txtProgressBar(min = 0, max = n_tasks,
                      initial = 0, width = 100, style = 3)
  return(pb)
}

# ---------------------------------------------------------
# Check for cluster errors, stop if any found
# ---------------------------------------------------------
stop_if_errors = function(pth, err, err_tol = 0, msg = NULL) {
  
  # Set default error message
  if (is.null(msg))
    msg = "Fatal errors when running cluster jobs"
  
  # Check the error file exists
  err_file = paste0(pth, err)
  if (file.exists(err_file)) {
    
    # Load errors from file and convert to vector
    errors = readLines(err_file)
    
    # Stop if any errors, and display them all
    if (length(errors) >= (err_tol + 1))
      stop(msg, " (", length(errors), " errors)\n\n", 
           paste(errors, collapse = "\n"))
  }
}

# ---------------------------------------------------------
# Convert comma-seperated string to vector of elements
# ---------------------------------------------------------
str2vec = function(x, v) {
  x[[v]] = x[[v]] %>% 
    str_split(",", simplify = TRUE) %>% 
    str_remove_all(" ")
  return(x)
}

# ---------------------------------------------------------
# Submit jobs to the cluster and wait until they are complete
# ---------------------------------------------------------
submit_cluster_jobs = function(o, n_jobs, bash_file, ...) {
  
  # Skip if no jobs to run
  if (n_jobs > 0) {
    
    # Check if user is currently running any cluster jobs (see myRfunctions.R)
    check_cluster_jobs(o, action = o$cluster_conflict_action)
    
    # Create a new log file for the cluster jobs (see myRfunctions.R)
    log_file = create_bash_log(o$pth$log, log = o$log_file, err = o$err_file)
    
    # Construct sbatch array command for running in parallel
    sbatch_array = paste0("--array=1-", n_jobs, "%", o$job_limit)
    
    # TODO: Perform some checks on these user defined options...
    
    # Format user-defined options into sbatch-interpretable options
    sbatch_options = c(paste0("--partition=", o$cluster_partition), 
                       paste0("--mem=",       o$job_memory), 
                       paste0("--qos=",       o$job_queue))
    
    # If using scicore partition, include run time option (unlimited on covid19 partition)
    if (o$cluster_partition == "scicore")
      sbatch_options = c(sbatch_options, paste0("--time=", o$job_time))
    
    # Collapse into a whitespace-seperated string
    sbatch_options = paste0(sbatch_options, collapse = " ")
    
    # Extract and format arguments to bash file (including log file)
    bash_inputs = paste(paste0(list(...), collapse = " "), log_file)
    
    # Concatenate system command
    sys_command = paste("sbatch", sbatch_options, sbatch_array, bash_file, bash_inputs)
    
    # Invoke this command
    system(sys_command)
    
    # Wait for all cluster jobs to complete (see myRfunctions.R)
    wait_for_jobs(o, log_file, n_jobs)
  }
}

# ---------------------------------------------------------
# Bi-directional setdiff - elements not in both x and y
# ---------------------------------------------------------
symdiff = function(x, y) setdiff(union(x, y), intersect(x, y))

# ---------------------------------------------------------
# Format a number with thousand mark separators
# ---------------------------------------------------------
thou_sep = function(val) {
  format_val = format(val, scientific = FALSE,
                      trim = TRUE, 
                      drop0trailing = TRUE, 
                      big.mark = ",")
  return(format_val)
}

# ---------------------------------------------------------
# Load an file if it exists, throw an error if not
# ---------------------------------------------------------
try_load = function(pth, file, msg = NULL, type = "rds", throw_error = TRUE, sep = FALSE) {
  
  # Initiate trivial output
  file_contents = NULL
  
  # Set default error message
  if (is.null(msg))
    msg = "Cannot load file"
  
  # Switch case for loading function
  loading_fnc = switch(
    tolower(type), 
    
    # Support both RDS and CSV
    "rds" = "readRDS", 
    "csv" = "read.csv",
    
    # Throw an error if anything else requested
    stop("File type '", type, "' not supported")
  )
  
  # Concatenate path and file name
  file_name = paste0(pth, ifelse(sep, file_sep(), ""), file, ".", type)
  
  # If file doesn't exist, throw an error if desired
  if (!file.exists(file_name) && throw_error == TRUE)
    stop(msg, " [missing: ", file_name, "]")
  
  # If file exists, try to load it
  if (file.exists(file_name)) {
    
    # Get the loading function and attempt to load file
    file_contents = tryCatch(
      get(loading_fnc)(file_name),
      
      # Catch the error - we may not want to throw it
      error = function(e) {
        
        # Throw descriptive error if desired
        if (throw_error == TRUE) 
          stop(msg, " [unreadable: ", file_name, "]")
      }
    )
  }
  
  return(file_contents)
}

# ---------------------------------------------------------
# Unlist and return names with seperator of choice
# ---------------------------------------------------------
unlist_format = function(x, sep = "$", ...) {
  y = unlist(x, ...)
  names(y) = gsub("\\.", sep, names(y))
  return(y)
}

# ---------------------------------------------------------
# Wait until all cluster jobs have finished
# ---------------------------------------------------------
wait_for_jobs = function(o, log_file, n_lines, 
                         wait_time = 1, pause_time = 5) {
  
  # Wait for log file to be created
  while (!file.exists(log_file)) Sys.sleep(wait_time)
  
  # Initiate a progress bar
  pb = start_progress_bar(n_lines)
  
  # Wait for at least one line to be written
  while (file.info(log_file)$size == 0) Sys.sleep(wait_time)
  
  # ---- Continuously update progress bar ----
  
  # Wait for all jobs to write to log file
  while (nrow(read.table(log_file)) < n_lines) {
    
    # Update progress bar
    setTxtProgressBar(pb, nrow(read.table(log_file)))
    
    # Check number of running and pending jobs
    n_jobs = n_slurm_jobs(user = o$user)
    
    # Have all jobs already finished?
    if (sum(unlist(n_jobs)) == 0) {
      
      # Pause - let any recently completed jobs write to log
      Sys.sleep(pause_time)
      
      break # Break out of while loop
    }
    
    # Wait before testing again
    Sys.sleep(wait_time)
  }
  
  # Finalise progress bar
  setTxtProgressBar(pb, n_lines)
  close(pb)
  
  # ---- Report any batch job failures ----
  
  # Number of jobs not reported in log file
  n_missing = n_lines - nrow(read.table(log_file))
  
  # A problem if this is non-zero
  if (n_missing > 0) {
    
    # Append this error message to an error log file
    err_message = paste0(" ! ", n_missing, " batch job(s) did not finish")
    write(err_message, file = paste0(o$pth$log, o$err_file), append = TRUE)
  }
}

