###########################################################
# MY R FUNCTIONS
#
# A series of helpful R functions.
#
# Written by A.J.Shattock
###########################################################

# ---------------------------------------------------------
# Check if user is currently running any cluster jobs
# ---------------------------------------------------------
check_cluster_jobs = function(o, action = "error") {
  
  # Check number of running and pending jobs
  n_jobs = n_slurm_jobs(user = o$user)
  
  # If this is non-zero then action needed
  if (sum(unlist(n_jobs)) > 0) {
    
    # Throw an error
    if (action == "error")
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
colour_scheme = function(map, pal = NULL, n = 1) {
  
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
    colours = get(pal)(n)	
  
  # A load of colour maps from the pals package
  #
  # See: https://www.rdocumentation.org/packages/pals/versions/1.6
  if (map == "pals")
    colours = get(pal)(n)
  
  # Colour Brewer colour schemes
  if (map == "brewer")
    colours = brewer_pal(palette = first_cap(pal))(n)
  
  # Viridis colour schemes
  if (map == "viridis")
    colours = viridis_pal(option = pal)(n)
  
  # Throw an error if colours not yetr defined
  if (is.null(colours))
    stop("Colour map '", map, "' not recognised (supported: base, pals, brewer, viridis)")
  
  return(colours)
}

# ---------------------------------------------------------
# Create a log file (see function wait_for_jobs)
# ---------------------------------------------------------
create_bash_log = function(pth, log = NULL, err = NULL) {
  for (this_name in c(log, err)) {
    this_file = file.path(pth, this_name)
    if (file.exists(this_file)) file.remove(this_file)
    Sys.sleep(0.1)
    file.create(this_file)
    Sys.sleep(0.1)
  }
  return(file.path(pth, log))
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
# Logistic curve
# ---------------------------------------------------------
logistic = function(x, slope, mid, lower = 0, upper = 1) {
  y = upper + (lower - upper) / (1 + (x / mid) ^ slope)
  return(y)
}

# ---------------------------------------------------------
# Parameter transformation: put a probability on the real number line
# ---------------------------------------------------------
logit <- function(p) {
  z <- log(p / (1 - p))
  return(z)
}

# ---------------------------------------------------------
# Parameter transformation: inverse of the above
# ---------------------------------------------------------
logit_inv <- function(p_logit) {
  z <- exp(p_logit) / (exp(p_logit) + 1)
  return(z)
}

# ---------------------------------------------------------
# Number of running and pending jobs on the cluster
# ---------------------------------------------------------
n_slurm_jobs = function(user) {
  
  # Base sq command for user
  sq = paste("squeue -u", user)
  
  # Concatenate full commands
  slurm_running = paste(sq, "-t running | wc -l")
  slurm_pending = paste(sq, "-t pending | wc -l")
  
  # System call to determine number of slurm processes
  n_running = system(slurm_running, intern = TRUE)
  n_pending = system(slurm_pending, intern = TRUE)
  
  # Convert to numeric and minus 1 to get number of jobs
  n_running = as.numeric(n_running) - 1
  n_pending = as.numeric(n_pending) - 1
  
  # Compile into list
  n_jobs = list(running = n_running, pending = n_pending)
  
  return(n_jobs)
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
sample_vec = function(x, ...) x[sample(length(x), ...)]

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
stop_if_errors = function(pth, err, msg = NULL) {
  
  # Set default error message
  if (is.null(msg))
    msg = "Fatal errors when running cluster jobs"
  
  # Check the error file exists
  err_file = file.path(pth, err)
  if (file.exists(err_file)) {
    
    # Load errors from file and convert to vector
    errors = readLines(err_file)
    
    # Stop if any errors, and display them all
    if (length(errors) > 0)
      stop(msg, " (", length(errors), " errors)\n\n", 
           paste(errors, collapse = "\n"))
  }
}

# ---------------------------------------------------------
# Submit jobs to the cluster and wait until they are complete
# ---------------------------------------------------------
submit_cluster_jobs = function(o, n_jobs, bash_file, ...) {
  
  # If user does 
  if (isFALSE(system("sinfo -V") == 0))
    stop("This OpenCOVID feature requires a slurm-based cluster connection") 
  
  # Skip if no jobs to run
  if (n_jobs > 0) {
    
    # Check if user is currently running any cluster jobs (see myRfunctions.R)
    check_cluster_jobs(o, action = o$cluster_conflict_action)
    
    # Create a new log file for the cluster jobs (see myRfunctions.R)
    log_file = create_bash_log(o$pth$code, log = o$log_file, err = o$err_file)
    
    # Construct sbatch array command for running in parallel
    sbatch_array = paste0("--array=1-", n_jobs, "%", o$job_limit)
    
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
  
  # Set default error message
  if (is.null(msg))
    msg = "Cannot find file"
  
  # Switch case for loading function
  loading_fnc = switch(tolower(type), 
                       
                       # Support both RDS and CSV
                       "rds" = "readRDS", 
                       "csv" = "read.csv",
                       
                       # Throw an error if anything else requested
                       stop("File type '", type, "' not supported")
  )
  
  # Concatenate path and file name
  file_name = paste0(pth, ifelse(sep, file_sep(), ""), file, ".", type)
  
  # Check whether file exists
  if (file.exists(file_name)) {
    
    # Get the loading function and load the file
    file_contents = get(loading_fnc)(file_name)
    
  } else {  # File doesn't exist, action needed
    
    # Either throw a descriptive error
    if (throw_error == TRUE) {
      stop(msg, " [missing file: ", file_name, "]")
      
    } else { # ... or provide a trivial value back
      file_contents = NULL
    }
  }
  
  return(file_contents)
}

# ---------------------------------------------------------
# Wait until all cluster jobs have finished
# ---------------------------------------------------------
wait_for_jobs = function(o, log_file, n_lines, 
                         wait_time = 1, pause_time = 2) {
  
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
  
  # Tell user how many jobs did not successfully complete
  if (n_missing > 0)
    message("!! Batch job finished with ", n_missing, " errors !!")
}

