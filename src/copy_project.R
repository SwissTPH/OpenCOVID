###########################################################
# COPY PROJECT
#
# Stand-alone script for uploading and downloading scenarios
# results to/from the unit's shared group folder on sciCORE.
# Useful for sharing results with others.
#
# This script is called by bash_copy.sh, and should only 
# be used via that bash script when logged into scicore's
# 'log-in transfer' node.
#
###########################################################

# ---- Default settings ----

# These defaults are overwritten if arguments provided via bash command:
#   sh bash_launch.sh copy_direction analysis_name

# Copy from (download) or to (upload) shared folder
copy_direction = "download" # OPTIONS: "download" or "upload"

# Name of analysis to copy
analysis_name = "demo"

# Directory name in shared group folder
dir_name = "file_transfer"

# ---- Extract inputs (if provided) ----

# Extract input from bash file
args = commandArgs(trailingOnly = TRUE)

# Any additional arguments provided?
if (length(args) > 0) {
  
  # Name arguments when provided from bash
  copy_direction = as.character(args[1])
  analysis_name  = as.character(args[2])
}

# ---- Construct command ----

# Base path to all files belonging to the unit
unit_path = file.path("", "scicore", "home", "penny")

# Get user ID
user = Sys.info()[["user"]]

# Name of your openCOVID repo
repo_name = basename(getwd())

# Construct paths to the user's files and shared group files
user_path  = file.path(unit_path, user, repo_name, "output", "2_scenarios", analysis_name, "")
group_path = file.path(unit_path, "GROUP", "OpenCOVID", dir_name, analysis_name, "")

# Set up rsync command
rsync_call = "rsync -aP"

# Copy access rights & exclude potentially huge simulations sub dir
set_access  = "--chmod=ugo=rwX"
set_exclude = "--exclude simulations"

# Construct upload command (copy TO shared folder)
if (tolower(copy_direction) == "upload")
  d = list(source = user_path, dest = group_path)

# Construct download command (copy FROM shared folder)
if (tolower(copy_direction) == "download")
  d = list(source = group_path, dest = user_path)

# ---- Sanity checks ----

# Check validity of copy direction argument
if (!exists("d"))
  stop("Invalid copy direction '", copy_direction, "' - must be 'upload' or 'download'")

# Check source directory exists
if (!dir.exists(d$source))
  stop("Directory '", d$source, "' does not seem to exist")

# Create destination directory if needed
if (!dir.exists(d$dest))
  dir.create(d$dest, recursive = TRUE)

# Concatenate rysnc command, options, and dirs
rsync_command = paste(rsync_call, set_access, d$source, d$dest, set_exclude)

# Execute rsync command
system(rsync_command)

