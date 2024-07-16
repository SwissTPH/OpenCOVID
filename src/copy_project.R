###########################################################
# COPY PROJECT
#
# Stand-alone script for uploading and downloading scenarios
# results to/from the unit's shared group folder on sciCORE.
#
# Useful for sharing results with others.
#
###########################################################

# ---- User settings ----

# Copy from (download) or to (upload) shared folder
copy_direction = "upload" # OPTIONS: "download" or "upload"

# Name of analysis to copy
analysis_name = "omicron"

# Name of your openCOVID repo
repo_name = "opencovid"

# ---- Construct command ----

# Base path to all files belonging to the unit
unit_path = file.path("", "scicore", "home", "penny")

# Get user ID
user = Sys.info()[["user"]]

# Construct paths to the user's files and shared group files
user_path  = file.path(unit_path, user, repo_name, "output", "2_scenarios", analysis_name, "")
group_path = file.path(unit_path, "GROUP", "OpenCOVID", "file_transfer", analysis_name, "")

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

