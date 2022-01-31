###########################################################
# NETWORKS
#
# Create contact edge list with desired network structure.
#
# The network_structure = "layers" functionality is currently
# under development. This will be available in version 3.0.
#
###########################################################

# ---------------------------------------------------------
# Parent function for create network
# ---------------------------------------------------------
create_network = function(o, p, ppl, do_plot, verbose) {
  
  if (verbose != "none")
    message("  > Creating contact network: ", p$network_structure)
  
  # Do not allow extremely low number of contacts - will likely lead to issues
  if (p$contacts < 1e-6)
    stop("Average number of daily contacts is unfeasibly low")
  
  # ---- Network layers that can be combined ----
  
  # Create layered network
  if (p$network_structure == "layers") {
    
    # This functionality will be available in version 3.0 - check back soon!
    stop("Network 'layers' are currently under development")
    
    # Other standard setup for all layered network types
    p = setup_layers(p)
    
    # Temporarily append an age_group variable in ppl
    ppl[, age_group := cut(age, breaks = c(-Inf, p$network$age_breaks), labels = FALSE)]
    
    # Initiate edge list - build this up with layers
    elist = NULL
    
    # Create household network layer
    if ("household" %in% p$network_layers)
      elist = create_household(o, p, ppl, elist)
    
    # Create school network layer
    if ("school" %in% p$network_layers)
      elist = create_school(o, p, ppl, elist)
    
    # Create workplace network layer
    if ("workplace" %in% p$network_layers)
      elist = create_workplace(o, p, ppl, elist)
    
    # Finally, create residual 'other' layer
    elist = create_residual(o, p, ppl, elist)
    
    # Plot network structure example (only if we have a household layer)
    if (do_plot == TRUE & "household" %in% p$network_layers)
      plot_network_structure(o, p, ppl, elist)
    
    # Now ready to remove temporary 'age groups' variable
    ppl[, age_group := NULL]
  }
  
  # ---- Standalone networks ----
  
  # Create basic small-world, age-structured network
  if (p$network_structure == "age")
    elist = create_age(o, p, ppl)
  
  # Create uber basic random network
  if (p$network_structure == "random")
    elist = create_random(o, p, ppl)
  
  # ---- Perform sanity checks ----
  
  # Throw an error if no edge list created
  if (is.null(elist))
    stop("No edge list created - investigation needed")
  
  # Throw an error if average number of contacts is drastically wrong
  mean_contacts = nrow(elist) / p$population_size
  if (abs(mean_contacts - p$contacts) > 0.1)
    stop("Requested an average of ", p$contacts, " contacts but generated ", round(mean_contacts, 3))
  
  # Calculate age correction factor
  network = age_correction_factor(p, ppl, elist)
  
  return(network)
}

# ---------------------------------------------------------
# Standard setup for all layered network types
# ---------------------------------------------------------
setup_layers = function(p, age_bin = 5) {
  
  # ---- IDs, names, and details or layers ----
  
  # Labels and construction properties for the different layers 
  p$network$layers = list(
    "h" = list(label = "household", 
               age_band      = NA,
               distrib       = "normal", 
               distrib_param = list(mean = 2.5, sd = 1)),
    "s" = list(label = "school", 
               age_band      = c(0, 20), 
               ppl_list      = numeric(0),
               distrib       = "normal", 
               distrib_param = list(mean = 8, sd = 1)),
    "w" = list(label = "workplace", 
               age_band      = c(21, 65), 
               average_size  = 20, 
               ppl_list      = numeric(0),
               distrib       = "normal", 
               distrib_param = list(mean = 6, sd = 1)),
    "c" = list(label = "community")
  )
  
  # ---- Initialize age breaks ----
  
  # Age breaks up to maximum modelled age, based on value of age_bin
  age_breaks = c(seq(from = age_bin, to = max(p$age$all), by = age_bin), Inf)
  
  # Associated labels (rename last label into something prettier)
  labels_breaks = levels(cut(0, breaks = c(0, age_breaks), include.lowest = TRUE))
  labels_breaks[length(age_breaks)] = paste0(age_breaks[length(age_breaks)-1], "+")
  
  # Append age breaks and labels to p list
  p$network$age_breaks = age_breaks
  p$network$labels_breaks = labels_breaks
  
  # ---- Import POLYMOD data ----
  
  # Load POLYMOD contact matrix
  age_bins = c(0, age_breaks[is.finite(age_breaks)])
  p$polymod = load_polymod(p, age_bins = age_bins)
  
  # Align age groups between POLYMOD and OpenCOVID
  list_age_groups = 1 : nrow(p$polymod$matrix)
  names(p$polymod$matrix) = as.character(list_age_groups)
  p$polymod$matrix = cbind(age_group = list_age_groups, p$polymod$matrix)
  
  return(p)
}

# ---------------------------------------------------------
# Create household network layer
# ---------------------------------------------------------
create_household = function(o, p, ppl, elist) {
  
  # Create layer template: household
  household_template = create_layer_template(p$network, ppl, "h")
  
  # Populate households
  households = sample_layer_members(p$network, "h", ppl, 
                                    household_template, 
                                    p$polymod$matrix)
  
  # Check whether all people are allocated to households (and not more)
  if(length(unique(unlist(households))) != nrow(ppl))
    stop("Number of people allocated to households differs from total population")
  
  # Construct contact matrix for households and append network layer ID
  household_matrix = span_contact_matrix(households, l = "h") %>%
    mutate(layer = p$network$layers$h$label)
  
  # Concatenate with all other layers
  elist = rbind(elist, household_matrix)
  
  return(elist)
}

# ---------------------------------------------------------
# Create school layer (incl. teachers)
# ---------------------------------------------------------
create_school = function(o, p, ppl, elist) {
  
  # Extract parameters for school network
  param_schools = p$network$layers[["s"]]
  
  # Remove people not relevant for this layer
  ppl_schools = ppl %>%
    select(id, age, age_group) %>%
    filter(age >= param_schools$age_band[1], 
           age <= param_schools$age_band[2])
  
  # Create layer template: school
  school_template = create_layer_template(p$network, ppl_schools, "s")
  
  # Populate school classrooms
  schools = sample_layer_members(p$network, "s", ppl_schools, 
                                 school_template, p$polymod$matrix)
  
  # Sanity check that we haven't lost anyone
  if(length(unique(unlist(schools))) != nrow(ppl_schools))
    stop("School members not sampled correctly")
  
  # Sample teacher from suitable age bracket
  ppl_teachers = ppl %>% 
    select(id, age, age_group) %>%
    filter(age > param_schools$age_band[2] & age < 65) %>% 
    sample_n(size = nrow(school_template), replace = FALSE)
  
  # Map teachers to classrooms
  schools = Map(c, schools, ppl_teachers$id)
  
  # Construct contact matrix for schools and append network layer ID
  school_matrix = span_contact_matrix(schools, l = "s") %>%
    mutate(layer = p$network$layers$s$label)
  
  # Concatenate with all other layers
  elist = rbind(elist, school_matrix)
  
  return(elist)
}

# ---------------------------------------------------------
# Create Workplace layer
# ---------------------------------------------------------
create_workplace = function(o, p, ppl, elist) {
  
  # Extract parameters for workplace network
  param_workplace = p$network$layers[["w"]]
  
  # Remove people not relevant for this layer
  ppl_workplace = ppl %>%
    select(id, age, age_group) %>%
    filter(age >= param_workplace$age_band[1], 
           age <= param_workplace$age_band[2])
  
  # Create layer template: school
  workplace_template = create_layer_template(p$network, ppl_workplace, "w")
  
  # Populate workplaces
  workplaces = sample_layer_members(p$network, "w", ppl_workplace, 
                                    workplace_template, p$polymod$matrix)
  
  # Sanity check that we haven't lost anyone
  if (length(unique(unlist(workplaces))) != nrow(ppl_workplace))
    stop("Workplace members not sampled correctly")
  
  # Construct contact matrix for workplace and append network layer ID
  workplace_matrix = span_contact_matrix(workplaces, l = "w") %>%
    mutate(layer = p$network$layers$w$label)
  
  # Concatenate with all other layers
  elist = rbind(elist, workplace_matrix)
  
  return(elist)
}

# ---------------------------------------------------------
# Create POLYMOD residual layer
# ---------------------------------------------------------
create_residual = function(o, p, ppl, elist) {
  
  # Calculation population statistics by age group
  pop_age_structure = ppl %>% group_by(age_group) %>% summarize(count = n())
  
  # Amend current edgelist by age groups
  current_elist = convert_edgelist(p, elist, ppl)
  
  # Scale POLYMOD matrix to simulated population (to work with absolute contacts)
  pm_matrix = unname(as.matrix(p$polymod$matrix[, -1]))
  pm_demo   = p$polymod$pop_size
  
  # Ensures age groups are maintained regardless of selected layers
  pm_demo_matrix = matrix(data = rep(pm_demo, times = length(pm_demo)),
                          ncol = length(pm_demo),
                          byrow = FALSE)
  
  # With the given population, how many contacts occur according to POLYMOD?
  pm_contacts_scaled = (pm_matrix * pm_demo_matrix) / sum(pm_demo) * 
    sum(pop_age_structure$count)
  
  # POLYMOD scaler to achieve average number of contacts
  if(!is.null(p$contacts)) {
    
    # Desired number of contacts
    pm_avg_original    = sum(pm_matrix * pm_demo_matrix) / sum(pm_demo)
    pm_contacts_scaled = round(pm_contacts_scaled * p$contacts / pm_avg_original)
  }
  
  # Reshape and plot scaled POLYMOD data (target state)
  pm_contacts_long = reshape2::melt(unname(pm_contacts_scaled), 
                                    varnames = c("from_age_group", "to_age_group"),
                                    value.name = "count")
  
  # ---- Calculate residual (POLYMOD - everything else) ----
  
  # How many contacts are still required?
  clean_current_elist = unname(as.matrix(current_elist$wide))
  clean_current_elist[is.na(clean_current_elist)] = 0
  
  # Age structure of contacts still to be sampled
  residual_matrix = pm_contacts_scaled - clean_current_elist
  
  # Due to symmetry of contact matrix
  residual_matrix[lower.tri(residual_matrix)] = NA
  
  # Ignore cases when household contacts exceed POLYMOD expectations
  residual_matrix[residual_matrix < 0] = 0
  
  # Contacts on the diagonal need to be sampled only once (mirrored automatically)
  diag(residual_matrix) = round(diag(residual_matrix) / 2)
  
  # Sample contacts for the residual layer
  residual_edgelist = sample_residual(residual_matrix, ppl, 
                                      current_elist$ext[, .(from, to)]) %>%
  mutate(layer = p$network$layers$c$label)
  
  # Concatenate with all other layers
  elist = rbind(elist, residual_edgelist)
  
  return(elist)
}

# ---------------------------------------------------------
# Create standalone, basic small-world, age-structured network
# ---------------------------------------------------------
create_age = function(o, p, ppl) {
  
  # Load POLYMOD contact matrix (single ages bins by default)
  polymod_matrix = load_polymod(p)$matrix
  
  # Extract the age groups that were created by socialmixr, usually the oldest individuals are lumped together
  polymod_age_groups = 1 : length(names(polymod_matrix))
  
  # Temporarily append an age_group variable in ppl, to later sample ID's according to age group
  ppl[, age_group := min(which(p$age$all == age[1]), length(polymod_age_groups)), by = list(age)]
  
  # How many individuals per age group have we created? Using data.tables .N variable and order by age_group
  created_demog = ppl[order(age_group), .N, by = list(age_group)]$N
  
  # Multiply the contact matrix by the demography counts and unlist to sample later. Unlist works by column
  contact_probs = unlist(polymod_matrix[, Map("*", .SD, created_demog)], use.names = FALSE)
  contact_probs[is.na(contact_probs)] = 0
  
  # Create ego_cell and alter_cell vectors to get cell identities
  ego_cell   = rep(1 : length(polymod_age_groups), each  = length(polymod_age_groups))
  alter_cell = rep(1 : length(polymod_age_groups), times = length(polymod_age_groups))
  
  # Pre-allocate edgelist
  n_contacts   = round(p$population_size * p$contacts / 2)
  elist = data.table(ego_age_group   = integer(n_contacts), 
                     alter_age_group = integer(n_contacts), 
                     from            = integer(n_contacts), 
                     to              = integer(n_contacts))
  
  # Sample contacts
  realized_contacts = sample(x = 1 : length(contact_probs), 
                             size = n_contacts, 
                             prob = contact_probs, 
                             replace = TRUE)
  
  # Assign the age groups per contacts to the edgelist
  elist[, ego_age_group   := ego_cell[realized_contacts]]
  elist[, alter_age_group := alter_cell[realized_contacts]]
  
  # Now sample from ppl$age_groups, using data.tables by = as well as .N (for outgoing and incoming contacts)
  elist[, from := sample(ppl[age_group == ego_age_group[1],   id], size = .N, replace = TRUE), by = list(ego_age_group)]
  elist[, to   := sample(ppl[age_group == alter_age_group[1], id], size = .N, replace = TRUE), by = list(alter_age_group)]
  
  # Remove self-loops and additional variables
  elist = elist[from != to, c("from", "to")]
  
  # Mirror to reflect two sidedness of contacts
  all_contacts_reverse = data.table(from = elist$to, to = elist$from)
  elist = rbindlist(list(elist, all_contacts_reverse), use.names = FALSE)
  
  # Remove temporary 'age groups' variable
  ppl[, age_group := NULL]
  
  # Append network layer ID
  elist$layer = "community"
  
  return(elist)
}

# ---------------------------------------------------------
# Create standalone, uber basic random network
# ---------------------------------------------------------
create_random = function(o, p, ppl) {
  
  # Number of nodes (people) and edges (contacts)
  n_edges = round(p$population_size * p$contacts / 2)
  
  # Generate simple network
  network = play_erdos_renyi(n = p$population_size, m = n_edges, directed = FALSE)
  
  # Format edges into a datatable - this is what we really need
  network_df = network %>% activate(edges) %>% as.data.table()
  
  # Repeat pair-wise so all contacts are double directed
  reverse_df = data.table(from = network_df$to,
                          to   = network_df$from)
  
  # Bind into single datatable
  elist = rbind(network_df, reverse_df)
  
  # Append network layer ID
  elist$layer = "community"
  
  return(elist)
}

# ---------------------------------------------------------
# Calculate age correction factor
# ---------------------------------------------------------
age_correction_factor = function(p, ppl, elist) {
  
  # Number of contacts per person 
  count_contact = table(elist$from)
  
  # Seperate into IDs and number of contacts
  id = as.numeric(names(unlist(count_contact)))
  n_contacts = as.numeric(unlist(count_contact))
  
  # Breaks to create age bins
  # 
  # NOTE: The 10 year age bins relates the the 10-year prognosis values we use
  age_breaks = seq(0, max(p$age$all) + 1, by = 10)
  
  # Total number of contacts per age group
  age_df = data.table(id = id, n_contacts = n_contacts) %>%
    full_join(ppl[, .(id, age)], by = "id") %>%
    mutate(age_group = cut(age, age_breaks, include.lowest = TRUE)) %>%
    group_by(age_group) %>%
    summarise(n_contacts = sum(n_contacts, na.rm = TRUE))
  
  # Relative probability of a contact (and therefore of infection - when all other factors are equal)
  age_contacts = age_df$n_contacts / sum(age_df$n_contacts)
  
  # Use this to correct prognosis probabilities to achieve required disease and death age-distributions
  age_correction = (1 - age_df$n_contacts / mean(age_df$n_contacts)) + 1  # See fn_prognosis to see this in action
  
  # Append output details to a list
  network = list(edge_list      = elist, 
                 age_contacts   = age_contacts, 
                 age_correction = age_correction)
  
  return(network)
}

# ---------------------------------------------------------
# Load POLYMOD contact matrix (single ages bins by default)
# ---------------------------------------------------------
load_polymod = function(p, age_bins = p$age$all) {
  
  # Create contact matrix from socialmixr
  #
  # NOTE: The contact_matrix function complains when all_ages are not consistent with it's data, however it performs
  #       interpolation between ages to prevent data loss. We therefore justify the use of suppressWarnings here.
  polymod_data = contact_matrix(survey     = polymod, 
                                countries  = p$contact_matrix_countries, 
                                age.limits = age_bins, 
                                symmetric  = TRUE, 
                                quiet      = TRUE) %>% suppressWarnings()
  
  # Convert matrix to datatable and store along with age group pop sizes
  polymod = list(matrix   = data.table(polymod_data$matrix), 
                 pop_size = polymod_data$demography$population)
  
  return(polymod)
}

