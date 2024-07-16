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
    # stop("Network 'layers' are currently under development")
    
    # Other standard setup for all layered network types
    p = setup_layers(p, age_bin = 1)
    
    # Temporarily append an age_group variable in ppl
    ppl[, age_group := min(which(p$age$all == age[1]), length(p$polymod$age_groups)), by = list(age)]
    
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

    # Finally, create residual 'other' layer if community_sampling is static
    # if(p$community_sampling == "static") {
    elist = create_residual(o, p, ppl, elist)
    # }
    
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
  if (abs(mean_contacts - p$contacts) > 0.1) {
    stop("Requested an average of ", p$contacts, " contacts but generated ", round(mean_contacts, 3))
  }
  
  # Calculate age correction factor
  network = age_correction_factor(p, ppl, elist)
  
  network$enhanced_p <- p
  
  return(network)
}

# ---------------------------------------------------------
# Standard setup for all layered network types
# ---------------------------------------------------------
setup_layers = function(p, ppl, age_bin = 5) {
  
  # ---- Initialize age breaks ----
  
  # Age breaks up to maximum modelled age, based on value of age_bin
  age_breaks = c(seq(from = age_bin, to = max(p$age$all), by = age_bin), Inf)
  
  # Associated labels (rename last label into something prettier)
  labels_breaks = levels(cut(0, breaks = c(0, age_breaks), include.lowest = TRUE))
  labels_breaks[length(age_breaks)] = paste0(age_breaks[length(age_breaks)-1], "+")
  
  # Append age breaks and labels to p list
  p$network$age_breaks    = age_breaks
  p$network$labels_breaks = labels_breaks
  p$network$age_bin       = age_bin
  
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
  
  message("    * Sampling household layer")
  
  # Create household template (layer id, layer size)
  household_template <- sample_household_template(p, ppl)
  
  if(!is.null(nrow(household_template))) {
    
    # Sample edge list for household layer
    elist_household <- sample_household_edgelist(p, ppl, household_template)
    
    # Concatenate household e-list with all other existing e-lists
    elist <- rbindlist(list(elist, elist_household))
    
  }

  return(elist)
}

# ---------------------------------------------------------
# Create school layer (incl. teachers)
# ---------------------------------------------------------
create_school = function(o, p, ppl, elist) {
  
  message("    * Sampling school layer")
  
  # Construct edge list for schools
  elist_school <- sample_school_edgelist(p, ppl)
  
  # Concatenate with all other layers
  elist <- rbindlist(list(elist, elist_school))
  
  return(elist)
  
}

# ---------------------------------------------------------
# Create Workplace layer (excl. teachers)
# ---------------------------------------------------------
create_workplace = function(o, p, ppl, elist) {
  
  message("    * Sampling workplace layer")
  
  # Construct edge list for workplaces
  elist_workplace <- sample_workplace_edgelist(p, ppl)
  
  # Concatenate with all other layers
  elist = rbindlist(list(elist, elist_workplace))
  
  return(elist)
}

# ---------------------------------------------------------
# Create POLYMOD residual layer
# ---------------------------------------------------------
create_residual = function(o, p, ppl, elist) {
 
  message("    * Sampling residual to POLYMOD")
  
  # Summarize existing edge list by age buckets
  current_elist = convert_edgelist(p, elist, ppl)

  # Load POLYMOD contact matrix
  polymod_matrix = p$polymod$matrix[, -1]
  
  polymod_matrix[is.na(polymod_matrix)] <- 0
  
  # Extract age groups from polymod matrix
  polymod_age_groups = 1:length(names(polymod_matrix))

  # Extract simulated people per age group
  created_demog = ppl[order(age_group), .N, by = list(age_group)]$N
  
  # Multiply the contact matrix by the demography counts to get absolute contact numbers
  contact_probs = polymod_matrix[, Map("*", .SD, created_demog)]
  
  # Resolve cases in which contact probs are NA's 
  contact_probs[is.na(contact_probs)] = 0
  
  # Determine scaler to align POLYMOD matrix with average number of calibrated contacts
  polymod_scaler <- p$contacts / (sum(contact_probs)/p$population_size)
  
  # Scale up/down POLYMOD matrix to achieve correct number of average contacts
  scaled_polymod <- (polymod_scaler * contact_probs) %>% as.matrix() %>% unname()
  
  # In case there is no current edge list, use scaled polymod as base for simulation
  if(is.null(current_elist)) {
    
    browser()
    
    residual_matrix <- scaled_polymod
    
    new_contacts <- p$contacts
    
  } else {
  
    # Extract contact matrix already simulated by other layers
    already_simulated <- current_elist$wide[, -1] %>% as.matrix() %>% unname()
    
    # Calculate raw version of community matrix (can contain negative values) 
    residual_matrix_raw <- scaled_polymod - already_simulated
    
    # Floor community matrix at 0
    residual_matrix_bounded <- pmax(residual_matrix_raw, 0)
    
    # Scale community matrix to account for already simulated contacts (that are now listed as negative)
    residual_matrix <- residual_matrix_bounded * 
      (sum(residual_matrix_bounded) / sum(residual_matrix_raw))
    
    new_contacts <- p$contacts - nrow(elist)/p$population_size
    
  }
  
  other_elist <- create_age_subnetwork(o, p, ppl, 
                                       mixing_matrix = residual_matrix,
                                       new_contacts = new_contacts)
  
  # Remove redundant community contacts
  if(!is.null(elist)) {
    
    other_elist <- fsetdiff(other_elist, elist[, .(from, to)])  
    
  }
  
  other_elist$layer <- "community"
  
  full_elist <- rbindlist(list(elist, other_elist))
  
  return(full_elist)
}

# ---------------------------------------------------------
# Create standalone, basic small-world, age-structured network
# ---------------------------------------------------------
create_age_subnetwork = function(o, p, ppl, mixing_matrix = NULL, new_contacts = NULL) {

  # Check whether age mixing pattern is provided, if not, use polymod
  if(is.null(mixing_matrix)) {
    
    # Load polymod matrix
    mixing_matrix <- load_polymod(p, age_bins = p$age$all)$matrix
    
    # Summarize age groups
    age_groups <- 1:length(names(mixing_matrix))  
    
    # How many individuals per age group have we created? Using data.tables .N variable and order by age_group
    created_demog = ppl[order(age_group), .N, by = list(age_group)]$N
    
    # Multiply the contact matrix by the demography counts and unlist to sample later. Unlist works by column
    contact_probs = unlist(mixing_matrix[, Map("*", .SD, created_demog)], use.names = FALSE)
    
  } else {
    
    age_groups <- 1:dim(mixing_matrix)[1]
    
    mixing_matrix <- mixing_matrix %>% as.data.table()
    
    contact_probs = unlist(mixing_matrix, use.names = FALSE)
    
  }
  
  contact_probs[is.na(contact_probs)] = 0
  
  # Create ego_cell and alter_cell vectors to get cell identities
  ego_cell   = rep(1 : length(age_groups), each  = length(age_groups))
  alter_cell = rep(1 : length(age_groups), times = length(age_groups))

  n_contacts <- round(p$population_size * new_contacts / 2)
  
  elist <- data.table(ego_age_group   = integer(n_contacts), 
                      alter_age_group = integer(n_contacts), 
                      from            = integer(n_contacts), 
                      to              = integer(n_contacts))

  # Sample contacts
  realized_contacts = sample(x       = 1:length(contact_probs), 
                             size    = n_contacts, 
                             prob    = contact_probs, 
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
                 pop_size = polymod_data$demography$population,
                 age_groups = 1:ncol(polymod_data$matrix))
  
  return(polymod)
}

