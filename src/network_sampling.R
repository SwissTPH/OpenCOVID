###########################################################
# NETWORK SAMPLING FUNCTIONS
#
# Layer sampling functions
#
###########################################################

# ---------------------------------------------------------
# Sample template for households
# ---------------------------------------------------------
sample_household_template <- function(p, ppl) {

  # Extract layer properties  
  param <- p$layer_properties$household
  
  if(param$mean > 0) {
    
    # Sample cluster sizes from Poisson distribution, floored at 1
    household_sizes <- list(size = pmax(1, rpois(p$population_size, lambda = param$mean))) %>%
      as.data.table()
  
    # Add column with cumulative cluster size
    household_sizes <- household_sizes[, "cumulative" := cumsum(size)][cumulative <= p$population_size]
    
    # Remove temporary cumulative sum of household sizes
    household_sizes[, cumulative := NULL]
    
    # Check whether all individuals are allocated
    if(sum(household_sizes$size) < p$population_size) {
      
      # If not all are allocated, catch remaining individuals and create one additional cluster
      household_sizes <- rbindlist(list(household_sizes, 
                                      list(size = p$population_size - sum(household_sizes$size))))
      
    }
    
    # Append layer_id and exclude cumulative sum of entity members
    household_sizes <- household_sizes[, "layer_id" := 1:.N][, .(layer_id, layer_size = size)]
    
  } else {
    
    household_sizes <- list()
    
  }
  
  return(household_sizes)
  
}

# ---------------------------------------------------------
# Create edge list for household layer
# ---------------------------------------------------------
sample_household_edgelist <- function(p, ppl, template) {
  
  # Import polymod matrix and exclude age_group column
  polymod_matrix <- p$polymod$matrix[, -1]
  
  # Use very low probabilities when entries are NA (TO-DO: could also be 0?)
  polymod_matrix[is.na(polymod_matrix)] <- 1e-9
  
  # How many individuals were extracted per age group
  created_demog = ppl[order(age_group), .N, by = list(age_group)]$N
  
  # Determine overall contacts per age groups (still matrix format)
  contact_probs = rowSums(polymod_matrix[, Map("*", .SD, created_demog)])
  
  # Summarize contact probabilities in data.table
  contact_probs <- data.table(age_group = 1:length(p$polymod$age_groups),
                              weight = contact_probs)
  
  # Ammend people data frame with weights
  ppl <- contact_probs[ppl, on = .(age_group)]
  
  # Allocate index household member
  template[, id := sample_vec(ppl[, id], 
                              size = nrow(template), 
                              prob = ppl[, weight], replace = FALSE)]
  
  # Join ppl on template (add age_group of index individual)
  template <- ppl[template, .(layer_id, layer_size, id, age_group), on = .(id)]
  
  # Create sampling pool 
  pool <- data.table::copy(ppl[, .(id, age_group)])
  
  # Exclude individuals sampled as index individuals from pool
  pool <- pool[!(id %in% template[, id]), .(id, age_group)]
  
  # Only sample households which have more than one member
  sampling_template <- template[layer_size > 1, ]

  # Preallocate elist
  elist <- NULL
  
  # Only proceed if any households left (circumvent case when household size = 0)
  if(nrow(sampling_template) > 0) {
  
    # Loop through sampling template
    for(i in 1:nrow(sampling_template)) {
      
      idx_age_group <- sampling_template[i, age_group]
      
      idx_polymod <- polymod_matrix[, .SD, .SDcol = as.character(idx_age_group)] %>%
        setnames(as.character(idx_age_group), "weight")
      
      idx_polymod[, age_group := .I]
      
      pool <- pool[idx_polymod, on = .(age_group)]
      
      if(any(is.na(pool[, weight])) || all(pool[, weight] == 0))
        browser()
      
      household_members <- sample_vec(pool[, id], 
                        size = sampling_template[i, layer_size] - 1, 
                        replace = FALSE,
                        prob = pool[, weight])
      
      pool[, weight := NULL]
      
      pool <- pool[!(id %in% household_members),]
      
      household_members <- c(sampling_template[i, id], household_members)
      
      this_elist <- CJ(from = household_members, to = household_members)[from != to,]
      
      elist <- rbindlist(list(elist, this_elist))
      
    }
  
  }
  
  elist$layer <- "household"
  
  return(elist)
  
}

# ---------------------------------------------------------
# Create edge list for school layer
# ---------------------------------------------------------
sample_school_edgelist <- function(p, ppl) {
  
  # ---- Preparations ----
  
  # Shortcut to school parameters
  param <- p$layer_properties$school
  
  # Extract IDs of students; serves as pool for sampling individuals
  students <- ppl[age >= param$age_min & age <= param$age_max, id]
  
  # Preallocate space for school edgelist
  school_elist <- NULL
  
  # ---- Build sampling template ----
  
  # Check whether sampling of class sizes should occur dynamically
  if(param$cluster_sampling$method == "static") {
    
    # How large should each class be?
    class_size <- param$cluster_sampling$value
    
    # Create sampling template for schools
    template <- data.table(class_id = 1:(floor(length(students)/class_size)),
                           class_size = class_size)
    
    if(sum(template$class_size) < length(students))
      template <- rbindlist(list(template, 
                                 list(class_id = nrow(template) + 1,
                                      class_size = length(students) - sum(template$class_size))))
    
  } else {
    
    # Not yet implemented
    browser()
    
  }
  
  # ---- Teacher status ----
  
  # Reset any pre-allocated attributes
  ppl[, teacher := FALSE]
  
  # Randomly sample number of teachers
  ppl[id %in% sample_vec(ppl[age >= p$layer_properties$workplace$age_min & age <= p$layer_properties$workplace$age_max, id], 
                         size = nrow(template), replace = FALSE), teacher := TRUE]
  
  # Extract IDs of teachers
  teachers <- ppl[teacher == TRUE, id]
  
  # ---- Building school edge list ----
  
  # Preallocate output variable
  elists <- NULL
  
  # Loop through template and construct edge list for each class
  for(i in 1:nrow(template)) {
    
    # Random sample of students according to template
    class_members <- sample_vec(students, template[i, class_size], replace = FALSE)
    
    # Ammend teacher ID to class
    class_members <- c(class_members, teachers[i])
    
    # Create age-structured small-world network within class
    elists[[i]] <- sample_subnetwork(p, ppl, class_members, p$layer_properties$school$contacts)
    
    # Remove individuals that were allocated from pool
    students <- setdiff(students, class_members)
    
  }
  
  # Test whether any students left unallocated
  if(length(students) > 0)
    stop("Sampling error in school layer")
  
  # Collapse list of single edge list into single edge list
  school_elist <- rbindlist(elists)
  
  # Append layer tag
  school_elist$layer <- "school"
  
  return(school_elist)
  
}

# ---------------------------------------------------------
# Create edge list for workplace layer
# ---------------------------------------------------------
sample_workplace_edgelist <- function(p, ppl) {
   
  # ---- Preparations ----
  
  # Shortcut to workplace parameters
  param <- p$layer_properties$workplace
  
  # Extract IDs from people who belong to workforce and are NOT teachers
  pool <- ppl[age >= param$age_min & age <= param$age_max & teacher == FALSE, id]
  
  # ---- Create sampling template ----
  
  # Check whether sampling algorithm is static or dynamic
  if(param$cluster_sampling$method == "static") {
    
    workplace_size <- param$cluster_sampling$value
    
    template <- data.table(workplace_id   = 1:(floor(length(pool)/workplace_size)),
                           workplace_size = workplace_size)
    
    if(sum(template$workplace_size) < length(pool))
      template <- rbindlist(list(template, 
                                 list(workplace_id = nrow(template) + 1,
                                      workplace_size = length(pool) - sum(template$workplace_size))))
    
  } else {
    
    browser()
    
  }
  
  # Number of workplace contacts for now - can be made dynamic too
  workplace_contacts <- param$contacts
  
  # ---- Building workplace edge list ----
  
  # Preallocating collection of individuals
  elists <- NULL
  
  # Loop through workplace template and create edge list
  for(i in 1:nrow(template)) {
    
    # Take random sample from workforce pool with size = workplace template
    members <- sample_vec(pool, template[i, workplace_size], replace = FALSE)
    
    # Create age-structured small-world network within workplace
    elists[[i]] <- sample_subnetwork(p, ppl, members, workplace_contacts)
    
    # Remove individuals that were allocated from pool
    pool <- setdiff(pool, members)
    
  }
  
  if(length(pool) > 0)
    stop("Sampling error in workplace layer")
  
  workplace_elist <- rbindlist(elists)
  
  workplace_elist$layer <- "workplace"
  
  return(workplace_elist)
  
}

# ---------------------------------------------------------
# Create standalone, basic small-world, age-structured network
# ---------------------------------------------------------
sample_subnetwork = function(p, ppl, pool_ids, n_contacts) {
  
  pool <- ppl[id %in% pool_ids, ]
  
  # How many individuals per age group have we created? Using data.tables .N variable and order by age_group
  pool_demog <- pool[order(age_group), .N, by = list(age_group)]
  
  polymod_matrix <- p$polymod$matrix[age_group %in% pool_demog$age_group, .SD, 
                                   .SDcol = c("age_group", as.character(pool_demog$age_group))][, -1]
  
  # Multiply the contact matrix by the demography counts and unlist to sample later. Unlist works by column
  contact_probs = unlist(polymod_matrix[, Map("*", .SD, pool_demog$N)], use.names = FALSE)
  contact_probs[is.na(contact_probs)] = 0
  
  # Create ego_cell and alter_cell vectors to get cell identities
  ego_cell   = rep(pool_demog$age_group, each  = length(pool_demog$N))
  alter_cell = rep(pool_demog$age_group, times = length(pool_demog$N))
  
  sample_contacts <- round(length(pool_ids)*n_contacts/2)
  
  # Pre-allocate edgelist
  elist = data.table(ego_age_group   = integer(sample_contacts), 
                     alter_age_group = integer(sample_contacts), 
                     from            = integer(sample_contacts), 
                     to              = integer(sample_contacts))
  
  # Sample contacts
  realized_contacts = sample(x = 1 : length(contact_probs), 
                             size = sample_contacts, 
                             prob = contact_probs, 
                             replace = TRUE)
  
  # Assign the age groups per contacts to the edgelist
  elist[, ego_age_group   := ego_cell[realized_contacts]]
  elist[, alter_age_group := alter_cell[realized_contacts]]
  
  # Now sample from ppl$age_groups, using data.tables by = as well as .N (for outgoing and incoming contacts)
  elist[, from := sample_vec(pool[age_group == ego_age_group[1],   id], size = .N, replace = TRUE), by = list(ego_age_group)]
  elist[, to   := sample_vec(pool[age_group == alter_age_group[1], id], size = .N, replace = TRUE), by = list(alter_age_group)]
  
  # Remove self-loops and additional variables
  elist = elist[from != to, c("from", "to")]
  
  # Mirror to reflect two sidedness of contacts
  all_contacts_reverse = data.table(from = elist$to, to = elist$from)
  elist = rbindlist(list(elist, all_contacts_reverse), use.names = FALSE)
  
  return(elist)
}

# ---------------------------------------------------------
# Auxiliary function: convert edgelists to different formats (wide, long)
# ---------------------------------------------------------
convert_edgelist <- function(p, elist, ppl, na.rm = FALSE) {
  
  # Prepare output format
  output <- list()
  
  if(!is.null(elist)) { 
    
    # Deep copy of input edge list
    elist_ext <- data.table::copy(elist)
    
    # Merge in age group of "from" column
    elist_ext <- elist_ext[ppl, on = list(from = id), from_age_group := i.age_group]
    
    # Merge in age group of "to" column
    elist_ext <- elist_ext[ppl, on = list(to = id), to_age_group := i.age_group]
    
    # Copy into output format
    output$ext <- elist_ext
    
    # Create 
    elist_long <- elist_ext[order(from_age_group, to_age_group), .N, by = .(from_age_group, to_age_group)]
    
    output$long <- elist_long
    
    age_groups <- p$polymod$matrix[["age_group"]]
    
    empty_elist <- CJ(from_age_group = age_groups, to_age_group = age_groups)
    empty_elist[, "N" := 0]
    
    elist_grouped <- rbindlist(list(elist_long, empty_elist))
    
    elist_grouped <- elist_grouped[elist_grouped[, .I[N == max(N)], by=.(from_age_group, to_age_group)]$V1][order(from_age_group, to_age_group)]
    
    wide_elist <- reshape2::dcast(elist_grouped, from_age_group ~ to_age_group, 
                                  value.var = "N")
    
    output$wide <- setDT(wide_elist)
    
    if(na.rm == TRUE) {
      output$wide[is.na(output$wide)] <- 0
    }
    
  } else {
    
    output <- NULL
    
  }
  
  return(output)
  
}