###########################################################
# MY RESULTS
#
# A sandbox for you to produce plots using model outcomes 
# and/or data. All plotting utility functions can be used
# here.
#
# IMPORTANT NOTE: 
# To use, save a new copy of this file as "my_results.R" 
# in the main code directory.
#
# If you find yourself often producing a particular type of 
# figure, you can request this be integrated as a default
# plotting function. Just ask the uber friendly dev team. 
# 
# Some helpful tips:
# - Check out the inputs of the figure_properties function. 
#   This can help you to slice model outcomes in all sorts
#   of different ways.
# - The functions format_data and format_results take any
#   scenario file and slice the data / model results according
#   to the arguments provided to figure_properties. This is
#   super useful to collate and slice into plotting dataframes.
# - Check out the save_fig function for easy saving and 
#   consistent-looking figures.
#
###########################################################

my_results = function(o) {
  
  message("  > Place your custom plotting functions here!")
  
  
  
  # Some examples using the "relax_npi" scenario to get you started...
  
  # NOTES: 
  #  - You will need to actually run the scenario 'relax_npi' for these to actually work
  #  - You could just as easily define a strategy, or even a subset of scenarios from a strategy
  
  
  
  
  
  # 1) ---- Plot using in-built plotting functions ----
  
  # Plot the impact of the vaccine scenario against the baseline
  fig_name = c("My awesome plot using in-built plotting functions")
  plot_temporal(o, "CH", fig_name, scenarios = "relax_npi", everything = TRUE)
  
  
  
  
  
  # 2) ---- A custom temporal plot using the helper functions ----
  
  # Load the scenario file - model output and data and contained within
  a = readRDS(paste0(o$pth$scenario, "CH_relax_npi.rds"))
  
  # Let's call figure_properties to define who we slice the model output and data
  f = figure_properties(o, plot_baseline = FALSE, scenarios = "relax_npi", 
                        plot_metrics = c("confirmed", "deaths"), 
                        plot_by = "age", plot_from = "2020-09-01")
  
  # Apply the slicing to the model results
  model_results = format_results(o, f, a)
  
  # Plot a simple line chart over time for each age group for the two metrics
  g = ggplot(model_results, aes(x = date, y = value, colour = group)) + 
    geom_line() +
    facet_wrap(~metric, scales = "free_y")
  
  # Save the plot to file
  save_fig(o, g, "My awesome plot")
  
  
  
  
  
  # 3) ---- Going it alone! Plotting without the helper functions ----
  
  # Load the scenario file - model output and data and contained within
  a = readRDS(paste0(o$pth$scenario, "CH_relax_npi.rds"))
  
  # Reduce down this output to something plottable
  model_results = a$model_output %>%
    filter(metric %in% c("confirmed", "deaths"), 
           grouping == "vaccine_priority")
  
  # Plot a simple line chart over time for each vaccine priority group for the two metrics
  g = ggplot(model_results, aes(x = date, y = value)) + 
    geom_area(aes(colour = group, fill = group)) +
    facet_wrap(~metric, scales = "free_y")
}

