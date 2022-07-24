###########################################################
# MY RESULTS
#
# A sandbox for you to produce plots using model outcomes.
# All plotting utility functions can be used here.
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
# - The function format_results can take any scenario file
#   and slice model outcomes according to the arguments
#   provided to figure_properties. This is super useful to
#   collate and slice into plotting dataframes.
# - Check out the save_fig function for easy saving and 
#   consistent-looking figures.
#
###########################################################

my_results = function(o) {
  
  message(" - Custom figures")
  
  
  # Some examples for the 'booster' scenario from 'demo' analysis to get you started.
  # You will need to run steps 1-3 for the 'demo' analysis for these to work.
  
  
  # 1) ---- Plot using in-built plotting functions ----
  
  # Plot the impact of the vaccine scenario against the baseline
  fig_name = c("My awesome plot using in-built plotting functions")
  plot_temporal(o, fig_name, scenarios = "booster")
  
  
  
  
  # 2) ---- A custom plot ----
  
  # Load the scenario file - model output and input contained within
  result = try_load(o$pth$scenarios, "booster")
  
  # Reduce down this output to something plottable
  model_output = result$output %>%
    filter(metric %in% c("all_new_infections", "deaths"), 
           grouping == "priority_group")
  
  # Plot lines over time for each vaccine prioriy group for the two metrics
  g = ggplot(model_output, aes(x = date, y = mean, colour = group)) + 
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = group), 
                linetype = 0, alpha = 0.3) + 
    geom_line(size = 2) +
    facet_wrap(~metric, scales = "free_y")
  
  # Save the plot to file
  fig_save(o, g, "My awesome custom plot")
}

