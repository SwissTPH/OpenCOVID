
source("auxiliary.R")

# Days to plot over
days = 0 : 600

# Data to fit to 
data_df = data.frame(days  = c(0, 56, 520),
                     value = c(0.98, 0.49, 0.1))

# Functional form to represent decay
use_fn = "exp_biphasic"  # OPTIONS: "exp_double" or "exp_biphasic"

# Starting points for optimisation
x0 = list(
  exp_double = c(
    b = 0.1, 
    g = 1, 
    h = 0.9, 
    d = 1.5), 
  exp_biphasic = c(
    peak  = 100, 
    p     = 0.83, 
    d1    = 46, 
    d2    = 572, 
    vmax  = 0.91,
    alpha = 0.92, 
    beta  = 24.5))

# Parameter lower bounds
lb = list(
  exp_double   = rep(1e-6, 4), 
  exp_biphasic = rep(1e-6, 7))

# Parameter upper bounds
ub = list(
  exp_double   = rep(Inf, 4), 
  exp_biphasic = c(Inf, 1, Inf, Inf, Inf, Inf, Inf))

# ---------------------------------------------------------
# Objective function for optimisation
# ---------------------------------------------------------
objective_fn = function(x, args_obj = NULL) {
  
  # Combine parmaters we're optimising with evaluation points into single list
  args_fn = c(list(days), as.list(x))
  
  # Call the function we're optimising for with all necessary parameter arguments
  immunity_curve = do.call(use_fn, args_fn)
  
  # Immunity at points we have data for
  immunity_pts = immunity_curve[pmax(data_df$days, 1)]
  
  # Squared distance between the curve and the data
  objective_val = sum((immunity_pts - data_df$value) ^ 2)
  
  if (is.na(objective_val))
    stop("Objective function error")
  
  # Return a list (as required by ASD)
  return(list(y = objective_val))
}

# Set a seed that gets a decent result
set.seed(200)

# Optimise parameters using ASD algorithm
x = asd(objective_fn,
        x0 = x0[[use_fn]], 
        lb = lb[[use_fn]],
        ub = ub[[use_fn]],
        max_iters  = 2000, 
        plot_iters = 500, 
        verbose    = TRUE)$x

# Combine optimised parmaters with evaluation points into single list
args_fn = c(list(days), as.list(x))

# Call the function we're optimising for with all necessary parameter arguments
immunity_fit = do.call(use_fn, args_fn)

# Construct a dataframe ready for plotting
immunity_df = data.frame(days = days, value = immunity_fit)

# Plot the optimised curve with the data
p = ggplot(immunity_df) + 
  aes(x = days, y = value) + 
  geom_line() + 
  geom_point(data = data_df) + 
  ylim(0, 1)

# Show the plot
print(p)

