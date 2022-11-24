
# ---- Data ----

# Data on deaths per million from www.bandolier.org.uk/booth/Risk/dyingage.html
deaths_per_million = c(172, 4881, 9375, 3020, 1852, 885, 350, 145, 54, 18, 7)

# Ages associated with each data point
years = c(0.5, 2.5, seq(10, 90, by = 10))

# Format into datatable and calculate daily rates
data_df = data.table(years  = years, 
                     deaths = deaths_per_million) %>%
  mutate(p_death = (1 / deaths) / 365)

# ---- Model ----

# Years to evaluate 
all_years = seq(0, 89, by = 1)

# Fit to discrete data (crudely)
model_df = logistic(x = all_years, 
                    slope = 10, 
                    mid   = 90, 
                    lower = 0, 
                    upper = 0.0008) %>%
  as_named_dt("p_death") %>%
  mutate(years = all_years)

# ---- Plot ----

# Plot model against the data
g = ggplot(data_df, aes(x = years, y = p_death)) + 
  geom_point(size = 2, colour = "red") + 
  geom_line(data = model_df, colour = "black")

# View the plot
print(g)

