###########################################################
# VALIDATION
#
# Plot validation of previous analyses. That is, overlay the
# latest data points on top of previous projections.
#
###########################################################

# ---------------------------------------------------------
# Plot a previous analysis with recent data on top
# ---------------------------------------------------------
run_validation = function(o) {
  
  # Only continue if specified by do_step (or forced)
  if (!is.element(5, o$do_step)) return()
  
  # Begin validation one day after originally-available data ends
  o$begin_validation = max(o$dates_data) + 1
  
  # Three types of validation figures
  validation1(o)  # Scenarios overlaid with latest data
  # validation2(o)  # Retrospective scenario
  # validation3(o)  # Age structure: model vs data
}

# ---------------------------------------------------------
# Scenarios overlaid with latest data
# ---------------------------------------------------------
validation1 = function(o) {
  
  message("* Producing validation set A")
  
  # Strategies to plot on temporal figure
  compare_strategies = "s1b" # "s1x"
  
  # Plot data up to this point
  data_to = "2021-08-01"
  
  # Metrics to plot
  metrics_lines = c("osi_level", "n_doses", "confirmed", "hospital_beds", "icu_beds", "deaths")
  metrics_bars  = c("osi_level", "n_doses", "confirmed", 
                    "hospital_admissions", "icu_admissions", "deaths")
  
  # Colour of validation data
  latest_colour = "navyblue" # "deeppink1"
  
  # ---- Load latest data and previous results ----
  
  # Re-create figures 2 and 3 from paper
  g2 = fig2(o, quiet = TRUE, strategies = compare_strategies, 
            override_metrics = metrics_lines, plot_legend = FALSE)
  g3 = fig3(o, quiet = TRUE, override_metrics = metrics_bars)
  
  # Dates to extract the latest data for
  o$dates_data = seq(format_date("2020-01-01"), format_date(data_to), by = "day")
  
  # First load the latest data (before we inherit options from analysis file)
  latest_data_list = load_data(o)[[o$cantons]]
  latest_epi = latest_data_list$epi %>%
    select(-grouping, -group)
  
  # ---- Compile validation data ----
  
  # Dates across which we wish the validate
  all_dates = sort(unique(latest_epi$date))
  validation_dates = all_dates[all_dates >= o$begin_validation]
  
  # Format OSI values (not in data$epi dataframe)
  latest_osi_level = latest_data_list$response %>% 
    filter(date >= o$begin_validation) %>%
    mutate(osi_level = osi_level * 100) %>% 
    fill(osi_level)
  
  # Format number of doses (not in data$epi dataframe)
  latest_n_doses = latest_data_list$vaccine %>%
    filter(date >= o$begin_validation) %>%
    select(date, n_doses = number_doses) %>%
    filter(n_doses > 0)
  
  # Fit a simple moving average model to number of doses
  ma_model = smooth::sma(latest_n_doses$n_doses, order = 7, h = 4, interval = "none")
  ma_vec   = c(ma_model$fitted[-(1 : 4)], ma_model$forecast)
  
  # Metric dictionary
  metric_dict  = o$metric_dict[metrics_lines]
  metric_names = unname(metric_dict)
  
  # Cumulative metric dictionary
  cumulative_dict   = o$cumulative_dict[metrics_bars]
  cumulative_names  = unname(cumulative_dict)
  
  # Combine it all into one dataframe, considering only validation dates
  validation_data = data.table(date = all_dates) %>%
    left_join(latest_osi_level, by = "date") %>%
    left_join(latest_n_doses, by = "date") %>%
    pivot_longer(cols = -date, 
                 names_to = "metric") %>%
    rbind(latest_epi) %>%
    filter(!is.na(value)) %>%
    mutate(when = ifelse(date >= o$begin_validation, 
                         "retrospective", "past")) %>%
    mutate(metric_id = metric, .before = "metric", 
           metric = recode(metric_id, !!!metric_dict), 
           metric = factor(metric, levels = metric_names))
  
  # Plot all validation data as points...
  validation_points = validation_data %>%
    filter(metric_id %in% metrics_lines) %>%
    rbind(mutate(validation_data[1, ], value = NA, when = "model")) %>%
    mutate(when = factor(when, levels = c("past", "retrospective", "model")))
  
  # ... Aside from number of doses, for which we'll plot a moving average
  validation_line = data.table(date = validation_dates, 
                               metric_id = "n_doses", 
                               value = ma_vec) %>%
    mutate(metric = recode(metric_id, !!!metric_dict), .after = "metric_id", 
           metric = factor(metric, levels = metric_names))
  
  # Sum cumulative metrics for comparing to impact bar graph
  validation_bars = validation_data %>%
    filter(date >= o$begin_validation) %>% 
    filter(metric_id %in% metrics_bars) %>%
    group_by(metric_id, metric) %>%
    summarise(value = sum(value)) %>%
    mutate(value = ifelse(metric_id == "osi_level", value / length(validation_dates), value), 
           metric = recode(metric_id, !!!cumulative_dict), 
           metric = factor(metric, levels = cumulative_names)) %>%
    as.data.table()
  
  # ---- Place data into manuscript figures ----
  
  # Scenario colours
  parent_strategies = paste0(str_sub(compare_strategies, 1, 2), "x")
  scenario_colours  = fig_colours(parent_strategies)

  # Plot the latest validation data
  g2 = g2 + new_scale("colour") + new_scale("fill") +
    geom_point(data = filter(validation_points, metric_id != "n_doses"), 
               mapping = aes(colour = when)) +
    geom_point(data = filter(validation_points, metric_id == "n_doses"), 
               mapping = aes(colour = when), alpha = 0.3) +
    scale_colour_manual(values = c(o$data_colour, latest_colour, scenario_colours), 
                        labels = c("Data available at time of analysis", 
                                   "Data available retrospectively", 
                                   "Model projections (scenario: relax 22 March, fast vaccine rollout)")) + 
    geom_line(data = validation_line, colour = latest_colour, size = 2) + 
    guides(color = guide_legend(override.aes = list(size = 5), ncol = TRUE)) +
    theme(strip.text  = element_text(size = 16), 
          legend.text = element_text(size = 18))
  
  # Save figure
  save_fig(o, g2, "Validation - fig 2")
  
  # Plot the latest validation data
  g3 = g3 + geom_hline(data = validation_bars, aes(yintercept = value),
                       colour = latest_colour, linetype = "dotdash", size = 2) +
    theme(strip.text = element_text(size = 15))

  # Save figure
  save_fig(o, g3, "Validation - fig 3")
}

# ---------------------------------------------------------
# Retrospective scenario
# ---------------------------------------------------------
validation2 = function(o) {
  
  message("* Producing validation set B")
  
  # Metrics to plot
  metrics = c("osi_level", "n_doses", "confirmed", "hospital_beds", "icu_beds", "deaths")
  
  # Metric colours
  colours = colour_scheme("brewer::dark2", n = length(metrics))
  
  # Validation data colour
  data_colour = "#777777"
  
  # ---- Load latest data ----
  
  # When to begin validation, and how long it runs for
  o$dates_data = seq(format_date("2020-01-01"), format_date("2021-09-01"), by = "day")
  
  # Only one scenario here: retrospective scenario using actual NPIs and vaccine roll out
  f = figure_properties(o, alt_baseline = "retrospective", plot_metrics = metrics)
  
  # First load the latest data (before we inherit options from analysis file)
  latest_data_list = load_data(o)$CH
  latest_df = latest_data_list$epi %>%
    filter(metric %in% metrics, 
           date >= o$begin_validation) %>%
    mutate(metric = factor(f$metric_names[metric], 
                           levels = f$metric_names)) %>%
    select(-grouping, -group)
  
  # ---- Load model outputs ----
  
  # Load retrospective scenario outcomes
  a = try_load(o$pth$scenario, paste0("CH_", f$scenarios))
  
  # Extract and format model output
  model_df = format_results(o, f, a) %>%
    select(date, metric, value) %>%
    filter(date %in% o$dates_data) %>%
    mutate(metric = factor(f$metric_names[metric], 
                           levels = f$metric_names))
  
  # Extract and format data used in scenario
  data_df = format_data(o, f, a) %>%
    select(date, metric, value) %>%
    filter(date %in% o$dates_data) %>%
    mutate(metric = factor(f$metric_names[metric], 
                           levels = f$metric_names))
  
  # ---- Produce plot ----
  
  # Plot model with available and retrospective data on top
  g = ggplot(model_df, aes(x = date, y = value, colour = metric)) + 
    geom_line(size = 2, show.legend = FALSE) + 
    geom_point(data = data_df, colour = data_colour) + 
    geom_point(data = latest_df, colour = "black") + 
    facet_wrap(~metric, scales = "free_y")
  
  # Add a vertical reference line where data end and prediction starts
  g = g + geom_vline(xintercept = o$begin_validation - 1, 
                     linetype = "dashed", colour = o$dash_colour)
  
  # Set colours
  g = g + scale_colour_manual(values = colours) + 
    scale_fill_manual(values = colours)
  
  # Set x axes to user-defined limits (see options.R)
  g = g + scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.05))) + 
    scale_x_date(date_breaks = o$x_tick_dates, date_labels = o$x_tick_format)
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(plot.title   = element_text(size = f$fontsize[1], hjust = 0.5),
          strip.text   = element_text(size = f$fontsize[2]),
          axis.text.y  = element_text(size = f$fontsize[3]),
          axis.text.x  = element_text(size = f$fontsize[3], hjust = 1, angle = 50),
          axis.title   = element_blank())
  
  # Save figure
  save_fig(o, g, "Validation - model")
}

# ---------------------------------------------------------
# Age structure: model vs data
# ---------------------------------------------------------
validation3 = function(o) {
  
  message("* Producing validation set C")
  
  # Age-disaggregated metrics to plot
  epi_metrics = c("confirmed", "deaths") # hospital_admissions
  all_metrics = c(epi_metrics, "total_vaccinated")
  
  # Name for vaccine metric
  vax_name = "Number vaccinated (at least one dose)"
  
  # Lower bound of all age groups
  o$plot_ages = seq(0, 80, by = 10)
  
  # Corresponding strings to represent these groups
  o$age_groups = c(paste0(seq(0, 70, by = 10), "-", seq(9, 79, by = 10)), "80+")
  
  # Scenario fo model outcomes
  scenario = "s1b.fastvax_npi1" # s1b.fastvax_npi1 s0b.fastvax_npi0 retrospective
  
  # Plot up to this date
  plot_to = "2021-09-01" # "2021-03-05" "2021-09-01"
  
  # Wrap strip names longer than n characters
  n_wrap = 22
  
  # ---- Model outputs by age ----
  
  # Load baseline scenario model output
  a = try_load(o$pth$scenario, paste0("CH_", scenario))
  
  # Figure properties - just this baseline, but by age
  f = figure_properties(o, alt_baseline = scenario, plot_metrics = all_metrics, 
                        plot_by = "age", plot_to = plot_to, override_colours = NA)
  
  # Over number vaccinated label
  f$metric_names[length(epi_metrics) + 1] = vax_name
  
  # Format model ouput
  model_df = format_results(o, f, a) %>%
    select(date, age = group, metric, value) %>%
    mutate(age = str_remove(age, "^Age "),
           metric = f$metric_names[metric],
           type = "Model")
  
  # ---- Data by age ----
  
  # Dates for data and model extraction
  o$dates_data = seq(format_date("2020-01-01"), format_date(plot_to), by = "day")
  
  # Load the latest data
  epi_df = load_data_epi(o, by_age = TRUE) %>% 
    filter(canton == "CH")
  
  # This data is weekly - we'll convert to daily by linear interpolation
  epi_weekly = epi_daily = select(epi_df, date, age, all_of(epi_metrics))
  
  # Loop through metrics
  for (metric in epi_metrics) {
    
    # Convert to wide format to interpolate one age group and one metric at a time
    epi_wide = epi_weekly %>% 
      select(date, age, all_of(metric)) %>%
      pivot_wider(id_cols = date, 
                  names_from = age,
                  values_from = metric)
    
    # For each age group, interpolate between weekly estimates
    for (age in o$age_groups)
      epi_wide[[age]] = na.approx(epi_wide[[age]], maxgap = 7, rule = 2) / 7
    
    # Re-convert back to long format
    epi_long = pivot_longer(epi_wide, cols = -date, names_to = "age")
    
    # Insert these values into daily dataframe
    epi_daily = epi_daily %>% 
      left_join(epi_long, by = c("date", "age")) %>%
      select(-all_of(metric)) %>%
      rename(!!metric := value)
  }
  
  # Convert all metrics back to long format
  epi_data_df = epi_daily %>%
    pivot_longer(cols = -c(date, age), 
                 names_to = "metric") %>%
    filter(!is.na(value)) %>%
    mutate(type = "Data") %>%
    arrange(metric, date, age) %>% 
    as.data.table()
  
  # A similar process for number vaccinated...
  
  # Load the latest data
  vax_df = load_data_vaccine(o, by_age = TRUE)
  
  # This data is weekly - we'll convert to daily by linear interpolation
  vax_weekly = select(vax_df, date, age, number_vaccines)
  
  # Convert to wide format to interpolate one age group at a time
  vax_wide = vax_weekly %>% 
    pivot_wider(id_cols = date, 
                names_from = age, 
                values_from = number_vaccines)
  
  # For each age group, interpolate between weekly estimates
  for (age in o$age_groups) {
    interp = na.approx(vax_wide[[age]], maxgap = 7, rule = 2)
    vax_wide[[age]] = round(cum_sum(interp / 7, na.rm = TRUE))
  }
  
  # Convert back to long format
  vax_data_df = vax_wide %>% 
    pivot_longer(cols = -date, 
                 names_to = "age") %>%
    filter(!is.na(value)) %>%
    mutate(metric = "total_vaccinated", .before = "value") %>%
    mutate(type = "Data") %>%
    arrange(age, date) %>% 
    as.data.table()
  
  # Combine epi and vaccine data
  data_df = rbind(epi_data_df, vax_data_df) %>%
    mutate(metric = f$metric_names[metric])
  
  # ---- Produce age plot for all metrics ----
  
  # Concatenate model and data dataframes
  plot_df = rbind(data_df, model_df) %>%
    mutate(type = as.factor(type), 
           metric = factor(metric, levels = f$metric_names))
  
  # Plot areas
  g = ggplot(plot_df, aes(x = date, y = value, color = age, fill = age)) + 
    geom_area() + 
    facet_grid(metric~type, scales = "free_y", 
               labeller = as_labeller(label_wrap_gen(n_wrap)))
  
  # Apply a title
  g = g + ggtitle("Epidemiological outcomes by age: model vs data")
  
  # Set x axes to user-defined limits (see options.R)
  g = g + scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.05))) + 
    scale_x_date(date_breaks = o$x_tick_dates, date_labels = o$x_tick_format)
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(plot.title   = element_text(size = 32, hjust = 0.5),
          strip.text   = element_text(size = 16),
          axis.text.y  = element_text(size = 12),
          axis.text.x  = element_text(size = 12, hjust = 1, angle = 50),
          axis.title   = element_blank(),
          legend.text  = element_text(size = 14),
          legend.title = element_blank(),
          legend.key   = element_blank(),
          legend.key.height = unit(2, "lines"),
          legend.box.background = element_rect(), 
          axis.line    = element_blank(), 
          panel.border = element_rect(colour = "black", fill = NA, size = 1))
  
  # Save figure
  save_fig(o, g, "Validation - age")
  
  # ---- Produce age plot for all metrics ----
  
  pop_data = load_data_pop(o) %>% 
    mutate(age = o$age_groups[age %/% 10 + 1]) %>% 
    group_by(age) %>% 
    summarise(pop = sum(value)) %>% 
    as.data.table()
  
  coverage_df = plot_df %>%
    filter(metric == vax_name, 
           value >= 1) %>%
    left_join(pop_data, by = "age") %>%
    mutate(value = 100 * value / pop, 
           age = factor(age, levels = rev(o$age_groups))) %>%
    select(-metric, -pop)
  
  g = ggplot(coverage_df, aes(x = date, y = value, linetype = type)) + 
    geom_line(size = 1.5) + 
    facet_wrap(~age)
  
  # Apply a title
  g = g + ggtitle("Vaccination coverage (%) by age: model vs data")
  
  # Set x axes to user-defined limits (see options.R)
  g = g + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + 
    scale_x_date(date_breaks = o$x_tick_dates, date_labels = o$x_tick_format)
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(plot.title   = element_text(size = 32, hjust = 0.5),
          strip.text   = element_text(size = 16),
          axis.text.y  = element_text(size = 12),
          axis.text.x  = element_text(size = 12, hjust = 1, angle = 50),
          axis.title   = element_blank(),
          legend.text  = element_text(size = 14),
          legend.title = element_blank(),
          legend.key   = element_blank(),
          legend.key.height = unit(2, "lines"),
          legend.key.width  = unit(2, "lines"),
          legend.box.background = element_rect(), 
          axis.line    = element_blank(), 
          panel.border = element_rect(colour = "black", fill = NA, size = 1))
  
  # Save figure
  save_fig(o, g, "Validation - vaccine")
}

