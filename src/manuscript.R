###########################################################
# MYRESULTS file to create most thesis-specific plots 
###########################################################

# ---------------------------------------------------------
# Main function to create custom plots
# ---------------------------------------------------------
plot_publication_figures = function(o) {
  
   message(" - Produce figures for manuscript")
   
   # Which plots should be produced?
   do_figure_2_s2_s3 <- TRUE
   do_figure_3 <- TRUE
   do_figure_4 <- TRUE
   do_figure_s5 <- TRUE
   
   # ---- Prepare data ----
   
   # Specify metrics that should be loaded from raw files
   metrics <- c("all_new_infections", "hospital_admissions")
   
   # Import or extract data for plotting   
   results <- load_manuscript_data(o, metrics)
   
   # Load dictionary stored in external csv file
   dict <- load_dictionary(o)
      
   # Add factor variables to some elements for easier filtering
   for(f in c("output", "cum_output", "stats", "output_raw", "stats_raw")) {
      results[[f]] <- add_factors(results[[f]], dict)
   }
   
   # Merge clustering data into statistics data
   results$stats <- results$network[results$stats, on = .(scenario)]
   results$stats_raw <- results$network[results$stats_raw, on = .(scenario)]
   
   # ---- Is this the scenario to determine the alternative beta? -----
   if(grepl("beta", o$analysis_name)) {
      
      # Produce diagnostic and time series plots
      assess_betas(o, results, dict)
      
      # End program here; purpose of this analysis is to determine adjusted betas
      return(NULL)
  
   }
      
   # ---- Plotting commands ----
   
   # ---- Figure 2: Temporal plot ----
   if(do_figure_2_s2_s3) {
      
      message("   * Plot Figures 2, S2, and S3")
      plot_figure_2_s2_s3(o, results, plot_metric = "all_new_infections", dict)
      
   }
      
   # ---- Figure 3: Cluster plots ----
   if(do_figure_3) {
      
      message("   * Plot Figure 3")
      plot_figure_3(o, results, plot_metric = "all_new_infections", dict)
      
   }
   
   # ---- Figure 4: Intervention plots ----
   if(do_figure_4) {
      
      message("   * Plot Figure 4")
      plot_figure_4(o, results, m, dict, sim_horizon)
      
   }   
   
   if(do_figure_s5) {
      
      message(paste0("   * Plot Figure S5"))
      
      plot_df = list()
      
      results_df = results$output_raw[(metric == "all_new_infections" & network_structure == "single_layer") | 
                                         (metric == "all_new_infections" & network_structure == "multi_layer" &
                                             grepl("altbeta",scenario) & contact_weight == "UNI" & household_size == 3)]
      
      results_cum_df = copy(results_df)
      
      results_cum_df[order(date), value := cumsum(value), by = list(scenario, network_structure, intervention, scenario_group, metric, seed)]
      
      plot_df[["new_infections"]] = results_df[, .(mean = mean(value), min = min(value), max = max(value)), 
                                                                by = list(scenario, network_structure, intervention, scenario_group, metric, date)]
      
      plot_df[["cumulative_infections"]] = results_cum_df[, .(mean = mean(value), min = min(value), max = max(value)), 
                                                      by = list(scenario, network_structure, intervention, scenario_group, metric, date)]
      
      plot_df = rbindlist(plot_df, idcol = "metric_detail")
      
      plot_df$metric_detail = factor(plot_df$metric_detail, levels = c("new_infections", "cumulative_infections"))
      
      dict$metric_detail = c(new_infections = "New infections",
                                cumulative_infections = "Cumulative infections")

      network_structure_colors = c("#000000", brewer.pal(3, "Set1")[2])
      
      figure_s5 = ggplot(plot_df, aes(x=date, y=mean, color=network_structure, 
                                      linetype = intervention, fill = network_structure)) +
         annotate(geom = "rect", xmin = 150, xmax = 300, 
                  ymin = 0, ymax = Inf, alpha = 0.1, fill = brewer.pal(3, "Set1")[1]) +
         geom_ribbon(aes(ymin = min, ymax = max), color = NA, alpha = .2) +
         geom_line(size = .8) + 
         scale_x_continuous(name = "Date", limits = c(0, 400), 
                            expand = expansion(mult = c(0,0))) +
         scale_y_continuous(name = "New infections", limits = c(0, NA),
                            expand = expansion(mult = c(0, 0.05)),
                            labels = function(x) format(x, big.mark = "'")) +
         scale_linetype_discrete(name = "Intervention", label = as_labeller(dict$intervention)) +
         scale_color_manual(name = "Network structure", label = as_labeller(dict$network_structure),
                            values = network_structure_colors) +
         scale_fill_manual(name = "Network structure", label = as_labeller(dict$network_structure),
                           values = network_structure_colors) +
         facet_wrap(~ metric_detail, scales = "free", labeller = as_labeller(dict$metric_detail)) +
         guides(color = guide_legend(order = 1),
                fill = guide_legend(order = 1),
                linetype = guide_legend(override.aes = list(fill = NA), order = 2)) +
         theme_bw() +
         theme(
            legend.position = "bottom",
            legend.title = element_text(size = 14, color = "black"),
            legend.text = element_text(size = 14, color = "black"),
            legend.box.margin = margin(0,0,0,0),
            axis.text = element_text(size = 14, color = "black"),
            axis.text.x = element_text(margin = margin(5,0,0,0)),
            axis.text.y = element_text(margin = margin(0, 5, 0, 0)),
            axis.title = element_text(size = 14, color = "black", margin = margin(5, 0, 0, 0)),
            axis.title.x = element_text(margin = margin(10, 0, 0, 0)),
            axis.title.y = element_blank(),
            plot.margin = margin(10,30,10,10),
            strip.background = element_rect(fill = NA, color = NA),
            strip.text = element_text(size = 14, color = "black", margin = margin(0, 0, 10, 0)))
      
      fig_save(o, figure_s5, title = c("Figure S5", "Temporal dynamics", "Intervention"),
               width = 10, height = 6)
      
   }
   
   return(NULL)
   
}  

# ---------------------------------------------------------
# Assess adjusted betas
# ---------------------------------------------------------
assess_betas <- function(o, results, dict) {
   
   # Use result data with traces instead - needed to calculate error of each seed
   trace_data <- load_manuscript_data_with_traces(o)
   
   # All networks should be compared against single-layer data
   comparator_scenario <- "SL_UNI"
   
   # Extract data of comparator scenario
   comparator_df <- results$output[scenario == comparator_scenario, .(scenario, date, metric, mean, median)]
      
   # list of remaining scenarios (ex baseline and comparator)
   analysis_scenarios <- results$output$scenario %>% unique() %>% setdiff(c("baseline", comparator_scenario))
   
   # Data of analysis scenarios
   analysis_df <- trace_data[scenario %in% analysis_scenarios, .(scenario, date, metric, seed, value)]
 
   # Merge comparator scenario into analysis scenarios
   analysis_df <- analysis_df[comparator_df, .(scenario, metric, date, seed, value, comparator_mean = mean, comparator_median = median), on = .(metric, date)]
   
   # Calculate RMSE of time series
   analysis_df[, "sq_res" := (value-comparator_mean)^2]
   
   # Extract contact weights and household sizes from scenario names
   analysis_df[grepl("UNI", scenario), contact_weight := "UNI"]
   analysis_df[grepl("CW1", scenario), contact_weight := "CW1"]
   analysis_df[grepl("CW2", scenario), contact_weight := "CW2"]
   analysis_df[grepl("CW3", scenario), contact_weight := "CW3"]
   analysis_df[grepl("CW4", scenario), contact_weight := "CW4"]
   analysis_df[grepl("CW5", scenario), contact_weight := "CW5"]
   
   analysis_df[grepl("HS0", scenario), household_size := 0]
   analysis_df[grepl("HS1", scenario), household_size := 1]
   analysis_df[grepl("HS2", scenario), household_size := 2]
   analysis_df[grepl("HS3", scenario), household_size := 3]
   analysis_df[grepl("HS4", scenario), household_size := 4]
   analysis_df[grepl("HS5", scenario), household_size := 5]

   # Extract beta parameter label from scenario name
   analysis_df[, scenario_group := gsub("_b[0-9]+$", "", scenario)]
   
   # Prepare data frame for plotting
   plot_df <- analysis_df %>%
      group_by(scenario, metric, seed, contact_weight, household_size, scenario_group) %>%
      summarize(sum_sq_res = sum(sq_res)) %>%
      ungroup() %>%
      group_by(scenario, metric, contact_weight, household_size, scenario_group) %>%
      summarize(rmse = sqrt(mean(sum_sq_res))) %>%
      ungroup() %>%
      as.data.table()
   
   # Merge numerical value of beta into data frame
   plot_df <- plot_df %>%
      left_join(results$stats %>% 
                   filter(metric == "all_new_infections" & variable == "peak_date") %>%
                   select(scenario, beta_value = beta), on = scenario)
   
   analysis_df_mean <- analysis_df %>%
      group_by(scenario, metric, date, scenario_group) %>%
      summarize(mean = mean(value)) %>%
      ungroup() %>%
      as.data.table()
   
   # ---- Plotting of data ----
   for(s in unique(plot_df$scenario_group)) {
      
      this_df <- plot_df[metric == "all_new_infections" & scenario_group == s & !(grepl("b0", scenario))]
      
      # Extract _b identifier
      this_df[, beta_label := str_replace(scenario, paste0(scenario_group, "_"), "")]
   
      # Prepare diagnostic plot for beta values
      diag_beta <- this_df %>%
         ggplot(aes(x=beta_value, y=rmse)) +
         geom_point() +
         geom_point(inherit.aes = FALSE,
                    data = this_df[rmse == min(rmse)], 
                    aes(x=beta_value, y=rmse),
                    color = "red") +
         geom_hline(data = this_df[rmse == min(rmse)],
                     aes(yintercept = rmse),
                    color = "red") + 
         geom_text_repel(aes(label = paste0(beta_label, ": ", beta_value)), size = 4) +
         theme_bw() +
         xlab("Adjusted beta") +
         ylab("RMSE") +
         # ggtitle(paste0(s, " - Adjusted beta")) +
         theme(
            axis.title = element_text(size = 14),
            axis.text = element_text(size = 12, color = "black"),
            plot.title = element_text(size = 14),
            plot.margin = margin(t = 10, r = 25, b = 10, l = 10)
         )
      
      fig_save(o, diag_beta, paste0("Diagnostic beta - ", s), width = 10, height = 6)
      
      # Extract lowest RMSE 
      min_obj <- this_df[rmse == min(rmse)]
   
      this_df_mean <- analysis_df_mean[metric == "all_new_infections" & scenario_group == s & !grepl("b0", scenario)]
      
      this_df_mean[, scenario := factor(scenario)]
      
      # Build color palette for time series by beta value
      color_values <- c("#000000", "#000000", 
                       colorRampPalette(brewer.pal(8, "Set1"))(this_df_mean[, scenario] %>% as.character() %>% unique() %>% length() - 1))
      
      color_labels <- c(comparator_scenario,
                        min_obj$scenario,
                        this_df_mean[scenario != min_obj$scenario, scenario] %>% as.character() %>% unique())
      
      names(color_values) <- color_labels
      
      linetype_values <- c("dashed", rep("solid", this_df_mean[, scenario] %>% as.character() %>% unique() %>% length()))
      
      linetype_labels <- c(comparator_scenario,
                           this_df_mean[, scenario] %>% as.character() %>% unique())
      
      names(linetype_values) <- linetype_labels
      
      linetype_values <- linetype_values[names(color_values)]
      
      # Prepare time series plot by beta
      diag_timeseries <- 
      ggplot(NULL) +
         geom_line(data = this_df_mean[scenario != min_obj$scenario],
                   aes(x=date, y=mean, 
                       color = scenario,
                       linetype = scenario),
                   alpha = .6, 
                   size = .8) +
         geom_line(inherit.aes = FALSE,
                   data = comparator_df[metric == "all_new_infections"],
                   aes(x=date, y=mean, 
                       color = scenario,
                       linetype = scenario)) +
         geom_line(inherit.aes = FALSE,
                   data = this_df_mean[scenario == min_obj$scenario, ],
                   aes(x=date, y=mean, 
                       color = scenario,
                       linetype = scenario)) +
         labs(x = "Date", y = "RMSE", 
              color  = "Simulation", linetype = "Simulation") +
         scale_color_manual(values = color_values) +
         scale_linetype_manual(values = linetype_values) +
         theme_bw() +
         theme(
            axis.title = element_text(size = 14),
            axis.text = element_text(size = 12, color = "black"),
            plot.title = element_text(size = 14),
            legend.text = element_text(size = 12),
            legend.position = "bottom"
         )
      
      # Save time series plot 
      fig_save(o, diag_timeseries, paste0("Diagnostic timeseries - ", s), width = 10, height = 6)
    
   }
   
}

# ---------------------------------------------------------
# Plotting function for Figure 4 (Intervention plots)
# ---------------------------------------------------------
plot_figure_4 <- function(o, results, plot_metric = "all_new_infections", dict, sim_horizon) {
   
   # Define labels for plotting
   metric_label <- qc(cum_value  = "Reduction in cumulative\ninfections/admissions (in %)",
                      peak_value = "Reduction in peak\ninfections/admissions (in %)")
   
   # Vector of all scenarios in the dataset
   all_scenarios <- results$stats$scenario %>% unique() %>% as.character()
   
   # Collate list of intervention scenarios and comparator scenarios
   plot_scenarios <- list(intervention_scenario = all_scenarios[grepl("_WPC", all_scenarios)]) %>%
      as.data.table()
   
   plot_scenarios[, comparator_scenario := sub("_WPC", "", intervention_scenario)]
   
   plot_df <- results$stats_raw[scenario %in% unlist(plot_scenarios) %>% unname() & variable %in% c("cum_value", "peak_value"), ] %>%
      copy()
   
   # Divide into data for interventions ...
   intervention_df <- plot_df[scenario %in% plot_scenarios[, intervention_scenario], ]
   
   intervention_df <- plot_scenarios[intervention_df, on = .(intervention_scenario = scenario) ]
   
   # ... and comparator scenarios
   comparator_df <- plot_df[scenario %in% plot_scenarios[, comparator_scenario], 
                            .(comparator_scenario = scenario, metric, generalized_cc, contact_weight, household_size, variable, comparator_value = value, seed)]
   
   plot_data <- comparator_df[intervention_df, on = .(comparator_scenario, metric, variable, seed)
   ][, .(intervention_scenario, comparator_scenario, generalized_cc, metric, contact_weight, household_size, variable, network_structure, value, comparator_value, seed)]
   
   # Calculate difference in metrics and error bars
   plot_data[, delta := 100*( (comparator_value - value)/comparator_value)]
   
   plot_data[, "altbeta" := "default"]
   plot_data[grepl("altbeta", intervention_scenario), altbeta := "adjusted"]
   
   plot_data = plot_data[, .(mean = mean(delta), lower = quantile(delta, 0.025), upper = quantile(delta, 0.975), min = min(delta), max = max(delta)), 
                         by = list(intervention_scenario, comparator_scenario, generalized_cc, metric, contact_weight, household_size, variable, network_structure, altbeta)]
   
   # Use only data for plotting with alternative beta OR single-layer network
   plot_data_subset <- plot_data[grepl("altbeta", intervention_scenario) | grepl("SL_", intervention_scenario), ]
   
   plot_data_subset = plot_data_subset[contact_weight == "UNI", ]
   
   # Main plot without fitted lines
   fig4 <- ggplot(plot_data_subset[!grepl("SL_", intervention_scenario), ], 
                  aes(x=generalized_cc, y = mean, color = contact_weight, shape = factor(household_size))) +
      facet_grid(variable ~ metric, switch = "y",
                 labeller = labeller(metric = dict$metric, variable = metric_label))
   
   fig4 <- fig4  +
      scale_color_brewer(palette = "Set1", name = "Contact weight", label = dict$contact_weight, guide = "none")
   
   fig4 <- fig4 +
      geom_errorbar(aes(ymin = min, ymax = max), size = .3) +
      scale_shape_discrete(name = "Household size") +
      geom_point(size = 3) +
      geom_point(data = plot_data_subset[grepl("SL_", intervention_scenario)], inherit.aes = FALSE,
                 aes(x = generalized_cc, y = mean), size = 3, shape = 18) +
      geom_errorbar(data = plot_data_subset[grepl("SL_", intervention_scenario)], inherit.aes = FALSE,
                    aes(x = generalized_cc, y = mean, ymin = min, ymax = max), size = .3) +
      geom_hline(aes(yintercept = 0), color = "blue", size = .5, linetype = "dashed") +
      scale_x_continuous(name   = "Generalised clustering coefficient",
                         limits = c(0, NA),
                         expand = expansion(mult = c(0.02, -0.5))) +
      scale_y_continuous(name   = "",
                         labels = function(x) format(x, big.mark = "'", scientific = FALSE),
                         expand = expansion(mult = c(0.05, 0.1)),
                         breaks = seq(0, 80, 20)) +
      theme_bw() +
      theme(legend.position = "bottom",
            legend.text = element_text(size = 14),
            legend.title = element_text(size = 14),
            legend.direction = "horizontal",
            legend.box.just = "left",
            legend.box = "vertical",
            # Facet styles
            strip.placement = "outside",
            strip.background = element_blank(),
            strip.text = element_text(size = 14),
            strip.text.x = element_text(margin = margin(b = 10)),
            # Axes
            axis.title = element_text(size = 14),
            axis.text = element_text(size = 12),
            axis.text.x = element_text(margin = margin(t = 5, b = 10)),
            strip.text.y = element_text(margin = margin(l = 0, r = 10)),
            plot.margin = margin(t = 10, r = 15, b = 10, l = -5, unit = "pt"))
   
   fig4 <- fig4 +
      stat_smooth(method      = "lm", 
                  formula     = y ~ poly(x, 1), 
                  inherit.aes = FALSE, 
                  aes(x = generalized_cc, y = mean), 
                  size        = .5,
                  color       = "black")
   
   ggsave(paste0("Figure 4 - Workplace closure.png"), fig4, path = o$pth$figures, width = 9, height = 8)   
   
   return(NULL)
   
}

# ---------------------------------------------------------
# Plotting function for Figure 3 (Clustering plot)
# ---------------------------------------------------------
plot_figure_3 <- function(o, results, plot_metric, dict, summary_metric = "mean") {

   # ---- Plot preparations ----
   
   # Define axis labels for Figure 3
   plot_labels <- list(
      all_new_infections  = qc(cum_value  = "Cumulative infections\n(per 100'000)",
                               peak_value = "Peak infections\n(per 100'000)",
                               peak_date  = "Time of peak infections\n(in days)"),
      hospital_admissions = qc(cum_value  = "Cumulative admissions\n(per 100'000)",
                               peak_value = "Peak admissions\n(per 100'000)",
                               peak_date  = "Time of peak admissions\n(in days)"))
   
   summary_metric_labels = list(mean   = "Mean",
                             median    = "Median")
   
   # ---- Prepare data for plot ----

   results_df = results$stats[scenario_group == "general" & scenario != "baseline" & metric == plot_metric, 
                        .(scenario, generalized_cc, network_structure, variable, contact_weight, household_size, metric, value = get(summary_metric))]
   
   # Define scenario which scenarios should be compared to
   comparator_scenario <- "SL_UNI"
   
   # Filter data for comparator scenario
   comparator_data <- results_df[scenario == comparator_scenario, .(comparator_scenario = scenario, variable, comparator_value = value)]
   
   # Merge in comparator data
   plot_df <- comparator_data[results_df, on = .(variable)]
   
   # Calculate basic impact metrics
   plot_df[, delta := value - comparator_value]
   
   # Convert some variables to factors (prevents errors from happening in plotting)
   plot_df <- plot_df[, scenario := factor(scenario)]
   
   plot_df[, variable := factor(variable, levels = names(plot_labels[[plot_metric]]))]
   
   # Consider only uniform contact weights and contact weights up to CW4
   plot_df <- plot_df[contact_weight != "CW5"]
   
   # ---- Start plotting ----
   
   plot_stats <- ggplot(plot_df[network_structure == c("multi_layer")], 
                        aes(x=generalized_cc, y=delta, color = contact_weight, shape = as.factor(household_size))) +
      geom_ribbon(stat="smooth", 
                  method = "lm", se=TRUE, alpha=0.1, size = 0, 
                  data = plot_df[network_structure == "multi_layer"], 
                  aes(x = generalized_cc, y = delta, color = contact_weight, fill = contact_weight),
                  formula = y ~ poly(x,1), 
                  inherit.aes = FALSE) +
      stat_smooth(method = "lm", formula = y ~ poly(x, 1), inherit.aes = FALSE, se = FALSE,
                  data = plot_df[network_structure == "multi_layer"], aes(x = generalized_cc, y = delta, color = contact_weight, fill = contact_weight),
                  show.legend = FALSE, alpha = .5) +
      scale_fill_brewer(palette = "Set1", name = "Contact weights", 
                        label = as_labeller(dict$contact_weight)) +
      geom_point(size = 4) + 
      geom_point(inherit.aes = FALSE, data = plot_df[network_structure == "single_layer"], aes(x=generalized_cc, y = delta), size = 4, shape = 18) +
      scale_color_brewer(palette = "Set1", name = "Contact weights", 
                         label = as_labeller(dict$contact_weight)) +
      facet_grid(vars(variable), scales = "free", switch = "y", labeller = as_labeller(plot_labels[[plot_metric]])) +
      scale_shape_manual(name = "Household size", values = c(15,16,17,0,1,2)) +
      scale_x_continuous(name = "Generalised clustering coefficient", 
                         expand = expansion(mult = c(0.02,0.02)),
                         limits = c(0, NA)) +
      scale_y_continuous(labels = function(x) { format(x, big.mark = "'") },
                         expand = expansion(mult = c(0.05, 0.05))) +
      geom_abline(intercept = 0, linetype = "dashed") +
      theme_bw() +
      guides(shape = guide_legend(nrow = 1, override.aes = list(color = "black")),
             color = guide_legend(nrow = 1, override.aes = list(size = 5)),
             fill = guide_legend(nrow = 1, override.aes = list(full = NULL))) +
      theme(legend.position = "bottom",
            legend.direction = "horizontal",
            legend.box.just = "left",
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 14),
            axis.title.y = element_blank(),
            axis.text = element_text(size = 12),
            axis.title.x = element_text(size = 14),
            strip.text.y = element_text(size = 14),
            legend.box = "vertical",
            strip.background = element_blank(),
            strip.placement = "outside") +
      guides(color = guide_legend(override.aes = list(shape = 15)))
   
   fig_save(o, plot_stats, title = c("Figure 3", "Clustering statistics", dict$metric[[plot_metric]], summary_metric_labels[[summary_metric]]), width = 10, height = 8)
   
   return(NULL)
   
}

# ---------------------------------------------------------
# Plotting function for Figure 2 (Temporal plot)
# ---------------------------------------------------------
plot_figure_2_s2_s3 <- function(o, results, plot_metric = "all_new_infections", dict) {
   
   # ---- Define scenarios and all labels for plotting ----
   
   # Define scenarios that should be plotted + assign labels 
   fig2_scenarios <- qc(SL_UNI     = "Single-layer",
                        ML_UNI_HS3 = "Multi-layer (uniform CW)",
                        ML_CW1_HS3 = "Multi-layer (CW1)",
                        ML_CW2_HS3 = "Multi-layer (CW2)",
                        ML_CW3_HS3 = "Multi-layer (CW3)",
                        ML_CW4_HS3 = "Multi-layer (CW4)")
   
   # Sort out n-1 colors (Single-layer results are shown as black lines)
   fig2_colors <- c("#555555", brewer.pal(length(fig2_scenarios)-1, "Set1"))
   
   # Labels for x-axes and y-axes as well as file names
   plot_labels <- list(
      all_new_infections = list(x  = "Date",
                                output      = list(y     = "New infections (per 100'000)",
                                              fname      = "New infections"),
                                cum_output  = list(y     = "Cumulative infections (per 100'000)",
                                                   fname = "Cumulative infections")),
      hospital_admissions = list(x = "Date",
                                 output     = list(y     = "New hospital admissions (per 100'000)",
                                                   fname = "Admissions"),
                                 cum_output = list(y     = "Cumulative hospital admissions (per 100'000)",
                                                  fname  = "Cumulative admissions")))
   
   series_filenames = list(
      all_new_infections = list(output = "Figure 2 - Time series - New infections",
                                cum_output = "Figure S1 - Time series - Cumulative infections"),
      hospital_admissions = list(output = "Figure X - Time series - New hospital admissions",
                                 cum_output = "Figure X - Time series - Cumulative hospital admissions"))
   
   trace_filenames = list(
      all_new_infections = list(output = "Figure S2 - Trace plot - New infections",
                                cum_output = "Figure S3 - Trace plot - Cumulative infections"),
      all_new_infections = list(output = "Figure X - Trace plot - New hospital admissions",
                                cum_output = "Figure X - Trace plot - Cumulative hospital admissions"))
   
   # ---- Prepare data ----
   
   # Extract network data for set of scenarios and output metric
   gcc_labels <- copy(results$stats[scenario %in% names(fig2_scenarios) & metric == plot_metric, 
                                   .(scenario, generalized_cc, variable, mean)])
   
   # Convert wide format to long format
   gcc_labels <- pivot_wider(gcc_labels, names_from = variable, values_from = mean) %>% 
      as.data.table()
   
   plot_set <- c("output", "cum_output")
   
   for(i in plot_set) {
      
      if(i == "output") {
         
         plot_df = results$output_raw[scenario %in% names(fig2_scenarios) & 
                               scenario_group == "general" & 
                               metric == plot_metric]
         
         plot_df = plot_df[, .(mean = mean(value), min = min(value), max = max(value)), by = list(scenario, date)]
         
         
      } else {
         
         plot_df = results$output_raw[scenario %in% names(fig2_scenarios) & 
                                         scenario_group == "general" & 
                                         metric == plot_metric]
         
         plot_df[order(date), value := cumsum(value), by = list(scenario, seed)]
         
         plot_df = plot_df[, .(mean = mean(value), min = min(value), max = max(value)), by = list(scenario, date)]
         
      }
      
      # Convert scenario variable from character to factor
      plot_df[, scenario := factor(scenario, levels = names(fig2_scenarios))]
      
      # Build baseline plot (misses only the label for the GCC)
      plot_series = ggplot(plot_df, aes(x = date, y = mean, color = scenario)) +
         geom_ribbon(aes(ymin = min, ymax = max, fill = scenario), alpha = .2, size = 0) +
         scale_color_manual(values = fig2_colors, name = "Network", label = as_labeller(fig2_scenarios)) +
         scale_fill_manual(values = fig2_colors, name = "Network", label = as_labeller(fig2_scenarios)) +
         scale_y_continuous(name = plot_labels[[plot_metric]][[i]][["y"]], limits = c(0, NA), expand = expansion(mult = c(0, 0.05)),
                            labels = function(x) { format(x, big.mark = "'") }) +
         scale_x_continuous(name = plot_labels[[plot_metric]][["x"]], limits = c(0, 400), expand = expansion(mult = c(0,0))) +
         geom_line(size = .8)
         
      if(i == "output") {
   
         # Add GCC next to the peak value
         plot_series = plot_series + 
            geom_label_repel(data = gcc_labels, 
                             aes(x = peak_date, y=peak_value, label = format(round(generalized_cc, 3), nsmall = 3), color = scenario), 
                             show.legend = FALSE, nudge_x = 30, nudge_y = 0,
                             min.segment.length = 100000, size = 5)
         
      } else {
         
         # Add GCC on the very right hand side
         plot_series = plot_series +
            geom_label_repel(data = gcc_labels, 
                             aes(x=max(plot_df[, .(date)]), y=cum_value, label = format(round(generalized_cc, 3), nsmall = 3), color = scenario), 
                             show.legend = FALSE, nudge_x = 30, nudge_y = 0,
                             min.segment.length = 100000, size = 5)
         
      }
         
      plot_series <- plot_series +
         guides(color=guide_legend(ncol=3,byrow=TRUE)) +
         theme_bw() +
         theme(
            legend.position = "bottom",
            legend.title = element_text(size = 14, color = "black"),
            legend.text = element_text(size = 14, color = "black"),
            legend.box.margin = margin(0,0,0,0),
            axis.text = element_text(size = 14, color = "black"),
            axis.text.x = element_text(margin = margin(5,0,0,0)),
            axis.text.y = element_text(margin = margin(0, 5, 0, 0)),
            axis.title = element_text(size = 14, color = "black", margin = margin(5, 0, 0, 0)),
            axis.title.x = element_text(margin = margin(10, 0, 0, 0)),
            axis.title.y = element_text(margin = margin(0, 10, 0, 0)),
            plot.margin = margin(10,30,10,10))
      
      fig_save(o, plot_series, title = series_filenames[[plot_metric]][[i]], width = 10, height = 8)
      
   }

   # ---- Trace plots ----
   
   # Filter for metric and relevant scenarios
   traces_df <- results$output_raw[metric == plot_metric & scenario %in% names(fig2_scenarios)]
   
   # Convert scenario and seed to factor variables
   traces_df[, scenario := factor(scenario, levels = names(fig2_scenarios))]
   traces_df[, seed := factor(seed)]
   
   # Preallocate empty variables to store data
   traces_list <- traces_labels <- traces_summary <- fig_traces <- list()
   
   # Store data in lists
   traces_list[["output"]] <- copy(traces_df[, .(metric, scenario, seed, date, value)])
   
   traces_list[["cum_output"]] <- copy(traces_df[order(date), .(date, value = cumsum(value)), by = list(metric, scenario, seed)])
   
   # Loop through output and cumulative output
   for(i in plot_set) {
      
      # Extract simulated series by seed
      plot_traces_df = traces_list[[i]]
      
      # Calculate mean over simulated series (for plotting)
      plot_mean_df = traces_list[[i]][, .(mean = mean(value)), by = list(metric, scenario, date)]
      
      # As intermediate step, calculate peak values by seed (necessary for plotting only, so that labels are well positioned)
      plot_traces_peak_df = traces_list[[i]][, .(mean = mean(value)), by = .(metric, scenario, date)] %>%
         group_by(metric, scenario) %>%
         filter(mean == max(mean)) %>%
         slice_max(order_by = -date, n = 1, with_ties = FALSE) %>%
         as.data.table()
      
      # Merge in GCC for each scenario
      traces_labels <- gcc_labels[plot_traces_peak_df, .(metric, scenario, date, mean, generalized_cc), on = .(scenario)]   
      
      # Main plotting command
      plot_traces <- ggplot(plot_traces_df, aes(x=date, y=value, linetype = seed, color = scenario)) +
         scale_color_manual(values = fig2_colors, name = "Network", label = as_labeller(fig2_scenarios)) +
         geom_line(size = .2, alpha = .3, show.legend = FALSE) + # Gives actual trace lines
         scale_linetype_manual(values = rep("solid", length(unique(traces_df$seed)))) + # Workaround to display all seeds; group attribute not working
         scale_y_continuous(name = plot_labels[[plot_metric]][[i]][["y"]], limits = c(0, NA), expand = expansion(mult = c(0, 0.05)),
                            labels = function(x) { format(x, big.mark = "'") }) +
         scale_x_continuous(name = plot_labels[[plot_metric]][["x"]], limits = c(0, 400), expand = expansion(mult = c(0,0))) +
         geom_line(data = plot_mean_df, aes(x=date, y=mean, color = scenario), inherit.aes = FALSE, show.legend = TRUE, size = .8) +
         guides(color=guide_legend(nrow=2,byrow=TRUE, override.aes = list(size = .8, alpha = 1)),
                linetype = "none") +
         geom_label_repel(data = traces_labels,  
                          aes(x=date, y=mean, color = scenario, label = format(round(generalized_cc, 3),nsmall = 3)), inherit.aes = FALSE,
                          nudge_x = 40, nudge_y = 60, min.segment.length = 100000, show.legend = FALSE, size = 5) +
         theme_bw() +
         theme(legend.position = "bottom",
               legend.title = element_text(size = 14, color = "black"),
               legend.text = element_text(size = 14, color = "black"),
               legend.box.margin = margin(0,0,0,0),
               axis.text = element_text(size = 14, color = "black"),
               axis.text.x = element_text(margin = margin(5,0,0,0)),
               axis.text.y = element_text(margin = margin(0, 5, 0, 0)),
               axis.title = element_text(size = 14, color = "black", margin = margin(5, 0, 0, 0)),
               axis.title.x = element_text(margin = margin(10, 0, 0, 0)),
               axis.title.y = element_text(margin = margin(0, 10, 0, 0)),
               plot.margin = margin(10,30,10,10))
      
      # Save plot into results folder
      fig_save(o, plot_traces, title = trace_filenames[[plot_metric]][[i]], width = 10, height = 8)
      
   }
   
   return(NULL)
   
}

# ---------------------------------------------------------
# Add descriptive factors based on scenario ID
# ---------------------------------------------------------
add_factors <- function(input_df, dict) {
   
   # Read CSV file which contains scenario details
   factor_df <- read.csv(paste0(o$pth$figures,"scenario_dict.csv"), sep = ";") %>% 
      as.data.table()
   
   # Convert scenario variable to factors (might prevent errors when merging)
   factor_df[, scenario := factor(scenario)]
   
   input_df[, scenario := factor(scenario)]
   
   # Loop through variables which are provided in dictionary
   for(d in names(dict)) {
      
      # Check whether dictionary label is in scenario variable
      if(d %in% names(factor_df)) {
         factor_df[[d]] <- factor(factor_df[[d]], levels = names(dict[[d]]))
      }
      
   }
   
   # Merge factor level into input dataframe
   output_df <- factor_df[input_df, on = .(scenario)]
   
   return(output_df)
   
}

# ------------------------------------------------------------------
# Load plotting labels (convert data table to list of named vectors)
# ------------------------------------------------------------------
load_dictionary <- function(o) {
   
   # Load label file stored as CSV in figures folder
   dict_df <- read.csv(paste0(o$pth$figures,"plot_dict.csv"), sep = ";") %>%
      as.data.table()
   
   # Preallocate space for output list
   dict_list <- NULL
   
   # Identify unique factors
   identifiers <- unique(dict_df$identifier)
   
   for(i in identifiers) {
      
      # Transfer factor levels
      dict_list[[i]] <- as.character(dict_df[identifier == i, label])
      
      # Add factor labels
      names(dict_list[[i]]) <- dict_df[identifier == i, level]
      
   }
   
   return(dict_list)
   
}

# ---------------------------------------------------------
# Read data files for cluster plotting
# ---------------------------------------------------------
load_manuscript_data <- function(o, metrics = c("all_new_infections", "hospital_admissions")) {
   
   # Check whether manuscript data file exists already
   if(file.exists(paste0(o$pth$figures, "manuscript_data.rds"))) {
      
      # If it exists, read data file
      results <- readRDS(paste0(o$pth$figures, "manuscript_data.rds"))
      
   } else {
      
      # Extract baseline yaml input
      baseline_yaml <- parse_yaml(o, "baseline")$parsed
      
      # Calculate population scaler to convert to per 100'000 figures
      pop_scaler <- 100000/baseline_yaml$population_size
      
      # Load list of simulations to extract scenario names
      sims <- try_load(o$pth$simulations, "all_simulations")
      
      # Extract list of scenarios
      scenarios = sims$scenario %>% unique()
      
      # Initialize progress bar
      pb = start_progress_bar(length(scenarios))
      i <- 1
      
      # Initialize some output lists
      output_raw_list = output_list = cum_output_list = network_list = stats_raw_list = stats_list = list()
      
      for(scenario in scenarios) {
         
         # Read in raw data by scenario
         results_raw_df <- try_load(o$pth$scenarios, paste0(scenario, "_raw"))[metric %in% metrics & grouping == "none", .(metric, date, scenario, seed, value)]
         
         output_raw_list[[scenario]] = results_raw_df
         
         # ---- Peak statistics ----
         
         # Extract data on peak value of metric
         peak_data = results_raw_df[results_raw_df[order(date), .I[value == max(value)], by=.(metric, scenario, seed)]$V1][order(date), .SD, by = .(metric, scenario, seed)]
         
         # Use first occurrence of peak value and rename columns
         peak_data = peak_data[peak_data[, .I[1], by = .(metric, scenario, seed)]$V1][order(metric, seed), .(metric, scenario, seed, peak_date = date, peak_value = value)]
         
         # ---- Cumulative values ----
         
         # Read simulation results
         cum_results = copy(results_raw_df)
         
         # Calculate cumulative values by metric and seed
         cum_results[order(date), cum_value := cumsum(value), by = .(metric, scenario, seed)]
         
         # Focus on last date captured in the simulation
         cum_results = cum_results[date == max(date), .SD, by = .(metric, scenario, seed)]
         
         # Compile full table of summary statistics
         stats_raw <- cum_results[peak_data, on = .(metric, scenario, seed)][, .(metric, scenario, seed, peak_date, peak_value, cum_value)]
         
         stats_raw <- pivot_longer(stats_raw, c("peak_date", "peak_value", "cum_value"), names_to = "variable") %>%
            as.data.table()
         
         stats_raw_list[[scenario]] = stats_raw
         
         # ---- Aggregated data (not raw file) ----
         
         results_df <- try_load(o$pth$scenarios, scenario)
         
         output_list[[scenario]] <- results_df$output[metric %in% metrics & grouping == "none"]
         
         cum_output_list[[scenario]] <- results_df$cum_output[metric %in% metrics & grouping == "none"]
         
         # Extract network summary statistics
         network_list[[scenario]] <- list(scenario       = scenario,
                                          cc             = results_df$clustering$cc,
                                          generalized_cc = results_df$clustering$cc_generalized,
                                          beta           = results_df$input$beta)

         stats_list[[scenario]] = results_df$clustering$stats
         
         setTxtProgressBar(pb, i)
         
         i <- i + 1
         
      }
      
      close(pb)
      
      # --- Consolidate stored data ----
      
      output = rbindlist(output_list)
      output_raw = rbindlist(output_raw_list)
      cum_output <- rbindlist(cum_output_list)
      network <- rbindlist(network_list)
      stats <- rbindlist(stats_list)
      stats_raw <- rbindlist(stats_raw_list)
      
      results <- list(output     = output,
                      output_raw = output_raw,
                      cum_output = cum_output,
                      network    = network,
                      stats      = stats,
                      stats_raw  = stats_raw)
      
      # ---- Scaling of statistics ----
      
      # First, raw statistics:
      variable_details <- data.table(variable = c("peak_date", "peak_value", "cum_value"),
                                     scale = c(FALSE, TRUE, TRUE))
      
      # Columns that are scaled with population scaler
      scalable_columns <- c("value")
      
      # Merge scaling details into stats
      results$stats_raw <- variable_details[results$stats_raw, on = .(variable)]
      
      # Multiply relevant columns by scaler
      results$stats_raw <- results$stats_raw %>%
         mutate(across(all_of(scalable_columns), function(x) ifelse(scale, x * pop_scaler, x))) %>%
         select(-scale)
      
      # Second, aggregated statistics
      
      # Columns that are scaled with population scaler
      scalable_columns <- c("mean", "median", "lower", "upper")
      
      # Merge scaling details into stats
      results$stats <- variable_details[results$stats, on = .(variable)]
      
      # Multiply relevant columns by scaler
      results$stats <- results$stats %>%
         mutate(across(all_of(scalable_columns), function(x) ifelse(scale, x * pop_scaler, x))) %>%
         select(-scale)
      
      # ---- Scaling of raw output ----
      
      # Multiply values by population scaler
      results$output_raw <- results$output_raw %>%
         mutate(value = value * pop_scaler)
      
      # ---- Scaling of aggregated data frames ----
      
      # In which list elements should columns be scaled
      idx <- c("output", "cum_output")
      
      # Which columns should be scaled
      scalable_columns <- c("mean", "median", "lower", "upper")
      
      for(i in idx) {
         
         results[[i]] <- results[[i]][metric %in% metrics & grouping == "none"][, .(date, metric, scenario, mean, median, lower, upper)]
         
         results[[i]] <- results[[i]] %>%
            mutate_at(scalable_columns, function(x) x * pop_scaler)
         
      }
      
      # Save manuscript data files
      saveRDS(results, paste0(o$pth$figures, "manuscript_data.rds"))
      
   }
   
   return(results)
   
} 

# ---------------------------------------------------------
# Read data files for cluster plotting
# ---------------------------------------------------------
load_manuscript_data_with_traces <- function(o, metrics = c("all_new_infections", "hospital_admissions")) {
   
   # Check whether manuscript data file exists already
   if(file.exists(paste0(o$pth$figures, "manuscript_data_traces.rds"))) {
      
      # If it exists, read in data
      output <- readRDS(paste0(o$pth$figures, "manuscript_data_traces.rds"))
      
   } else {
      
      # Load list of simulations to extract scenario names
      scenarios <- try_load(o$pth$simulations, "all_simulations")[, scenario] %>% 
         unique()
      
      # Pre-allocate space for result lists
      output_list <- NULL
      
      pb = start_progress_bar(length(scenarios))
      
      i <- 1
      
      # Iterate through all scenario files
      for(s in scenarios) {
         
         # Read in summary data by scenario and seed
         this_result <- try_load(o$pth$scenarios, paste0(s, "_raw"))[metric %in% metrics & grouping == "none", 
                                                                 .(metric, date, scenario, seed, value)]
         
         this_result[order(date), cumulative_value := cumsum(value), by = .(metric, scenario, seed)]
         
         output_list[[s]] <- copy(this_result)
         
         setTxtProgressBar(pb, i)
         
         i <- i + 1
      }
      
      close(pb)
      
      # Consolidate data and prepare for export
      output <- rbindlist(output_list)
      
      # Extract baseline yaml input
      baseline_yaml <- parse_yaml(o, "baseline")$parsed
      
      # Simulation horizon
      sim_horizon <- baseline_yaml$n_days
      
      # Calculate population scaler to convert to per 100'000 figures
      pop_scaler <- 100000/baseline_yaml$population_size
      
      # Columns that are scaled with population scaler
      scalable_columns <- c("value", "cumulative_value")
      
      # Multiply relevant columns by scaler
      output <- output %>% 
         mutate_at(scalable_columns, function(x) x * pop_scaler)
      
      # Save manuscript data files
      saveRDS(output, paste0(o$pth$figures, "manuscript_data_traces.rds"))
      
   }
   
   return(output)
   
} 
