###########################################################
# MEDIA IMAGE
#
# One nice tidy subplot for media purposes.
#
###########################################################

# ---------------------------------------------------------
# Plot media image
# ---------------------------------------------------------
run_media_image = function(o) {
  
  # Only continue if specified by do_step (or forced)
  if (!is.element(4, o$do_step)) return()
  
  message("* Plotting media image")
  
  # These results are for 'omicron' analysis
  if (o$analysis_name != "omicron")
    stop("\nTo produce this figure, you must run steps 1-3 for 'omicron' analysis.", 
         "\nNote that step 2 will require connection to a cluster")
  
  # ---- Figure properties ----
  
  # Name of heatmap to work with
  heatmap_name = c("Heatmap", "ICU beds", "Analysis 2")
  
  # Just a subset of severities will do
  severities = c("Severity=0.5", "Severity=1", "Severity=2")
  
  # Better names for layman figure
  sev_names = c("Less severe than Delta", "Equally severe", "More severe than Delta")
  
  # Colour bin values
  icu_bins = seq(0, 30, by = 5)
  
  # ---- Generate heatmap ----
  
  # Load previously saved heatmap from file
  fig_pth = heatmap_id(o, heatmap_name)
  plot_df = readRDS(fig_pth)$plot_df %>%
    filter(w %in% severities)  %>%
    mutate(w = recode(w, !!!setNames(sev_names, severities)), 
           w = factor(w, levels = rev(sev_names)))
  
  # Construct suitable colour bar labels
  icu_labels    = rev(icu_bins)[-1]
  icu_labels[1] = paste0(">", icu_labels[1])
  
  # No label for max value
  icu_labels = c(rev(icu_labels), "")
  
  # Generate heat map from subsetted plot dataframe
  plot_info = 
    plot_heatmap(o, "Media image", NULL,
                 plot_df        = plot_df,
                 plot_metrics   = "icu_beds",
                 colour_palette = "-Greys", 
                 colour_limits  = c(min(icu_bins), max(icu_bins)),
                 colour_bins    = icu_bins, 
                 x_lab          = "Infectivity", 
                 y_lab          = "Immune evading", 
                 x_scale_fn     = label_number(accuracy = 0.1),
                 y_scale_fn     = percent,
                 z_labels       = icu_labels, 
                 plot_title     = "Peak COVID-19-related ICU occupancy", 
                 plot_subtitle  = "(per 100,000 people)", 
                 n_wrap         = 16)
  
  # Set text sizes for the media image
  g = plot_info$g + theme(
    plot.title       = element_text(size = 35, hjust = 0.5), 
    plot.subtitle    = element_text(size = 20, hjust = 0.5), 
    strip.background = element_blank(),
    strip.text.y     = element_text(size = 20),
    strip.text.x     = element_blank(),
    axis.title       = element_text(size = 35),
    axis.text        = element_text(size = 14),
    panel.spacing.y  = unit(1.5, "lines"))
  
  # Save to file
  fig_save(o, g, "Media image", width = 11, height = 11)
}

