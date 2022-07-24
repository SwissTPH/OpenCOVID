###########################################################
# MANUSCRIPT FIGURES
#
# Results for Omicron manuscript (December 2021).
#
# Written by A.J.Shattock (andrewjames.shattock@swisstph.ch)
###########################################################

# ---------------------------------------------------------
# Main plot-calling function
# ---------------------------------------------------------
manuscript_figures = function(o) {
  
  message(" - All figures for Omicron manuscript")
  
  # These results are for 'omicron' analysis
  if (o$analysis_name != "omicron")
    stop("\nTo produce these results, you must run steps 1-2 for 'omicron' analysis.", 
         "\nNote that step 2 will require connection to a cluster")
  
  # ---- Plotting flags ----
  
  # Set flags for what to plot
  do_temporal      = TRUE  # Temporal plots
  do_n_doses       = TRUE  # Number of doses
  do_heat_variant  = TRUE  # Heat maps of Omicron prevalence
  do_heat_hosp     = TRUE  # Heat maps for peak cases in hospital
  do_heat_metrics  = TRUE  # Heat maps of cases & deaths
  do_heat_compare  = TRUE  # Heat maps comparing these metrics
  
  # Force re-create expensive heatmaps
  force_replot = FALSE
  
  # ---- Figure options ----
  
  # Full names for the two analyses
  analysis_names = c("No future vaccination", 
                     "Expanded vaccination")
  
  # Full names for cases/deaths averted plots
  avert_names = c("Infections averted", 
                  "Deaths averted")
  
  # Define some key scenarios so we can reference by name
  s0 = s1 = list(base = ".inf11.esc1.sev3",
                 inf1 = ".inf21.esc1.sev3",
                 inf2 = ".inf21.esc1.sev5",
                 esc  = ".inf11.esc21.sev3",
                 sev1 = ".inf13.esc5.sev4",
                 sev2 = ".inf16.esc11.sev5",
                 high = ".inf21.esc21.sev5")
  
  # Start and end days of analysis
  plot_from = 61
  plot_to   = 240
  
  # Initaite date for x-axis of temporal plots
  date_from = "2021-12-01"
  
  # Facet indices for Delta/Omicron dominance text                   
  facet_text = c(5, 2) # (row, column)
  
  # Define colourbar bins for hospital peak
  hosp_bins  = seq(0, 350, by = 60)
  avert_bins = seq(0, 0.6, by = 0.1)
  
  # Number of points (across each axis) for interpolation of heatmap 
  heatmap_interpolate1 = 100
  heatmap_interpolate2 = 1000
  
  # Heatmap fontsize: title, subtitle, axes, x-strip, y-strip, legend, ticks
  heatmap_fontsize = c(24, 20, 24, 18, 15, 13, 12)
  
  # Number of characters per line in heatmap strips
  heatmap_wrap = 30
  
  # ---- Array set up ----
  
  # Read all scenario names and remove baseline
  all_arrays = parse_yaml(o, "*read*")[-1]
  n_arrays   = length(all_arrays)
  
  # Analysis file should now contain only the two arrays
  if (n_arrays != 2)
    stop("Analysis file should contain exactly 2 array scenarios")
  
  # Format IDs and names into list for easy access
  a = list(id = names(all_arrays),
           name = unname(all_arrays))
  
  # Dictionary for full analysis names
  analysis_dict = setNames(analysis_names, a$name)
  
  # ---- Temporal plots ----
  
  # Check plotting flag
  if (do_temporal == TRUE) {
    
    message("  > Temporal plots")
    
    # Variant prevalence over time
    plot_temporal(o, "Delta dynamics", 
                  plot_baseline = FALSE,
                  scenarios     = paste0(a$id, s0$base), 
                  plot_metrics  = c("all_new_infections", "hospital_beds", "deaths", 
                                    "seasonality", "pop_susceptibility"), 
                  plot_from     = plot_from, 
                  date_from     = date_from)
    
    # Loop through arrays
    for (i in 1 : n_arrays) {
      
      # Add array ID to scenario ID for easy referencing
      s1[] = paste0(a$id[i], s0)
      
      # Variant prevalence over time
      plot_temporal(o, c("Temporal", a$name[i]), 
                    plot_baseline = FALSE, 
                    scenarios     = c(s1$base, s1$inf1, s1$inf2, s1$esc, s1$sev1, s1$sev2, s1$high),
                    plot_metrics  = c("all_new_infections", "hospital_beds", "variant_prevalence"),
                    plot_by       = "variant", 
                    plot_geom     = "area",
                    plot_from     = plot_from, 
                    date_from     = date_from, 
                    n_wrap        = 20, 
                    scenario_name_fn  = scenario_name_fn, 
                    override_colours  = c("#D35495", "#4D9221"), # Extremes of PiYG palette
                    override_fontsize = c(NA, 12, 14, 12)) # title, facets, legend, ticks
    }
  }
  
  # ---- Number of doses ----
  
  # Only do number of doses for analysis 2 (with boosters)
  if (do_n_doses == TRUE) {
    
    message("  > Number of doses")
    
    # Add array ID to scenario ID for easy referencing
    s1[] = paste0(a$id[2], s0)
    
    # Cumulative doses over time for key scenarios
    fig_name = c("Number of doses", "Area", a$name[2])
    g = plot_temporal(o, fig_name,
                      alt_baseline  = s1$base,
                      plot_metrics  = "n_doses",
                      cumulative    = TRUE,
                      plot_by       = "vaccine_group",
                      plot_geom     = "area",
                      plot_from     = plot_from, 
                      date_from     = date_from)
    
    # Remove strip text and add a title
    g = g + labs(title    = "Number of vaccine doses administered", 
                 subtitle = "(per 100,000 people)") + 
      theme(plot.subtitle = element_text(size = 24, hjust = 0.5), 
            strip.text    = element_blank())
    
    # Save the updated figure with a smaller height
    fig_save(o, g, fig_name)
  }
  
  # ---- Heat maps: variant prevalece ----
  
  message("  > Extracting variant prevalence details")
  
  # Initiate list to store variant prevalence contour info
  vp_contour = vp_plots = list()
  
  # Loop through arrays
  for (i in 1 : n_arrays) {
    
    # Variant prevalence figure name and path to plot info
    fig_name = c("Variant prevalence", a$name[i])
    vp_pth   = heatmap_id(o, fig_name)
    
    # Check plotting flag
    if (do_heat_variant == TRUE) {
      
      # Load plotting info is appropriate
      if (!file.exists(vp_pth) || force_replot == TRUE) {
        
        # Heat map of moving parts: Omicron prevalence over analysis period
        plot_heatmap(o, fig_name, a$id[i],
                     plot_metrics   = "variant_prevalence",
                     plot_by        = "variant",
                     group_idx      = 2,  # 1 := Delta, 2:= Omicron
                     summarise      = "max",
                     plot_from      = plot_from,
                     plot_to        = plot_to,
                     colour_limits  = c(0, 100),
                     n_interpolate  = heatmap_interpolate1, 
                     contour_val    = 50,             # Countours at 50% Omicron prevalence
                     contour_fn     = "y ~ exp(-x)",  # ... using exponential decay curves
                     contour_by     = "v",            # Different linetypes for each analysis
                     w_percent      = FALSE,
                     v_percent      = TRUE)
      }
      
      # Store plot dataframe in list
      vp_plots[[i]] = readRDS(vp_pth)$plot_df %>%
        mutate(v = a$name[i])
    }
    
    # If variant prevalence file has been created, extract the countour dataframe
    if (file.exists(vp_pth)) {
      vp_contour[[i]] = readRDS(vp_pth)[qc(contour_df, contour_fn)]
      
      # Append analysis name as 4th dim
      vp_contour[[i]]$contour_df$v = analysis_names[i]
    }
  }
  
  # Check plotting flag
  if (do_heat_variant == TRUE) {
    
    message("  > Plotting variant prevalence")
    
    # Manuscript figure number
    fig_name = "Variant prevalence"
    
    # Combine plotting dataframes
    plot_df = rbindlist(vp_plots) %>%
      mutate(w = factor(w, levels = rev(levels(w))), 
             v = unname(analysis_dict[v]), 
             v = factor(v, levels = rev(analysis_names)))
    
    # Plot manuscript figure
    plot_info =
      plot_heatmap(o, fig_name, NULL,
                   plot_df        = plot_df,
                   plot_metrics   = "variant_prevalence",
                   plot_by        = "variant",
                   colour_palette = "-PiYG",  # Pink -> Green
                   colour_limits  = c(0, 100),
                   contour_val    = 50,
                   contour_fn     = "y ~ exp(-x)",
                   contour_by     = "v",
                   x_lab          = "Infectivity", 
                   y_lab          = "Immune evading", 
                   x_scale_fn     = label_number(accuracy = 0.1),
                   y_scale_fn     = percent,
                   z_labels       = function(x) paste0(x, "%"),
                   plot_title     = "Omicron variant prevalence (after 6 months)",
                   n_wrap         = heatmap_wrap,
                   override_fontsize = heatmap_fontsize)
    
    # Add suplot labels to figure
    add_subplot_labels(plot_info$g, plot_info$plot_df, fig_name)
  }
  
  # ---- Heat maps: peak of hospital cases ----
  
  # Check plotting flag
  if (do_heat_hosp == TRUE) {
    
    message("  > Plotting hospital peak")
    
    # Initiate list to store variant prevalence contour info
    plot_list = list()
    
    # Loop through arrays
    for (i in 1 : n_arrays) {
      
      # Define figure name
      fig_name = c("Heatmap", "Hospitalisations", a$name[i])
      fig_pth  = heatmap_id(o, fig_name)
      
      # Create new plotting dataframe only if needed
      if (!file.exists(fig_pth) || force_replot) {
        
        # Heat map of moving parts: max hospital beds
        plot_heatmap(o, fig_name, a$id[i],
                     plot_metrics   = c("hospital_beds", "icu_beds"),
                     summarise      = "max",
                     plot_from      = plot_from,
                     plot_to        = plot_to,
                     colour_limits  = c(min(hosp_bins), max(hosp_bins)), 
                     n_interpolate  = heatmap_interpolate2,
                     w_percent      = FALSE,
                     v_percent      = TRUE)
      }
      
      # Store plot dataframe in list
      plot_list[[i]] = readRDS(fig_pth)$plot_df %>%
        mutate(v = a$name[i])
    }
    
    # Manuscript figure number
    fig_name = "Figure 1"
    
    # Combine plotting dataframes
    plot_df = rbindlist(plot_list) %>%
      mutate(w = factor(w, levels = rev(levels(w))), 
             v = unname(analysis_dict[v]), 
             v = factor(v, levels = rev(analysis_names)))
    
    # Construct suitable colour bar labels
    hosp_labels    = rev(hosp_labels)[-1]
    hosp_labels[1] = paste0(">", hosp_labels[1])
    
    # No label for max value
    hosp_labels = c(rev(hosp_labels), "")
    
    # Plot manuscript figure
    plot_info =
      plot_heatmap(o, fig_name, NULL,
                   plot_df        = plot_df,
                   plot_metrics   = c("hospital_beds", "icu_beds"), 
                   colour_palette = "-Greys", 
                   colour_limits  = c(min(hosp_bins), max(hosp_bins)),
                   colour_bins    = hosp_bins, 
                   x_lab          = "Infectivity", 
                   y_lab          = "Immune evading", 
                   x_scale_fn     = label_number(accuracy = 0.1),
                   y_scale_fn     = percent,
                   z_labels       = hosp_labels, 
                   plot_title     = "Peak hospitalisation occupancy (per 100,000 people)", 
                   n_wrap         = heatmap_wrap,
                   override_fontsize = heatmap_fontsize)
    
    # Create label dataframe - letter for each facet
    text_df = expand_grid(x = c(0.5, 1.5), y = 0.25, 
                          w = levels(plot_df$w)[facet_text[1]], 
                          v = levels(plot_df$v)[facet_text[2]]) %>%
      mutate(w = factor(w, levels = levels(plot_df$w)), 
             v = factor(v, levels = levels(plot_df$v)),         
             text = c("Delta dominant", "Omicron dominant")) %>%
      as.data.table()
    
    # Add text to existing figure
    plot_info$g = plot_info$g + 
      geom_text(data = text_df, mapping = aes(label = text), size = 4, 
                hjust = "middle", vjust = "centre", colour = "darkred")
    
    # Append variant prevalence contour lines (and subplot labels)
    add_vp_contour(o, vp_contour, plot_info, fig_name, colour = "darkred")
  }
  
  # ---- Heat maps: cases and deaths ----
  
  # Check plotting flag
  if (do_heat_metrics == TRUE) {
    
    # Plot cases ...
    
    message("  > Plotting cumulative cases")
    
    # Initiate list to store variant prevalence contour info
    plot_list = list()
    
    # Loop through arrays
    for (i in 1 : n_arrays) {
      
      # Define figure name
      fig_name = c("Heatmap", "Cases", a$name[i])
      fig_pth  = heatmap_id(o, fig_name)
      
      # Create new plotting dataframe only if needed
      if (!file.exists(fig_pth) || force_replot) {
        
        # Heat map of moving parts: cumulative cases
        plot_heatmap(o, fig_name, a$id[i],
                     plot_metrics   = "all_new_infections",
                     summarise      = "sum",
                     plot_from      = plot_from,
                     plot_to        = plot_to,
                     colour_limits  = c(0, NA),
                     n_interpolate  = heatmap_interpolate1,
                     w_percent      = FALSE,
                     v_percent      = TRUE)
      }
      
      # Store plot dataframe in list
      plot_list[[i]] = readRDS(fig_pth)$plot_df %>%
        mutate(v = a$name[i])
    }
    
    # Manuscript figure number
    fig_name = "Cumulative cases"
    
    # Combine plotting dataframes
    plot_df = rbindlist(plot_list) %>%
      mutate(w = factor(w, levels = rev(levels(w))), 
             v = unname(analysis_dict[v]), 
             v = factor(v, levels = rev(analysis_names)))
    
    # Plot manuscript figure
    plot_info =
      plot_heatmap(o, fig_name, NULL,
                   plot_df        = plot_df,
                   plot_metrics   = "all_new_infections",
                   colour_limits  = c(0, NA),
                   x_lab          = "Infectivity", 
                   y_lab          = "Immune evading",  
                   x_scale_fn     = label_number(accuracy = 0.1),
                   y_scale_fn     = percent,
                   plot_title     = "SARS-CoV-2 infections",
                   plot_subtitle  = "(per 100,000 people over 6 months)", 
                   n_wrap         = heatmap_wrap,
                   override_fontsize = heatmap_fontsize)
    
    # Append variant prevalence contour line
    add_vp_contour(o, vp_contour, plot_info, fig_name)
    
    # Plot deaths ...
    
    message("  > Plotting cumulative deaths")
    
    # Initiate list to store variant prevalence contour info
    plot_list = list()
    
    # Loop through arrays
    for (i in 1 : n_arrays) {
      
      # Define figure name
      fig_name = c("Heatmap", "Deaths", a$name[i])
      fig_pth  = heatmap_id(o, fig_name)
      
      # Create new plotting dataframe only if needed
      if (!file.exists(fig_pth) || force_replot) {
        
        # Heat map of moving parts: cumulative deaths
        plot_heatmap(o, fig_name, a$id[i],
                     plot_metrics   = "deaths",
                     summarise      = "sum",
                     plot_from      = plot_from,
                     plot_to        = plot_to,
                     colour_limits  = c(0, NA),
                     n_interpolate  = heatmap_interpolate1,
                     w_percent      = FALSE,
                     v_percent      = TRUE)
      }
      
      # Store plot dataframe in list
      plot_list[[i]] = readRDS(fig_pth)$plot_df %>%
        mutate(v = a$name[i])
    }
    
    # Manuscript figure number
    fig_name = "Cumulative deaths"
    
    # Combine plotting dataframes
    plot_df = rbindlist(plot_list) %>%
      mutate(w = factor(w, levels = rev(levels(w))), 
             v = unname(analysis_dict[v]), 
             v = factor(v, levels = rev(analysis_names)))
    
    # Plot manuscript figure
    plot_info =
      plot_heatmap(o, fig_name, NULL,
                   plot_df        = plot_df,
                   plot_metrics   = "deaths",
                   colour_limits  = c(0, NA),
                   x_lab          = "Infectivity", 
                   y_lab          = "Immune evading", 
                   x_scale_fn     = label_number(accuracy = 0.1),
                   y_scale_fn     = percent,
                   plot_title     = "COVID-19-related deaths",
                   plot_subtitle  = "(per 100,000 people over 6 months)", 
                   n_wrap         = heatmap_wrap,
                   override_fontsize = heatmap_fontsize)
    
    # Append variant prevalence contour line
    add_vp_contour(o, vp_contour, plot_info, fig_name)
  }
  
  # ---- Comparison heat maps: cases, deaths, and ICU ----
  
  # Check plotting flag
  if (do_heat_compare == TRUE) {
    
    message("  > Plotting cases and deaths averted")
    
    # Initiate list to store variant prevalence contour info
    plot_list = list()
    
    # Plot cases...
    
    # Define figure name
    fig_name = c("Comparison heatmap", "Cases")
    fig_pth  = heatmap_id(o, fig_name)
    
    # Create new plotting dataframe only if needed
    if (!file.exists(fig_pth) || force_replot) {
      
      # Comparison: cumulative cases
      plot_heatmap(o, fig_name, a$id[1],
                   dif_array      = a$id[2],
                   dif_relative   = TRUE,
                   plot_metrics   = "all_new_infections",
                   summarise      = "sum",
                   plot_from      = plot_from,
                   plot_to        = plot_to,
                   colour_limits  = c(0, 1),
                   n_interpolate  = heatmap_interpolate2,
                   w_percent      = FALSE,
                   v_percent      = TRUE)
    }
    
    # Store plot dataframe in list
    plot_list[[1]] = readRDS(fig_pth)$plot_df %>%
      mutate(v = avert_names[1])
    
    # Plot deaths...
    
    # Define figure name
    fig_name = c("Comparison heatmap", "Deaths")
    fig_pth  = heatmap_id(o, fig_name)
    
    # Create new plotting dataframe only if needed
    if (!file.exists(fig_pth) || force_replot) {
      
      # Comparison: cumulative deaths
      plot_heatmap(o, fig_name, a$id[1],
                   dif_array      = a$id[2],
                   dif_relative   = TRUE,
                   plot_metrics   = "deaths",
                   summarise      = "sum",
                   plot_from      = plot_from,
                   plot_to        = plot_to,
                   colour_limits  = c(0, 1),
                   n_interpolate  = heatmap_interpolate2,
                   w_percent      = FALSE,
                   v_percent      = TRUE)
    }
    
    # Store plot dataframe in list
    plot_list[[2]] = readRDS(fig_pth)$plot_df %>%
      mutate(v = avert_names[2])
    
    # Combine cases and deaths...
    
    # Manuscript figure number
    fig_name = "Figure 2"
    
    # Combine plotting dataframes
    plot_df = rbindlist(plot_list) %>%
      mutate(w = factor(w, levels = rev(levels(w))), 
             v = factor(v, levels = avert_names))
    
    # Construct suitable colour bar labels
    avert_labels    = paste0(rev(avert_bins)[-1] * 100, "%")
    avert_labels[1] = paste0(">", avert_labels[1])
    
    # No label for max value
    avert_labels = c(rev(avert_labels), "")
    
    # Combined plots: cases and deaths averted
    plot_info =
      plot_heatmap(o, fig_name, NULL,
                   plot_df        = plot_df,
                   plot_metrics   = "all_new_infections",
                   colour_palette = "-Spectral",
                   colour_limits  = c(min(avert_bins), max(avert_bins)),
                   colour_bins    = avert_bins,
                   x_lab          = "Infectivity", 
                   y_lab          = "Immune evading", 
                   x_scale_fn     = label_number(accuracy = 0.1),
                   y_scale_fn     = percent,
                   z_labels       = avert_labels,
                   n_wrap         = 50,
                   plot_title     = "Infections and deaths averted by expanded vaccination", 
                   override_fontsize = heatmap_fontsize)
    
    # Create label dataframe - letter for each facet
    text_df = expand_grid(x = c(0.5, 1.5), y = 0.25, 
                          w = levels(plot_df$w)[facet_text[1]], 
                          v = levels(plot_df$v)[facet_text[2]]) %>%
      mutate(w = factor(w, levels = levels(plot_df$w)), 
             v = factor(v, levels = levels(plot_df$v)),         
             text = c("Delta dominant", "Omicron dominant")) %>%
      as.data.table()
    
    # Add text to existing figure
    plot_info$g = plot_info$g + 
      geom_text(data = text_df, mapping = aes(label = text), size = 4, 
                hjust = "middle", vjust = "centre", colour = "black")
    
    # A bit of additonal work required to get both contours on all facets
    for (i in 1 : length(vp_contour)) {
      this_df = vp_contour[[i]]$contour_df
      
      # Set 4th dimension to be metrics names 
      vp_contour[[i]]$contour_df =
        rbind(mutate(this_df, v = factor(avert_names[1], levels = avert_names)), 
              mutate(this_df, v = factor(avert_names[2], levels = avert_names)))
    }
    
    # Append second variant prevalence contour for dif array - as dashed line
    add_vp_contour(o, vp_contour, plot_info, fig_name)
  }
}

# ---------------------------------------------------------
# Add Omicron domination threshold contour to existing plot
# ---------------------------------------------------------
add_vp_contour = function(o, vp_contour, plot_info, fig_name, colour = "black") {
  
  # Throw a warning in no variant prevalence info
  if (length(vp_contour) == 0) {
    warning("No variant prevalence contour to add")
    
    return()  # Return out
  }
  
  # Extract ggplot
  g = plot_info$g
  
  # Loop through contours to add
  for (i in seq_along(vp_contour)) {
    
    # Filter for severities which have been plotted
    this_contour = vp_contour[[i]]$contour_df %>%
      filter(w %in% unique(plot_info$plot_df$w))
    
    # Append contour to original plot
    g = g + geom_smooth(data      = this_contour, 
                        mapping   = aes(x = x, y = y), 
                        method    = "lm", 
                        formula   = vp_contour[[i]]$contour_fn, 
                        fullrange = TRUE, 
                        se        = FALSE, 
                        linetype  = length(vp_contour) - i + 1,
                        size      = 2, 
                        colour    = colour)
  }
  
  # Extract limits of data
  data_lims = list(x_min = min(plot_info$plot_df$x), 
                   x_max = max(plot_info$plot_df$x), 
                   y_min = min(plot_info$plot_df$y), 
                   y_max = max(plot_info$plot_df$y))
  
  # Setting limits within scale_x_continuous doesn't work - this is a workaround
  g = g + coord_cartesian(xlim = c(data_lims$x_min, data_lims$x_max), 
                          ylim = c(data_lims$y_min, data_lims$y_max))
  
  # Add suplot labels to figure
  add_subplot_labels(g, plot_info$plot_df, fig_name)
}

# ---------------------------------------------------------
# Add subplot labels to existing plot
# ---------------------------------------------------------
add_subplot_labels = function(g, plot_df, fig_name, x_gap = 0.03, y_yap = 0.05) {
  
  # Determine x and y positions from extremes of plotting data
  x_val = min(plot_df$x) + (max(plot_df$x) - min(plot_df$x)) * x_gap
  y_val = max(plot_df$y) - (max(plot_df$y) - min(plot_df$y)) * y_yap
  
  # Create label dataframe - letter for each facet
  label_df = expand_grid(x = x_val, y = y_val, 
                         w = levels(plot_df$w), 
                         v = levels(plot_df$v)) %>%
    mutate(w = factor(w, levels = levels(plot_df$w)), 
           v = factor(v, levels = levels(plot_df$v)),         
           label = LETTERS[row_number()]) %>%
    as.data.table()
  
  # Add text to existing figure
  g = g + geom_text(data = label_df, mapping = aes(label = label), 
                    size = 8 ,fontface = "bold", colour = "black", 
                    hjust = "middle", vjust = "top")
  
  # Save the updated figure
  fig_save(o, g, fig_name)
}

# ---------------------------------------------------------
# Convert full scenario name to something concise
# ---------------------------------------------------------
scenario_name_fn = function(x) {
  
  # Store IDs before culling them with str_replace
  id = names(x)
  
  # Reduce scenario names down
  x1 = x %>%
    str_replace_all("Analysis [0-9]+ \\(", "") %>%
    str_replace_all(", Coverage *.+", "") %>%
    str_replace_all(": ", "=") %>%
    str_replace_all(", ", "\n") %>%
    str_replace_all("\\)", "")
  
  # Now need to change Immune escape to percentage...
  
  # Etract Immune evading substring and assocaited value
  x2 = str_extract(x1, "Immuno.*+")
  x3 = str_extract(x2, "[0-9]+.*+")
  
  # Reform string into percentage
  x4 = paste0("Immune evading=", as.numeric(x3) * 100, "%")
  
  # Replace proportion with percentage and rename vector
  y = str_replace_all(x1, x2, x4) %>% setNames(id)
  
  return(y)
}

# ---------------------------------------------------------
# Convert fig name to path of associated plot_info file
# ---------------------------------------------------------
heatmap_id = function(o, fig_name) {
  
  # Remove characters
  file_id = fig_name %>%
    tolower() %>%
    paste(collapse = "_") %>%
    str_replace_all(" ", "_")
  
  # Concatenate with file path and extension
  file_pth = paste0(o$pth$figures, file_id, ".rds")
  
  return(file_pth)
}

