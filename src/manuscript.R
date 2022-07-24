###########################################################
# MANUSCRIPT FIGURES
#
# Results for booster manuscript (June 2022).
#
###########################################################

# ---------------------------------------------------------
# Main plot-calling function
# ---------------------------------------------------------
manuscript_figures = function(o) {
  
  # These results are for 'booster' analysis
  if (o$analysis_name != "booster")
    stop("To produce these results, you must run steps 1-2 for 'booster' analysis")
  
  # ---- Plotting flags ----
  
  # Plot everything
  do_all = FALSE
  
  # Plotting flags: main figures
  do_fig1 = FALSE  # Booster strategies
  do_fig2 = FALSE  # Cumulative strategies (including 6m variants)
  do_fig3 = FALSE  # Booster coverage
  do_fig4 = FALSE  # Timing bars
  do_fig5 = TRUE  # Sensitivity bars
  do_fig6 = FALSE  # Impact by priority group
  
  # Plotting flags: supplement figures
  do_sup1  = FALSE  # Booster strategies (no variant)
  do_sup2  = FALSE  # Booster strategies (1y variant)
  do_sup3  = FALSE  # Booster strategies (6m variant)
  do_sup4  = FALSE  # Cumulative strategies (including 6m variants)
  do_sup5  = FALSE  # Impact of key strategies
  do_sup6  = FALSE  # Timing curves
  do_sup7  = FALSE  # Priority groups
  do_sup8  = FALSE  # Doses by priority group
  do_sup9  = FALSE  # Susceptibility
  do_sup10 = FALSE  # Variant prevalence
  do_sup11 = FALSE  # Variant dominance
  do_sup12 = FALSE  # Number of infections histogram
  do_sup13 = FALSE  # Sensitivity curves
  do_sup14 = FALSE  # Sensitivity timing
  do_sup15 = FALSE  # Sensitivity bars per 100k doses
  do_sup16 = FALSE  # Immunity profile
  do_sup17 = FALSE  # Seasonality profile
  
  # Table flags
  do_tab1 = FALSE  # Cases averted relative to no boosters
  
  # ---- Metrics, scenarios, and colours ----
  
  # Key metrics for plotting
  metrics = list(
    key   = c("all_new_infections", "hospital_beds", "total_doses"), 
    cov   = c("all_new_infections", "hospital_beds", "booster_coverage_12m"),
    cum   = c("all_new_infections", "hospital_admissions", "n_doses"), 
    sens  = c("hospital_admissions"), # c("all_new_infections", "hospital_admissions"),
    time1 = c("all_new_infections", "hospital_admissions"), 
    time2 = c("all_new_infections", "hospital_beds", "R_effective", 
              "total_doses", "pop_susceptibility", "seasonality"))
  
  # Descriptions for each scenario type
  scenarios = 
    c(x0 = "No boosters",
      x1 = "Booster vulnerable every 12 months", 
      x2 = "Booster all eligible every 12 months", 
      x3 = "Booster all eligible every 12 months\n plus vulnerable every 6 months", 
      x4 = "Booster all eligible every 6 months", 
      xA = "Annual boosters", 
      xB = "Biannual boosters", 
      v0 = "No emerging variant", 
      v1 = "Novel variant emerging every 12 months", 
      v2 = "Novel variant emerging every 6 months", 
      c1 = "universal booster coverage", 
      c2 = "lower booster coverage", 
      s1 = "Seasonality effect",
      s2 = "Vaccine infection block",
      s3 = "Variant timing",
      s4 = "Variant severity",
      s5 = "Variant infectivity",
      s6 = "Variant immune evasion",
      m1 = "6 months before", 
      m2 = "5 months before", 
      m3 = "4 months before", 
      m4 = "3 months before", 
      m5 = "2 months before", 
      m6 = "1 month before", 
      m7 = "Peak winter", 
      m8 = "1 month after", 
      m9 = "2 months after", 
      t1 = "Seasonality effect", 
      t2 = "Number of doses per day", 
      b  = "Parameter best estimates",
      u  = "Parameter upper bound", 
      l  = "Parameter lower bound")
  
  # Define colour schemes
  colours = list(
    key    = "brewer::set1",
    cov    = "darkblue",
    sens   = "viridis::inferno", 
    timing = "viridis::viridis")
  
  # ---- Figure properties ----
  
  # Set up scenarios, seasonality, and variants for all figures
  f = fig_setup(o, metrics, scenarios, colours)
  
  # Seasonality colourbar name and palette
  f$season_name = "Season"
  f$season_cols = "RdBu"
  
  # Stretch the colourbar limits for more 'mid' colour
  f$season_stretch = 0.04
  
  # Figure fontsize (title, facets, legend, ticks)
  f$fontsize = c(9, 8, 7, 7) 
  
  # Wrap after every n chars
  f$text_wrap = c(24, 20)
  
  # Nudge variant emergence date text
  f$nudge = list(x = 10, y = 2)
  
  # ---- Figure 1 ----
  
  # Check plotting flag
  if (do_fig1 || do_all) {
    
    message(" - Manuscript figure 1")
    
    # Booster strategies
    fig_name = "Fig 1 - Booster strategies"
    plot_temporal(o, fig_name   = NULL,
                  plot_baseline = FALSE,
                  scenarios     = f$scenarios$key[, 1:2],
                  plot_metrics  = f$metrics$key,
                  legend_rows   = 1,
                  line_width    = 0.75, 
                  n_wrap        = f$text_wrap, 
                  override_colours = f$colours$key, 
                  override_fontsize = f$fontsize) %>%
      plot_variant_emergence(f, generic_names = FALSE) %>%
      plot_season_colourbar(f, save = fig_name)
  }
  
  # ---- Figure 2 ----
  
  # Check plotting flag
  if (do_fig2 || do_all) {
    
    message(" - Manuscript figure 2")
    
    # Define a custom faceting function list
    facet_custom = list(fn = "ggh4x::facet_grid2", 
                        cols   = "vars(scenario_group)",
                        rows   = "vars(metric)", 
                        scales = "'free'", 
                        independent = "'y'", 
                        labeller    = "f$label_grid")
    
    # Custom labels (needed as ggh4x messes up labelling)
    label_custom = qc(A, D, G, B, E, H, C, F, I)
    
    # Cumulative strategies, including 6m variants
    fig_name = "Fig 2 - Cumulative strategies"
    plot_temporal(o, fig_name   = NULL,
                  plot_baseline = FALSE,
                  scenarios     = t(f$scenarios$key),
                  plot_metrics  = f$metrics$cum,
                  cumulative    = TRUE,
                  aes_reverse   = TRUE,
                  aes_linetype  = FALSE,
                  line_width    = 0.75, 
                  n_wrap        = f$text_wrap, 
                  facet_custom  = parse_fn(facet_custom, evaluate = FALSE),
                  facet_labels  = FALSE,
                  override_colours = f$colours$key, 
                  override_fontsize = f$fontsize) %>%
      facet_labels(tag_pool = label_custom) %>%
      plot_variant_emergence(f) %>%
      plot_season_colourbar(f, facet = "grid2", save = fig_name)
  }
  
  # ---- Figure 3 ----
  
  # Check plotting flag
  if (do_fig3 || do_all) {
    
    message(" - Manuscript figure 3")
    
    # Ignore the first handful of days
    date_shift = 60
    
    # Define a custom faceting function list
    facet_custom = list(fn = "ggh4x::facet_grid2", 
                        cols   = "vars(scenario_group)",
                        rows   = "vars(metric)", 
                        scales = "'free'", 
                        independent = "'y'", 
                        labeller    = "f$label_grid")
    
    # Custom labels (needed as ggh4x messes up labelling)
    label_custom = qc(A, D, B, E, C, F)
    
    # Booster coverage
    fig_name = "Fig 3 - Booster coverage"
    plot_temporal(o, fig_name   = NULL,
                  plot_baseline = FALSE,
                  scenarios     = f$scenarios$cov[1:2, ],
                  plot_metrics  = f$metrics$key,
                  aes_reverse   = TRUE,
                  aes_linetype  = FALSE,
                  plot_from     = date_shift,
                  legend_rows   = 1,
                  line_width    = 0.75, 
                  n_wrap        = f$text_wrap[1], 
                  facet_custom  = parse_fn(facet_custom, evaluate = FALSE),
                  facet_labels  = FALSE,
                  override_colours = f$colours$cov, 
                  override_fontsize = f$fontsize) %>%
      facet_labels(tag_pool = label_custom) %>%
      plot_variant_emergence(f) %>%
      plot_season_colourbar(f, save = fig_name, 
                            facet   = "grid2", 
                            shift   = date_shift)
  }
  
  # ---- Figure 4 ----
  
  # Check plotting flag
  if (do_fig4 || do_all) {
    
    message(" - Manuscript figure 4") 
    
    # Timing bars
    fig_name = "Fig 4 - Timing bars"
    plot_impact(o, fig_name, 
                plot_baseline = FALSE, 
                scenarios     = f$scenarios$timing, 
                plot_metrics  = f$metrics$time1, 
                person_days   = f$person_days,
                n_wrap        = 40, 
                override_colours = f$colours$timing, 
                override_fontsize = f$fontsize)
  }
  
  # ---- Figure 5 ----
  
  # Check plotting flag
  if (do_fig5 || do_all) {
    
    message(" - Manuscript figure 5")
    
    # Sensitivity bars
    fig_name = "Fig 5 - Sensitivity bars"
    plot_sensitivity(f, fig_name)
  }
  
  # ---- Figure 6 ----
  
  # Check plotting flag
  if (do_fig6 || do_all) {
    
    message(" - Manuscript figure 6")
    
    # Impact by priority group
    fig_name = "Fig 6 - Impact by priority group"
    g = plot_impact(o, fig_name   = NULL,
                    plot_baseline = FALSE,
                    scenarios     = c("x0v1c1", "x1v1c1", "x2v1c1"),
                    plot_metrics  = f$metrics$cum[-1], 
                    plot_by       = "priority_group",
                    person_days   = f$person_days,
                    group_stack   = FALSE,
                    legend_rows   = 1,
                    n_wrap        = c(24, 32),
                    x_wrap        = 21,
                    x_rotate      = 0,
                    y_expand      = 0.1,
                    override_fontsize = f$fontsize)
    
    # Details of dose data plotted
    doses_df = ggplot_build(g)$plot$data %>%
      filter(grepl("vaccine", metric))
    
    # Create text label for empty 'no booster' facet
    text_df = doses_df %>%
      filter(grepl("No boosters", scenario)) %>%
      slice_head(n = 1) %>%
      mutate(label = "No vaccine booster doses\n administered")
    
    # Place text label in center of facet
    g = g + geom_text(data    = text_df, 
                      mapping = aes(label = label), 
                      x       = (length(levels(doses_df$group)) + 1) / 2, 
                      y       = max(doses_df$upper) / 2, 
                      vjust   = 0.5, 
                      hjust   = 0.5, 
                      size    = 3)
    
    # Save figure
    fig_save(o, g, fig_name)
  }
  
  # ---- Supplementary figure 1 ----
  
  # Check plotting flag
  if (do_sup1 || do_all) {
    
    message(" - Supplement figure 1")
    
    # Booster strategies: no variant
    fig_name = "Fig S01a - Booster strategies (no variant)"
    plot_temporal(o, fig_name   = NULL,
                  alt_baseline  = "x0v0c1",
                  scenarios     = f$scenarios$key[-1, 1, drop = FALSE],
                  plot_metrics  = f$metrics$key[1 : 2],
                  legend_rows   = 1,
                  line_width    = 1.5, 
                  n_wrap        = f$text_wrap, 
                  override_colours = f$colours$key, 
                  override_fontsize = f$fontsize) %>%
      plot_season_colourbar(f, save = fig_name)
    
    # Booster strategies: no variant
    fig_name = "Fig S01b - Booster strategies (no variant)"
    plot_temporal(o, fig_name   = NULL,
                  plot_baseline = FALSE,
                  scenarios     = f$scenarios$key[, 1],
                  plot_metrics  = f$metrics$key[1 : 2],
                  legend_rows   = 5,
                  line_width    = 1.5, 
                  n_wrap        = f$text_wrap, 
                  override_colours = f$colours$key, 
                  override_fontsize = f$fontsize) %>%
      plot_season_colourbar(f, facet = "wrap", save = fig_name)
  }
  
  # ---- Supplementary figure 2 ----
  
  # Check plotting flag
  if (do_sup2 || do_all) {
    
    message(" - Supplement figure 2")
    
    # Booster strategies: with variant
    fig_name = "Fig S02a - Booster strategies (1y variant)"
    plot_temporal(o, fig_name   = NULL,
                  alt_baseline  = "x0v1c1",
                  scenarios     = f$scenarios$key[-1, 2, drop = FALSE],
                  plot_metrics  = f$metrics$key[1 : 2],
                  legend_rows   = 1,
                  line_width    = 1.5, 
                  n_wrap        = f$text_wrap, 
                  override_colours = f$colours$key, 
                  override_fontsize = f$fontsize) %>%
      plot_variant_emergence(f) %>%
      plot_season_colourbar(f, save = fig_name)
    
    # Booster strategies: with variant
    fig_name = "Fig S02b - Booster strategies (1y variant)"
    plot_temporal(o, fig_name   = NULL,
                  plot_baseline = FALSE,
                  scenarios     = f$scenarios$key[, 2],
                  plot_metrics  = f$metrics$key[1 : 2],
                  legend_rows   = 5,
                  line_width    = 1.5, 
                  n_wrap        = f$text_wrap, 
                  override_colours = f$colours$key, 
                  override_fontsize = f$fontsize) %>%
      plot_variant_emergence(f) %>%
      plot_season_colourbar(f, facet = "wrap", save = fig_name)
  }
  
  # ---- Supplementary figure 3 ----
  
  # Check plotting flag
  if (do_sup3 || do_all) {
    
    message(" - Supplement figure 3")
    
    # Booster strategies: with variant
    fig_name = "Fig S03a - Booster strategies (6m variant)"
    plot_temporal(o, fig_name   = NULL,
                  alt_baseline  = "x0v2c1",
                  scenarios     = f$scenarios$key[-1, 3, drop = FALSE],
                  plot_metrics  = f$metrics$key[1 : 2],
                  legend_rows   = 1,
                  line_width    = 1.5, 
                  n_wrap        = f$text_wrap, 
                  override_colours = f$colours$key, 
                  override_fontsize = f$fontsize) %>%
      plot_variant_emergence(f) %>%
      plot_season_colourbar(f, save = fig_name)
    
    # Booster strategies: with variant
    fig_name = "Fig S03b - Booster strategies (6m variant)"
    plot_temporal(o, fig_name   = NULL,
                  plot_baseline = FALSE,
                  scenarios     = f$scenarios$key[, 3],
                  plot_metrics  = f$metrics$key[1 : 2],
                  legend_rows   = 5,
                  line_width    = 1.5, 
                  n_wrap        = f$text_wrap, 
                  override_colours = f$colours$key, 
                  override_fontsize = f$fontsize) %>%
      plot_variant_emergence(f) %>%
      plot_season_colourbar(f, facet = "wrap", save = fig_name)
  }
  
  # ---- Supplementary figure 4 ----
  
  # Check plotting flag
  if (do_sup4 || do_all) {
    
    message(" - Supplement figure 4")
    
    # Cumulative booster strategies: grouped by policy
    fig_name = "Fig S04 - Cumulative booster strategies"
    plot_temporal(o, fig_name   = NULL,
                  plot_baseline = FALSE,
                  scenarios     = f$scenarios$key,
                  plot_metrics  = f$metrics$cum,
                  cumulative    = TRUE,
                  legend_rows   = 1,
                  line_width    = 1.5, 
                  n_wrap        = f$text_wrap, 
                  override_colours = f$colours$key, 
                  override_fontsize = f$fontsize) %>%
      plot_variant_emergence(f, generic_names = TRUE) %>%
      plot_season_colourbar(f, save = fig_name)
  }
  
  # ---- Supplementary figure 5 ----
  
  # Check plotting flag
  if (do_sup5 || do_all) {
    
    message(" - Supplement figure 5")
    
    # Impact of key strategies
    fig_name = "Fig S05 - Impact of key strategies"
    plot_impact(o, fig_name,
                plot_baseline = FALSE,
                scenarios     = t(f$scenarios$all[1:2, ]),
                plot_metrics  = f$metrics$cum,
                person_days   = f$person_days,
                bar_legend    = TRUE,
                legend_rows   = 3,
                n_wrap        = f$text_wrap,
                override_colours = f$colours$all,
                override_fontsize = f$fontsize)
  }
  
  # ---- Supplementary figure 6 ----
  
  # Check plotting flag
  if (do_sup6 || do_all) {
    
    message(" - Supplement figure 6")
    
    # Timing curves
    fig_name = "Fig S06 - Timing curves"
    plot_temporal(o, fig_name = NULL,
                  plot_baseline = FALSE,
                  scenarios     = f$scenarios$timing,
                  plot_metrics  = f$metrics$time2,
                  n_wrap        = f$text_wrap,
                  override_colours = f$colours$timing,
                  override_fontsize = f$fontsize) %>%
      plot_season_colourbar(f, save = fig_name,
                            facet   = "wrap",
                            season  = 2)
  }
  
  # ---- Supplementary figure 7 ----
  
  # Check plotting flag
  if (do_sup7 || do_all) {
    
    message(" - Supplement figure 7")
    
    # Priority groups
    fig_name = "Fig S07 - Priority groups"
    plot_impact(o, fig_name,
                plot_baseline = FALSE,
                scenarios     = f$scenarios$variant_no,
                plot_metrics  = f$metrics$cum, 
                plot_by       = "priority_group",
                person_days   = f$person_days,
                group_stack   = FALSE,
                legend_rows   = 1,
                n_wrap        = c(21, 22), # f$text_wrap, 
                scenario_name_fn = simple_scenarios, 
                override_fontsize = f$fontsize)
  }
  
  # ---- Supplementary figure 8 ----
  
  # Check plotting flag
  if (do_sup8 || do_all) {
    
    message(" - Supplement figure 8")
    
    # Doses by priority group
    fig_name = "Fig S08 - Doses by priority group"
    plot_temporal(o, fig_name   = NULL,
                  scenarios     = f$scenarios$variant_6m[-1],
                  plot_baseline = FALSE,
                  plot_metrics  = "n_doses",
                  cumulative    = TRUE,
                  plot_by       = "priority_group", 
                  plot_geom     = "area", 
                  n_wrap        = 40, # f$text_wrap,
                  scenario_name_fn = partial_scenarios, 
                  override_fontsize = f$fontsize) %>%
      plot_variant_emergence(f) %>%
      plot_season_colourbar(f, save = fig_name, 
                            facet   = "wrap")
  }
  
  # ---- Supplementary figure 9 ----
  
  # Check plotting flag
  if (do_sup9 || do_all) {
    
    message(" - Supplement figure 9")
    
    # Specialist metricd for this supplementary figure
    sus_metrics = c("total_doses", "all_new_infections", "pop_susceptibility")
    
    # Define a custom faceting function list
    facet_custom = list(fn = "ggh4x::facet_grid2", 
                        cols   = "vars(scenario_group)",
                        rows   = "vars(metric)", 
                        scales = "'free'", 
                        independent = "'y'", 
                        labeller    = "f$label_grid")
    
    # Custom labels (needed as ggh4x messes up labelling)
    label_custom = qc(A, D, G, B, E, H, C, F, I)
    
    # Susceptibility
    fig_name = "Fig S09 - Susceptibility"
    plot_temporal(o, fig_name = NULL,
                  plot_baseline = FALSE,
                  scenarios     = f$scenarios$all,
                  plot_metrics  = sus_metrics,
                  aes_reverse   = TRUE,
                  aes_linetype  = FALSE,
                  legend_rows   = 3,
                  line_width    = 1.5,
                  n_wrap        = f$text_wrap,
                  facet_custom  = parse_fn(facet_custom, evaluate = FALSE),
                  facet_labels  = FALSE,
                  override_colours = f$colours$all,
                  override_fontsize = f$fontsize) %>%
      facet_labels(tag_pool = label_custom) %>%
      plot_variant_emergence(f) %>%
      plot_season_colourbar(f, facet = "grid2", save = fig_name)
  }
  
  # ---- Supplementary figure 10 ----
  
  # Check plotting flag
  if (do_sup10 || do_all) {
    
    message(" - Supplement figure 10")
    
    # Variant prevalence
    fig_name = "Fig S10a - Variant prevalence"
    plot_temporal(o, fig_name   = NULL,
                  plot_baseline = FALSE,
                  scenarios     = f$scenarios$variant_1y,
                  plot_geom     = "area",
                  plot_by       = "variant",
                  n_wrap        = 20, 
                  n_wrap        = f$text_wrap, 
                  scenario_name_fn = partial_scenarios, 
                  override_fontsize = f$fontsize) %>%
      plot_variant_emergence(f) %>%
      plot_season_colourbar(f, save = fig_name)
    
    # Variant prevalence
    fig_name = "Fig S10b - Variant prevalence"
    plot_temporal(o, fig_name   = NULL,
                  plot_baseline = FALSE,
                  scenarios     = f$scenarios$variant_6m,
                  plot_geom     = "area",
                  plot_by       = "variant",
                  n_wrap        = 20, 
                  n_wrap        = f$text_wrap, 
                  scenario_name_fn = partial_scenarios, 
                  override_fontsize = f$fontsize) %>%
      plot_variant_emergence(f) %>%
      plot_season_colourbar(f, save = fig_name)
  }
  
  # ---- Supplementary figure 11 ----
  
  # Check plotting flag
  if (do_sup11 || do_all) {
    
    message(" - Supplement figure 11")
    
    # Variant dominance
    fig_name = "Fig S11a - Variant dominance"
    plot_impact(o, fig_name,
                plot_baseline = FALSE,
                scenarios     = f$scenarios$variant_1y,
                plot_metrics  = f$metrics$cum,
                plot_by       = "variant", 
                person_days   = f$person_days,
                group_stack   = FALSE, 
                bar_width     = 1, 
                legend_rows   = 1,
                n_wrap        = c(23, 20), # f$text_wrap, 
                scenario_name_fn = partial_scenarios, 
                override_fontsize = f$fontsize)
    
    # Variant dominance
    fig_name = "Fig S11b - Variant dominance"
    plot_impact(o, fig_name,
                plot_baseline = FALSE,
                scenarios     = f$scenarios$variant_6m,
                plot_metrics  = f$metrics$cum,
                plot_by       = "variant", 
                person_days   = f$person_days,
                group_stack   = FALSE, 
                bar_width     = 1, 
                legend_rows   = 1,
                n_wrap        = f$text_wrap, 
                scenario_name_fn = partial_scenarios, 
                override_fontsize = f$fontsize)
  }
  
  # ---- Supplementary figure 12 ----
  
  # Check plotting flag
  if (do_sup12 || do_all) {
    
    message(" - Supplement figure 12")
    
    # Number of infections histogram: grouped by scenario
    fig_name = "Fig S12a - Number of infections"
    plot_num_infections(o, fig_name, 
                        plot_baseline = FALSE,
                        scenarios     = f$scenarios$variant_1y,
                        facet_by      = "scenario", 
                        n_wrap        = 38, 
                        scenario_name_fn = partial_scenarios, 
                        override_colours = f$colours$key, 
                        override_fontsize = f$fontsize)
    
    # Number of infections histogram: grouped by num infections
    fig_name = "Fig S12b - Number of infections"
    plot_num_infections(o, fig_name, 
                        plot_baseline = FALSE,
                        scenarios     = f$scenarios$variant_1y,
                        facet_by      = "group", 
                        legend_rows   = 5, 
                        scenario_name_fn = partial_scenarios, 
                        override_colours = f$colours$key, 
                        override_fontsize = f$fontsize)
  }
  
  # ---- Supplementary figure 13 ----
  
  # Check plotting flag
  if (do_sup13 || do_all) {
    
    message(" - Supplement figure 13")
    
    # Define a custom faceting function list
    facet_custom = list(fn = "facet_grid", 
                        cols   = "vars(scenario_group)",
                        rows   = "vars(scenario_type)", 
                        labeller = "f$label_grid")
    
    # Custom labels (needed as have more than 26 facets)
    label_custom = c(LETTERS, paste0("A", LETTERS))
    
    # Sensitivity analysis: curves
    fig_name = "Fig S13 - Sensitivity curves"
    plot_temporal(o, fig_name,
                  plot_baseline = FALSE,
                  scenarios     = f$scenarios$sens_mat,
                  aes_linetype  = FALSE,
                  plot_metrics  = f$metrics$sens,
                  legend_rows   = 1,
                  line_width    = 0.5, 
                  n_wrap        = c(20, 15),
                  facet_custom  = parse_fn(facet_custom, evaluate = FALSE),
                  facet_labels  = FALSE,
                  override_colours = f$colours$sens_mat,
                  override_fontsize = f$fontsize) %>% 
      apply_theme(panel.grid.major.y = element_line(), 
                  panel.grid.minor.y = element_line()) %>%
      facet_labels(tag_pool = label_custom) %>%
      plot_variant_emergence(f) %>%
      plot_season_colourbar(f, save = fig_name,
                            height  = 0.12,
                            legend  = FALSE)
  }
  
  # ---- Supplementary figure 14 ----
  
  # Check plotting flag
  if (do_sup14 || do_all) {
    
    message(" - Supplement figure 14")
    
    # Timing sensitivity
    fig_name = "Fig S14a - Timing sensitivity"
    plot_impact(o, fig_name,
                plot_baseline = FALSE,
                scenarios     = f$scenarios$timing_mat,
                plot_metrics  = f$metrics$time1,
                person_days   = f$person_days,
                bar_legend    = TRUE,
                legend_rows   = 3,
                n_wrap        = f$text_wrap,
                override_colours = f$colours$timing,
                override_fontsize = c(18, 14, 14, 14))
    
    # Timing sensitivity
    fig_name = "Fig S14b - Timing sensitivity"
    plot_timing_sensitivity(f, fig_name)
  }
  
  # ---- Supplementary figure 15 ----
  
  # Check plotting flag
  if (do_sup15 || do_all) {
    
    message(" - Supplement figure 15")
    
    # Sensitivity bars
    fig_name = "Fig S15 - Per dose"
    plot_sensitivity(f, fig_name,
                     plot_metric  = "avert_dose",
                     group_policy = FALSE)
  }
  
  # ---- Supplementary figure 16 ----
  
  # Check plotting flag
  if (do_sup16 || do_all) {
    
    message(" - Supplement figure 16")
    
    # Immunity profile: acquired and vaccine-induced
    fig_name = "Fig S16 - Immunity profile"
    plot_immunity_profiles(o, fig_name, f$model_input)
  }
  
  # ---- Supplementary figure 17 ----
  
  # Check plotting flag
  if (do_sup17 || do_all) {
    
    message(" - Supplement figure 17")
    
    # Load scenarios with seasonal extremes
    lower = try_load(o$pth$scenarios, "s1l_x0")$input
    upper = try_load(o$pth$scenarios, "s1u_x0")$input
    
    # Define colours and descriptive labels
    colours = c("gold", o$baseline_colour, "forestgreen")
    labels  = c("Lower bound (scaler = 0.2)", 
                "Best estimate (scaler = 0.3)", 
                "Upper bound (scaler = 0.4)")
    
    # Immunity profile: acquired and vaccine-induced
    fig_name = "Fig S17 - Seasonality profile"
    plot_seasonality_profile(o, NULL, 
                             lower, f$model_input, upper, 
                             colours = colours, 
                             labels  = labels) %>%
      plot_season_colourbar(f, save = fig_name,
                            facet   = "wrap", 
                            height  = 0.2)
    
  }
  
  # ---- Cases averted table ----
  
  # Check flag
  if (do_tab1 || do_all) {
    
    message(" - Cases averted table")
    
    # Preallocate matrix 
    data_mat = t(f$scenarios$all[, -1])
    data_mat[] = NA
    
    # Remove newline chars from scenario descriptions
    rownames(data_mat) = str_remove(rownames(data_mat), "\n")
    
    # Define two additonal matrices with same structure
    lower_mat = upper_mat = data_mat
    
    # Loop through variant assumptions
    for (i in 1 : nrow(f$scenarios$all)) {
      
      # Scenario IDs for this variant assumption
      baseline  = f$scenarios$all[i, 1]
      scenarios = f$scenarios$all[i, -1]
      
      # Load ouput for the baseline comparator
      base_result = try_load(o$pth$scenarios, paste0(baseline, "_raw"))
      
      # Extract cumulative value(s) for no booster baseline
      base_value = base_result %>%
        filter(metric %in% f$metrics$sens, 
               is.na(group)) %>%
        group_by(scenario, metric, seed) %>%
        summarise(total = sum(value)) %>%
        pull(total)
      
      # Loop through each alternative scenario
      for (j in 1 : length(scenarios)) {
        scen = scenarios[j]
        
        # Load ouput for this scenario
        scen_output = try_load(o$pth$scenarios, paste0(scen, "_raw"))
        
        # Again, extract cumulative value for each metric
        scen_value = scen_output %>%
          filter(metric %in% f$metrics$sens, 
                 is.na(group)) %>%
          group_by(scenario, metric, seed) %>%
          summarise(total = sum(value)) %>%
          pull(total)
        
        # Take the proportioinal difference
        averted = (base_value - scen_value) / base_value
        
        # Mean value across all seeds
        data_mat[j, i] = mean(averted)
        
        # Min and max values across all seeds
        lower_mat[j, i] = quantile(averted, probs = 0.025)
        upper_mat[j, i] = quantile(averted, probs = 1 - 0.025)
      }
    }
    
    # Save the populated matrix to file
    file_name = paste0(o$pth$figures, "Cases averted - ")
    write.table(data_mat, file = paste0(file_name, "Mean.csv"), sep = ";")
    
    # Also save bounds (as seperate files)
    write.table(lower_mat, file = paste0(file_name, "Lower.csv"), sep = ";")
    write.table(upper_mat, file = paste0(file_name, "Upper.csv"), sep = ";")
  }
}

# ---------------------------------------------------------
# Set up scenarios, seasonality, and variants for all figures
# ---------------------------------------------------------
fig_setup = function(o, metrics, scenarios, colours) {
  
  # Initiate figure properties and copy list of metric details
  f = list(metrics = metrics)
  
  # ---- Scenarios ----
  
  # Initiate scenario list and copy descriptions
  s = list(name = scenarios)
  
  # IDs of all scenarios defined in this analysis
  s_id = names(parse_yaml(o, "*read*", read_array = TRUE)[-1])
  
  # Differentiate scenarios with and without emerging variants
  s$variant_no = s_id[grepl("v0c1", s_id)]
  s$variant_1y = s_id[grepl("v1c1", s_id)]
  s$variant_6m = s_id[grepl("v2c1", s_id)]
  
  # Create a matrix of scenario IDs so we can plot by scenario group
  s$key = t(rbind(s$variant_no, s$variant_1y, s$variant_6m))
  
  # Define grouping and element names (rows and columns respectively)
  rownames(s$key) = s$name[grepl("x[0-9]+", names(s$name))]
  colnames(s$key) = s$name[grepl("v[0-9]+", names(s$name))]
  
  # Differentiate scenarios by booster coverage
  s$cov = rbind(s_id[grepl("x2v0", s_id)], 
                s_id[grepl("x2v1", s_id)],
                s_id[grepl("x2v2", s_id)])
  
  # Define grouping and element names (rows and columns respectively)
  rownames(s$cov) = s$name[grepl("v[0-9]+", names(s$name))]
  colnames(s$cov) = s$name[grepl("c[0-9]+", names(s$name))]
  
  # Initiate a matrix that'll combine 'key' and 'cov'
  s_all = t(s$cov)
  
  # Rename the coverage scenarios so they can be differentiated
  rownames(s_all) = paste0(s$name[["x2"]], " (", rownames(s_all), ")")
  
  # Bind with s$key and remove duplicate
  s$all = unique(rbind(s_all, s$key))
  
  # Retain order and transpose
  s$all = t(s$all[order(s$all[, 1]), ])
  
  # Vector of all sensitivity analysis figures
  s$sens_vec = s_id[grepl("s[0-9]+", s_id)]
  
  # Indices for parameter & policy descriptions
  row_idx = grepl("s[0-9]+", names(s$name))
  col_idx = grepl("x[0-9]+", names(s$name))
  
  # Reshape vector into matrix: param x policy
  s$sens_mat = t(matrix(s$sens_vec, nrow = sum(col_idx)))
  
  # Set row and column names
  rownames(s$sens_mat) = s$name[row_idx] %>% rep(each = 2)
  colnames(s$sens_mat) = s$name[col_idx]
  
  # Timing scenarios
  s$timing   = s_id[grepl("m[0-9]+_[a-z]+", s_id)]
  all_timing = s_id[grepl("m[0-9]+", s_id)]
  
  # Reshape timing sensitivity scenarios into a matrix
  s$timing_mat = matrix(all_timing, nrow = length(s$timing))
  
  timing_names = s$name[grepl("m[0-9]+", names(s$name))]
  
  param_names = s$name[grepl("t[0-9]+", names(s$name))]
  bound_names = s$name[grepl("[u,l]", names(s$name))]
  
  param_bounds = param_names %>% 
    expand_grid(bound_names) %>%
    unite("x", sep = "\n ") %>%
    pull(x)
  
  rownames(s$timing_mat) = timing_names
  colnames(s$timing_mat) = c(s$name[["b"]], param_bounds)
  
  # Append values to f list
  f$scenarios = s
  
  # ---- Colours ----
  
  # Colours of key scenarios (as defined by palette_scenario)
  colours_key   = colour_scheme(colours$key, n = nrow(s$key) - 1)
  f$colours$key = c(o$baseline_colour, colours_key)
  
  # Index colour of scenario we explore lower coverage for
  col_idx = which(grepl("x2", s$variant_no))
  col_cov = c(f$colours$key[col_idx], colours$cov)
  
  # Combine key and cov scenarios
  col_all = c(f$colours$key, col_cov)
  
  # Temporarily convert to named vector to set appropriate ordering
  col_all = setNames(col_all, c(s$variant_no, s$cov[1, ]))
  col_all = unique(col_all[order(names(col_all))])
  
  # Append these colours to f list
  f$colours$cov = col_cov
  f$colours$all = col_all
  
  # Set unique colour scheme for timing scenarios
  f$colours$timing = colour_scheme(colours$timing, n = length(s$timing))
  
  # Colours represent strategies for main sensitivity plot
  f$colours$sens_vec = f$colours$key[-1]
  
  # Temporal sensitivity plot uses unique colour scheme
  col_sens = colour_scheme(colours$sens, n = nrow(s$sens_mat) / 2 + 2)
  
  # Remove the two most extreme colours (too dark/light respectively)
  col_sens = rev(rev(col_sens[-1])[-1])
  
  # Repeat for lower and upper bounds
  f$colours$sens_mat = rep(col_sens, each = 2)
  
  # ---- Seasonality ----
  
  # Load a scenario that models emerging variants
  f$model_input = try_load(o$pth$scenarios, "baseline")$input
  timing_input  = try_load(o$pth$scenarios, s$timing[1])$input
  
  # Seasonality curve 1: All key future scenarios
  f$season1 = data.table(date  = 1 : f$model_input$n_days, 
                         value = f$model_input$seasonality)
  
  # Seasonality curve 2: Timing scenarios
  f$season2 = data.table(date  = 1 : timing_input$n_days, 
                         value = timing_input$seasonality)
  
  # Use model input to calculate 'per 100k people per year' scaler
  f$person_days = 1e5 / (f$model_input$n_days / 365)
  
  # ---- Variants ----
  
  # Load scenarios that model emerging variants
  variant_1y = try_load(o$pth$scenarios, s$variant_1y[1])$input
  variant_6m = try_load(o$pth$scenarios, s$variant_6m[1])$input
  
  # We want to know the date variants were imported
  f$variant = list(
    name     = variant_1y$variants$name[-1], 
    annual   = variant_1y$variants$import_day[-1], 
    biannual = variant_6m$variants$import_day[-1])
  
  return(f)
}

# ---------------------------------------------------------
# Simplify scenario names by removing all details within brackets
# ---------------------------------------------------------
simple_scenarios = function(x) {
  
  # Remove everything incased within ( and )
  y = str_remove(x, " \\(.+\\)") %>%
    setNames(names(x))
  
  return(y)
}

# ---------------------------------------------------------
# Simplify scenario names by removing some details within brackets
# ---------------------------------------------------------
partial_scenarios = function(x) {
  
  # Remove first set of text within brackets
  y = str_replace(x, " \\(.+\\, ", " (") %>%
    setNames(names(x))
  
  return(y)
}

# ---------------------------------------------------------
# Apply key-value paairs to a ggplot object theme
# ---------------------------------------------------------
apply_theme = function(g, ...) {
  
  # Pass key-value pairs directly into theme
  g = g + theme(...)
  
  return(g)
}

# ---------------------------------------------------------
# Add arrows to plots to show when variants emerge
# ---------------------------------------------------------
plot_variant_emergence = function(g, f, generic_names = FALSE, save = NULL) {
  
  # Extract original plotting dataframe from ggplot object
  plot_df = ggplot_build(g)$plot$data
  
  # All variant scenarios considered in this plot
  variant_scenarios = plot_df$scenario %>%
    str_extract("[a-z]+ variants") %>%
    unique() %>%
    setdiff("no variants")
  
  # Use specific variant names by default
  variant_names = f$variant$name
  
  # Also able to use generic names
  if (generic_names == TRUE) {
    
    # Split into a list of vectors of single chars
    name_chars = strsplit(f$variant$name, NULL)
    
    # Extract chars that are common to all names
    common_chars = Reduce(function(x, y) y[match(x, y, 0)], name_chars)
    
    # Collapse into string and repeat this name for each variant
    generic_name  = paste(common_chars, collapse = "")
    variant_names = rep(generic_name, length(f$variant$name))
  }
  
  # Loop through the different variant scenarios
  for (variant_scenario in variant_scenarios) {
    
    # Import dates for all variants in this variant scenario
    variant_import = f$variant[[word(variant_scenario, 1)]]
    
    # Loop through each variant
    for (i in seq_along(variant_import)) {
      
      # Ignore trivial variants (ie those not imported)
      if (variant_import[i] > 0) {
        
        # Regular expression needs to match the whole phrase
        reg_exp = paste0("\\<", variant_scenario, "\\>")
        
        # Reduce dataframe down to this date for each variant facet
        variant_df = plot_df %>%
          filter(grepl(reg_exp, scenario)) %>%
          select(-lower, -upper) %>%
          mutate(date  = variant_import[i],
                 value = Inf) %>%
          unique()
        
        # Plot dashed lines at the point variants come in
        g = g + geom_vline(data     = variant_df, 
                           mapping  = aes(xintercept = date), 
                           color    = "grey", 
                           linetype = "dashed", 
                           size     = 0.25)
        
        # Add some white space to the RHS to nudge text into a decent position
        nudge_name   = c(variant_names[i], rep(" ", f$nudge$y))
        variant_name = paste0(nudge_name, collapse = "")
        
        # Add novel variant text next to the dashed line
        g = g + geom_text(data  = variant_df, 
                          label = variant_name, 
                          color = "grey",
                          size  = 2,
                          angle = 90,
                          vjust = 1,
                          hjust = 1,
                          nudge_x = f$nudge$x)
      }
    }
  }
  
  # Save updated figure if desired
  if (!is.null(save))
    fig_save(o, g, save)
  
  return(g)
}

# ---------------------------------------------------------
# Add arrows to plots to show when variants emerge
# ---------------------------------------------------------
plot_season_colourbar = function(g1, f, season = 1, shift = 0, height = 0.15, 
                                 facet = "grid", legend = TRUE, save = NULL) {
  
  # Number of facet columns in original plot
  facet_cols = facet_dims(g1)[2]
  
  # Seasonality profile to plot
  season_df = f[[paste0("season", season)]]
  
  # Repeat seasonality details for each so we can produce a consistent facet_grid
  plot_df = season_df %>%
    filter(date >= shift) %>%
    expand_grid(cols = 1 : facet_cols) %>%
    mutate(metric = f$season_name, y = 0) %>%
    as.data.table()
  
  # Plot seasonality values as tiles, scenario_group times
  g2 = ggplot(plot_df, aes(x = date, y = y, fill = value)) + 
    geom_tile(show.legend = FALSE)
  
  # Consider facetting only if necessary
  if (facet != "none") {
    
    # Construct faceting function
    facet_fn = paste0("facet_", facet)
    
    # Apply faceting function
    if (facet == "wrap") g2 = g2 + get(facet_fn)(~cols)
    if (facet != "wrap") g2 = g2 + get(facet_fn)(metric~cols)
  }
  
  # In order to stretch the colour map, we want extreme values of the data
  min_season = min(plot_df$value)
  max_season = max(plot_df$value)
  
  # Use these extremes to find the mid point
  mid_season = min_season + (max_season - min_season) / 2
  
  # Use colourbar limits to stretch out the colourbar
  lims_season = c(min_season - mid_season * f$season_stretch, 
                  max_season + mid_season * f$season_stretch)
  
  # Set desired colour map and apply the stretching
  g2 = g2 + scale_fill_distiller(palette   = f$season_cols, 
                                 direction = 1, 
                                 limits    = lims_season)
  
  # Remove spacing around the facet
  g2 = g2 + scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0))
  
  # Turn off everything so only left with the tiles
  g2 = g2 + theme_classic() +
    theme(axis.line    = element_blank(),
          axis.title   = element_blank(),
          axis.text    = element_blank(),
          axis.ticks   = element_blank(),
          strip.text.x = element_blank(),
          strip.text   = element_text(size = f$fontsize[2]),
          strip.background = element_blank(), 
          plot.margin  = margin(b = 2, t = 2),
          panel.border = element_rect(size = 0.5, colour = "black", fill = NA))
  
  # Turn off ticks on the main figure we're adding to
  g1 = g1 + theme(plot.margin  = margin(2, 2, 2, 2), 
                  axis.text.x  = element_blank(),
                  axis.ticks.x = element_blank())
  
  # Extract legend from initial plot (if desired)
  if (legend == TRUE)  use_legend = get_legend(g1)
  if (legend == FALSE) use_legend = NULL
  
  # Combine the seasonality colourbars to the main figure
  g = ggarrange(g1, g2,
                nrow    = 2,
                ncol    = 1,
                heights = c(1, height),
                align   = "v",
                legend  = "bottom",
                legend.grob = use_legend) 
  
  # This results in a black background behind the legend - make it white
  g = g + theme(plot.margin       = margin(b = 2),
                plot.background   = element_rect(fill = "white", colour = "white"))
  
  # Save updated figure if desired
  if (!is.null(save))
    fig_save(o, g, save)
  
  return(g)
}

# ---------------------------------------------------------
# Specialist figure: sensitivty analysis bars
# ---------------------------------------------------------
plot_sensitivity = function(f, fig_name, plot_metric = "avert_diff", group_policy = TRUE) {
  
  # Descriptive sensitivity scenario labels
  bounds = c(s1u = "value: 0.4",      s1l = "value: 0.2", 
             s2u = "value: 95%",      s2l = "value: 60%", 
             s3u = "value: +1 month", s3l = "value: -1 month", 
             s4u = "value: 1.2",      s4l = "value: 0.8", 
             s5u = "value: 1.4",      s5l = "value: 1.1", 
             s6u = "value: 50%",      s6l = "value: 0%")
  
  # Bar width (1 => bars touch)
  bar_width = 0.96
  
  # Consistent y limit for all plots
  y_lim = list(avert_diff = 0.9, 
               avert_dose = 0.15)
  
  # Amonut to nudge text from upper y limit
  text_size = list(avert_diff = 1.75, 
                   avert_dose = 1.7)
  
  # Amonut to nudge text from upper y limit
  text_nudge = list(avert_diff = 0.08, 
                    avert_dose = 0.01)
  
  # ---- Policy grouping ----
  
  # For ungrouped policies extract x0-style strings
  if (group_policy == FALSE)
    pattern = "x[0-9]+"
  
  # Grouping policy requires different reg exp strings
  if (group_policy == TRUE) {
    pattern = "x[0,A,B]+"
    
    # We'll also want a different (and unique) colour scheme
    f$colours$sens_vec = c("peru", "gold")
  }
  
  # ---- Extract model outcomes ----
  
  # Initiate results list
  result_list = list()
  
  # Combine all sensitivity scenarios with the best estimate equivalents
  scenarios = c(f$scenarios$sens_vec, f$scenarios$variant_1y)
  
  # Loop through all sensitivity scenarios
  for (scenario in scenarios) {
    
    # Load data file for for this scenario
    result = try_load(o$pth$scenarios, scenario)
    
    # Scale aggregated metrics to per 100k per year
    scaler = f$person_days / result$input$population_size
    
    # Extract values of all desired metrics
    result_list[[scenario]] = result$cum_output %>%
      filter(metric %in% c(f$metrics$sens, "n_doses"),
             is.na(group)) %>%
      group_by(scenario, metric) %>%
      summarise(value = max(mean)  * scaler, 
                lower = max(lower) * scaler, 
                upper = max(upper) * scaler) %>%
      ungroup() %>%
      as.data.table()
  }
  
  # Combine results into single datatable
  result_df = rbindlist(result_list)
  
  # A little more work to do if grouping policies
  if (group_policy) {
    
    # We want xA-style rather than x0-style IDs, which we'll then group
    result_df$scenario = str_replace(result_df$scenario, "x[1,2]", "xA")
    result_df$scenario = str_replace(result_df$scenario, "x[3,4]", "xB")
  }
  
  # Easy access scenario ID names (also remove newline chars)
  id_names = setNames(str_remove(f$scenarios$name, "\n"), names(f$scenarios$name))
  
  # Convert these ID-names into a dataframe that can be joined
  id_df = data.table(id = names(id_names), name = id_names)
  
  # Combine results, and use ID-names to group params, bounds, and policies
  result_df = result_df %>%
    mutate(policy_id = str_extract(scenario, pattern), 
           param_id  = str_extract(scenario, "s[0-9]+"), 
           bound_dir = str_extract(scenario, "[l,u]+"), 
           bound_id  = paste(param_id, bound_dir, sep = "")) %>%
    left_join(id_df[, .(id, policy = name)], by = c("policy_id" = "id")) %>%
    left_join(id_df[, .(id, param  = name)], by = c("param_id"  = "id")) %>%
    left_join(data.table(bound_id = names(bounds), 
                         bound    = bounds), 
              by = "bound_id") %>%
    replace_na(list(param = id_names[["b"]])) %>%
    select(param, bound, bound_dir, policy, metric, value, lower, upper)
  
  # Number of non-baseline scenarios sensitivity scenarios
  n_scenarios = sum(grepl("s[0-9]+", names(id_names)))
  
  # ---- Calculate gains from booster strategies ----
  
  # All policy (or policy groups) we're plotting
  all_policy = id_names[grepl(pattern, names(id_names))]
  
  # Difference in metric(s) averted by more intensive booster policies
  data_df = result_df %>%
    group_by(param, bound, bound_dir, policy, metric) %>%
    summarise(val = min(value), 
              lb  = min(lower), 
              ub  = min(upper)) %>%
    ungroup() %>%
    # Total number of hospitalisations
    arrange(param, bound, factor(policy, levels = all_policy)) %>%
    group_by(param, bound, metric) %>%
    mutate(N_val = thou_sep(round(max(val), -2)), 
           N_lb  = thou_sep(round(max(lb), -1)),
           N_ub  = thou_sep(round(max(ub), -1)),
           n = paste0("n = ", N_val, " [", N_lb, "-", N_ub, "]")) %>%
    ungroup() %>%
    select(-N_val, -N_lb, -N_ub) %>%
    # Metric set 1: hospitalisations averted
    group_by(param, bound, metric) %>%
    mutate(avert_diff   = pmax(0, lag(val, default = first(val)) - val) / max(val),
           avert_diff_l = pmax(0, lag(lb,  default = first(lb))  - val) / max(val),
           avert_diff_u = pmax(0, lag(ub,  default = first(ub))  - val) / max(val), 
           to   = cumsum(avert_diff),
           from = lag(to)) %>%
    ungroup() %>%
    filter(policy != id_names[["x0"]]) %>%
    # Metric set 2: hospitalisations averted per dose
    group_by(param, bound, policy) %>%
    mutate(n_doses = val[metric == "n_doses"], 
           avert_dose   = avert_diff   / (n_doses / 1e4), 
           avert_dose_l = avert_diff_l / (n_doses / 1e4),
           avert_dose_u = avert_diff_u / (n_doses / 1e4)) %>%
    ungroup() %>%
    filter(metric != "n_doses") %>%
    mutate(policy = factor(policy, levels = all_policy[-1])) %>%
    # Select metric of interest
    select(param, bound, bound_dir, policy, metric, n, from, to, 
           value = one_of(plot_metric), 
           lower = one_of(paste0(plot_metric, "_l")), 
           upper = one_of(paste0(plot_metric, "_u"))) %>%
    as.data.table()
  
  # Set text positioning variables
  data_df$x = mean(as.integer(data_df$policy))
  data_df$y = y_lim[[plot_metric]] - text_nudge[[plot_metric]]
  
  # Extract baseline results
  base_df = data_df[is.na(bound), ] %>%
    mutate(y = y_lim[[plot_metric]] - text_nudge[[plot_metric]] / 2)
  
  # Parameters ordered by influence
  worst_case = data_df %>%
    filter(param  != id_names[["b"]]) %>%
    group_by(param, bound, bound_dir) %>%
    summarise(avert_total = sum(value)) %>%
    ungroup() %>%
    group_by(param) %>%
    slice_max(avert_total, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(avert_total) %>%
    select(param, bound, bound_dir) %>%
    as.data.table()
  
  # Best-case equivalent of worst case order
  best_case = data_df %>%
    filter(param  != id_names[["b"]],
           !bound %in% worst_case$bound) %>%
    arrange(factor(param, levels = worst_case$param)) %>%
    select(param, bound, bound_dir) %>%
    unique()
  
  # Apply this ordering to plotting data
  plot_df = rbind(worst_case, best_case) %>%
    left_join(data_df, by = qc(param, bound, bound_dir)) %>%
    filter(!is.na(bound)) %>%
    mutate(param  = paste0("<b>", param, "</b>"), 
           param  = paste(param, bound, sep = "<br>"))
  
  # Set factors so ordering is retained
  plot_df$param = factor(plot_df$param, levels = unique(plot_df$param))
  
  # ---- Create sensitivity plot ----
  
  # Standard plot as a waterfall
  if (plot_metric == "avert_diff") {
    
    # Create waterfall plot
    g1 = ggplot(plot_df) +  
      aes(x     = as.integer(policy),
          xmin  = as.integer(policy) - bar_width/2, 
          xmax  = as.integer(policy) + bar_width/2,
          y     = from + value / 2,
          ymin  = from, 
          ymax  = to,
          fill  = policy, 
          label = paste0(round(100 * value), "% [", 
                         round(100 * lower), "-", 
                         round(100 * upper), "]")) +
      geom_rect(colour = "black", size = 0.25) + 
      geom_text(size = text_size$avert_diff)
  }
  
  # Alternative 'per dose' plot as bars
  if (plot_metric == "avert_dose") {
    
    # Create bar plot
    g1 = ggplot(plot_df) +  
      aes(x     = as.integer(policy),
          y     = value,
          fill  = policy, 
          label = paste0(round(100 * value, 1), "%\n[", 
                         round(100 * lower, 1), "-", 
                         round(100 * upper, 1), "]")) +
      geom_bar(stat = "identity", colour = "black", size = 0.25) +
      geom_text(size = text_size$avert_dose, 
                nudge_y = 0.01)
  }
  
  # Set total text above each set of bars
  g1 = g1 + geom_text(aes(x = x, y = y, label = n), size = 2)
  
  # Set of labels for each facet (skipping 'A')
  facet_labels_idx = 1 : (n_scenarios * 2) + 1
  facet_labels = toupper(letters[facet_labels_idx])
  
  # Apply faceting and set labels
  g1 = g1 + facet_wrap(~param, nrow = 2) +
    tag_facets(tag_pool = facet_labels, tag_suffix = "")
  
  # Set colour scheme
  g1 = g1 + scale_fill_manual(values = f$colours$sens_vec)
  
  # Prettify y axis 
  g1 = g1 + scale_y_continuous(labels = percent, 
                               limits = c(0, y_lim[[plot_metric]]), 
                               expand = expansion(mult = c(0, 0)))
  
  # Prettify legend
  g1 = g1 + guides(fill = guide_legend(
    nrow  = ifelse(group_policy, 1, length(all_policy) - 1), 
    title = "Booster strategies (increasing intensity)", 
    title.position = "top", 
    title.hjust    = 0.5))
  
  # Prettify theme
  g1 = g1 + theme_classic() + 
    theme(strip.text    = element_text(size = 6, margin = margin(2, 2, 2, 2)),
          strip.text.x  = element_markdown(),
          axis.text.y   = element_text(size = 5),
          axis.text.x   = element_blank(),
          axis.ticks    = element_line(size = 0.25),
          axis.ticks.length = unit(0.1, "lines"),
          axis.ticks.x  = element_blank(),
          axis.title    = element_blank(),
          axis.line     = element_blank(),
          plot.margin   = margin(2, 2, 2, 2), 
          panel.grid.major.y = element_line(), 
          panel.grid.minor.y = element_line(),
          panel.border  = element_rect(size = 0.25, colour = "black", fill = NA),
          panel.spacing = unit(0.25, "lines"),
          strip.background = element_blank(), 
          legend.text   = element_text(size = 7),
          legend.title  = element_text(size = 8, face = "bold"),
          legend.margin = margin(2, 2, 2, 2),
          legend.position = "bottom", 
          legend.key.height = unit(1, "lines"),
          legend.key.width  = unit(1, "lines"),
          legend.box.background = element_rect(size = 0.25))
  
  # Remove legend title if and only if grouping - it's obvious for only 2 bars
  if (group_policy == TRUE)
    g1 = g1 + theme(legend.title = element_blank())
  
  # ---- Create baseline equivalent ----
  
  # Standard plot as a waterfall
  if (plot_metric == "avert_diff") {
    
    # Create baseline waterfall plot
    g2 = ggplot(base_df) +  
      aes(x    = as.integer(policy),
          xmin = as.integer(policy) - bar_width/2, 
          xmax = as.integer(policy) + bar_width/2,
          y    = from + value / 2,
          ymin = from, 
          ymax = to,
          fill = policy, 
          label = paste0(round(100 * value), "% [", 
                         round(100 * lower), "-", 
                         round(100 * upper), "]")) +
      geom_rect(colour = "black", size = 0.25) + 
      geom_text(size = text_size$avert_diff)
    
    # Use special colour set
    colours = f$colours$sens_vec
  }
  
  # Alternative 'per dose' plot as bars
  if (plot_metric == "avert_dose") {
    
    # Create baseline bar plot
    g2 = ggplot(base_df) +  
      aes(x = as.integer(policy), 
          y = value,
          fill = rev(policy), 
          label = paste0(round(100 * value, 1), "%\n[", 
                         round(100 * lower, 1), "-", 
                         round(100 * upper, 1), "]")) +
      geom_bar(stat = "identity", colour = "black", size = 0.25, show.legend = FALSE) + 
      geom_text(size = text_size$avert_dose, nudge_y = 0.005)
    
    # Use default policy colours
    colours = rev(f$colours$sens_vec)
  }
  
  # Add text to baseline plot
  g2 = g2 + geom_text(aes(x = x, y = y, label = n), size = 2)
  
  # Set and tag facets
  g2 = g2 + facet_wrap(~param, labeller = label_wrap_gen(16)) + 
    tag_facets(tag_levels = "A", tag_suffix = "")
  
  # Set colour scheme
  g2 = g2 + scale_fill_manual(values = colours)
  
  # Prettify y axis 
  g2 = g2 + scale_y_continuous(labels = percent, 
                               limits = c(0, y_lim[[plot_metric]]), 
                               expand = expansion(mult = c(0, 0)))
  
  # Set descriptive y axis label
  g2 = g2 + ylab("Incremental hospitalisations averted relative to no boosters")
  
  # Prettify theme
  g2 = g2 + theme_classic() + 
    theme(strip.text    = element_text(size = 6, margin = margin(2, 2, 2, 2)),
          strip.text.x  = element_text(face = "bold"),
          axis.text.y   = element_text(size = 5),
          axis.text.x   = element_blank(),
          axis.ticks    = element_line(size = 0.25),
          axis.ticks.length = unit(0.1, "lines"),
          axis.ticks.x  = element_blank(),
          axis.title.y  = element_text(size = 9),
          axis.title.x  = element_blank(),
          axis.line     = element_blank(),
          panel.grid.major.y = element_line(), 
          panel.grid.minor.y = element_line(),
          panel.border  = element_rect(size = 0.25, colour = "black", fill = NA),
          plot.margin   = margin(2, 2, 2, 2), 
          strip.background = element_blank())
  
  # ---- Shared legend ----
  
  # Copy the legend - we'll use this in combined plot
  shared_legend = get_legend(g1)
  
  # Now remove it from the sensitivity plot
  g1 = g1 + theme(legend.position = "none")
  
  # Use a caption to describe what 'n' refers to 
  g1 = g1 + labs(caption = "n = hospitalisations with no boosters") +
    theme(plot.caption = element_text(size = 6))
  
  # A blank equivalent for the baseline, simply so things align
  g2 = g2 + labs(caption = " ") +
    theme(plot.caption = element_text(size = 6))
  
  # ---- Combine plots ----
  
  # Combine the seasonality colourbars to the main figure
  g = ggarrange(g2, g1,
                nrow    = 1,
                ncol    = 2,
                widths  = c(0.24, 1),
                legend  = "bottom",
                legend.grob = shared_legend)
  
  # This results in a black background behind the legend - make it white
  g = g + theme(plot.margin     = margin(b = 2), 
                plot.background = element_rect(fill   = "white", 
                                               colour = "white"))
  
  # Save these figures to file
  fig_save(o, g, fig_name, width = 8)
}

# ---------------------------------------------------------
# Specialist figure: timing sensitivty analysis
# ---------------------------------------------------------
plot_timing_sensitivity = function(f, fig_name) {
  
  # Descriptive names for leegnd entries
  stat_dict = c("peak"  = "Objective: minimise peak hospitalisations", 
                "total" = "Objective: minimise total hospitalisations")
  
  # Two custom colours for optimal by peak and total
  colours = c("gold", "forestgreen")
  
  # Report optimal if within tolernace of true optimal
  optim_tol = 0.025
  
  # ---- Extract model outcomes ----
  
  # Short hand for sceanario matix
  scenarios   = f$scenarios$timing_mat
  n_scenarios = dim(scenarios)
  
  # Initiate results list
  optim_list = list()
  
  # Loop through all timing sensitivity scenarios
  for (i in 1 : n_scenarios[2]) {
    
    # Initiate results list
    result_list = list()
    
    # Loop through
    for (j in 1 : n_scenarios[1]) {
      
      # Load raw data file for for this scenario
      raw_file   = paste0(scenarios[j, i], "_raw")
      raw_result = try_load(o$pth$scenarios, raw_file)
      
      # Extract values of all desired metrics
      result_list[[j]] = raw_result %>%
        filter(metric %in% qc(hospital_beds, hospital_admissions), 
               is.na(group), 
               date > 90) %>%
        pivot_wider(names_from = metric)
    }
    
    # Determine peak and total, and filter out the non-optimal per seed
    optim_list[[i]] = rbindlist(result_list) %>%
      group_by(scenario, seed) %>%
      summarise(peak  = max(hospital_beds), 
                total = sum(hospital_admissions)) %>%
      ungroup() %>%
      pivot_longer(cols = c(peak, total), 
                   names_to = "stat") %>%
      group_by(stat, seed) %>%
      mutate(min = min(value) * (1 - optim_tol), 
             max = min(value) * (1 + optim_tol)) %>%
      ungroup() %>%
      filter(value >= min, 
             value <= max) %>%
      mutate(param = colnames(scenarios)[i]) %>%
      arrange(seed, stat, scenario)
  }
  
  # Total number of simulations (to calculate proportions)
  n_seeds = length(unique(raw_result$seed))
  n_sims  = n_seeds * n_scenarios[2]
  
  # Combine all optimal timings across all scenarios
  optim_df = rbindlist(optim_list) %>%
    mutate(month = str_extract(scenario, "m[0-9]+"), 
           month = as.numeric(substring(month, 2))) %>%
    group_by(stat, month) %>%
    summarise(n_optimal = n()) %>%
    ungroup() %>%
    mutate(p_optimal = n_optimal / n_sims) %>%
    as.data.table()
  
  # Expand the data to include trivial values: helps the smoothing
  plot_df = expand_grid(stat   = c("peak", "total"), 
                        month  = 1 : n_scenarios[1]) %>%
    left_join(optim_df, by = qc(stat, month)) %>%
    replace_na(list(n_optimal = 0, 
                    p_optimal = 0)) %>%
    mutate(stat = recode(stat, !!!stat_dict), 
           stat = factor(stat, stat_dict),
           diff_optim = lag(n_optimal, default = 0)) %>%
    filter(n_optimal + diff_optim > 0) %>%
    select(-diff_optim) %>%
    as.data.table()
  
  # ---- Create sensitivity plot ----
  
  # Degrees of freedom for the smoothing splines
  n_df = sum(table(plot_df$month) > 1) - 1
  
  # Create smooth density plot of optimal months
  g = ggplot(plot_df) +  
    aes(x = month, y = p_optimal, colour = stat, fill = stat) + 
    stat_smooth(geom     = "area", 
                position = "identity",
                method   = "glm", 
                formula  = y ~ ns(x, n_df), 
                alpha    = 0.3, 
                size     = 2)
  
  # Set colour scheme
  g = g + scale_fill_manual(values = colours) + 
    scale_colour_manual(values = colours)
  
  # Set descriptive title and subtitle
  g = g + labs(title    = "Optimal start for annual booster campaign",
               subtitle = "(considering parameter and stochastic uncertainty)")
  
  # Prettify y axis 
  g = g + scale_y_continuous(name   = "Percentage of simulations \n", 
                             labels = percent, 
                             limits = c(0, NA), 
                             expand = expansion(mult = c(0, 0.1)))
  # Prettify x axis 
  g = g + scale_x_continuous(limits = c(1, n_scenarios[1]), 
                             breaks = 1 : n_scenarios[1], 
                             labels = rownames(scenarios), 
                             expand = expansion(mult = c(0, 0)))
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(plot.title    = element_text(size = 32, hjust = 0.5),
          plot.subtitle = element_text(size = 20, hjust = 0.5),
          plot.margin   = margin(15, 15, 15, 15, "pt"),
          axis.text     = element_text(size = 16),
          axis.text.x   = element_text(hjust = 1, angle = 50),
          axis.title    = element_text(size = 24),
          axis.title.x  = element_blank(),
          axis.line     = element_blank(),
          panel.border  = element_rect(size = 1, colour = "black", fill = NA),
          panel.spacing = unit(2, "lines"),
          legend.text   = element_text(size = 16),
          legend.title  = element_blank(),
          legend.key    = element_blank(),
          legend.position = "bottom", 
          legend.key.height = unit(2, "lines"),
          legend.key.width  = unit(2, "lines"),
          legend.box.background = element_rect())
  
  # Save these figures to file
  fig_save(o, g, fig_name)
}

