# Input required: prok_latdist_tax, euk_latdist_tax (ALON Stage 2 final outputs including taxonomy);
# derep_prok_sw_epi, derep_euk_sw_epi (abundance observations -- must include lat and long);
# proks_gam, euks_gam (GAM predictions)


library(data.table)

## Find most abundant OTUs per lat distribution

## Proks
# make sure both are data.tables
setDT(prok_latdist_tax)
setDT(derep_prok_sw_epi)

# Preferential high-latitude taxa
prok_pref_highlat <- unique(
  prok_latdist_tax[
    excl_status == "preferential" &
      final_geo == "high_latitude",
    taxon_id
  ]
)

# Rank them by total abundance
prok_highlat_abundance <- derep_prok_sw_epi[
  sequence %in% prok_pref_highlat,
  .(
    total_abundance = sum(mean_norm_coverage, na.rm = TRUE),
    n_occurrences = .N
  ),
  by = sequence
][order(-total_abundance)]

prok_highlat_abundance

prok_highlat_abundance[1]

prok_pref_tropical <- unique(
  prok_latdist_tax[
    excl_status == "preferential" &
      final_geo == "tropical",
    taxon_id
  ]
)

prok_tropical_abundance <- derep_prok_sw_epi[
  sequence %in% prok_pref_tropical,
  .(
    total_abundance = sum(mean_norm_coverage, na.rm = TRUE),
    n_occurrences = .N
  ),
  by = sequence
][order(-total_abundance)]

prok_tropical_abundance

prok_tropical_abundance[1]


## Euks
# make sure both are data.tables
setDT(euk_latdist_tax)
setDT(derep_euk_sw_epi)

# Preferential high-latitude taxa
euk_pref_highlat <- unique(
  euk_latdist_tax[
    excl_status == "preferential" &
      final_geo == "high_latitude",
    taxon_id
  ]
)

# Rank them by total abundance
euk_highlat_abundance <- derep_euk_sw_epi[
  taxon %in% euk_pref_highlat,
  .(
    total_abundance = sum(mean_norm_coverage, na.rm = TRUE),
    n_occurrences = .N
  ),
  by = taxon
][order(-total_abundance)]

euk_highlat_abundance

euk_highlat_abundance[1]

euk_pref_tropical <- unique(
  euk_latdist_tax[
    excl_status == "preferential" &
      final_geo == "tropical",
    taxon_id
  ]
)

euk_tropical_abundance <- derep_euk_sw_epi[
  taxon %in% euk_pref_tropical,
  .(
    total_abundance = sum(mean_norm_coverage, na.rm = TRUE),
    n_occurrences = .N
  ),
  by = taxon
][order(-total_abundance)]

euk_tropical_abundance

euk_tropical_abundance[1]


### Plot per OTU
# ============================================================
# Plot ONE exact prokaryotic taxon_id / sequence
# ============================================================

plot_exact_taxon_map <- function(
  derep_df,
  taxon_id,
  coverage_col = "mean_norm_coverage",
  taxon_col = "sequence",
  run_col = "rep_Run",
  lat_col = "manual_latitude",
  lon_col = "manual_longitude",
  filter_col = NULL,
  filter_value = NULL,
  aggregate_fun = sum,
  world_fill = "grey95",
  world_line = "grey70",
  ring_color = "grey30",
  point_color = "steelblue",
  point_alpha = 0.5,
  ring_size = 1,
  ring_stroke = 0.3,
  show_zeroes = FALSE,
  size_limits = NULL,
  size_breaks = NULL,
  size_range = c(1, 8),
  plot_title = NULL
) {

  library(dplyr)
  library(ggplot2)


  # ----------------------------------------------------------
  # Optional filtering
  # e.g. filter_col = "layer", filter_value = "epi"
  # ----------------------------------------------------------

  derep_df_filt <- derep_df

  if (!is.null(filter_col) && !is.null(filter_value)) {
    derep_df_filt <- derep_df_filt %>%
      dplyr::filter(.data[[filter_col]] %in% filter_value)
  }


  # ----------------------------------------------------------
  # Check that taxon exists
  # ----------------------------------------------------------

  if (!taxon_id %in% derep_df_filt[[taxon_col]]) {
    stop(
      "Taxon ID not found in column '",
      taxon_col,
      "': ",
      taxon_id
    )
  }


  # ----------------------------------------------------------
  # All sampled runs
  # Used to show zero/absent locations
  # ----------------------------------------------------------

  all_runs <- derep_df_filt %>%
    dplyr::select(
      all_of(c(run_col, lat_col, lon_col))
    ) %>%
    dplyr::rename(
      rep_Run = all_of(run_col),
      manual_latitude = all_of(lat_col),
      manual_longitude = all_of(lon_col)
    ) %>%
    dplyr::distinct() %>%
    dplyr::filter(
      !is.na(manual_latitude),
      !is.na(manual_longitude)
    )


  # ----------------------------------------------------------
  # Select exact taxon
  # ----------------------------------------------------------

  matched_df <- derep_df_filt %>%
    dplyr::filter(
      .data[[taxon_col]] == taxon_id
    ) %>%
    dplyr::filter(
      !is.na(.data[[lat_col]]),
      !is.na(.data[[lon_col]])
    )


  # ----------------------------------------------------------
  # Aggregate coverage per run/location
  # ----------------------------------------------------------

  detected_runs <- matched_df %>%
    dplyr::group_by(
      rep_Run = .data[[run_col]],
      manual_latitude = .data[[lat_col]],
      manual_longitude = .data[[lon_col]]
    ) %>%
    dplyr::summarise(
      plot_coverage = aggregate_fun(
        .data[[coverage_col]],
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    dplyr::filter(
      is.finite(plot_coverage),
      plot_coverage > 0
    ) %>%
    dplyr::arrange(
      dplyr::desc(plot_coverage)
    )


  # ----------------------------------------------------------
  # Sampled locations without this taxon
  # ----------------------------------------------------------

  absent_runs <- all_runs %>%
    dplyr::anti_join(
      detected_runs %>%
        dplyr::select(rep_Run),
      by = "rep_Run"
    )


  # ----------------------------------------------------------
  # World map
  # ----------------------------------------------------------

  world <- ggplot2::map_data("world")


  p <- ggplot() +

    geom_polygon(
      data = world,
      aes(
        x = long,
        y = lat,
        group = group
      ),
      fill = world_fill,
      color = world_line,
      linewidth = 0.2
    ) +

    geom_hline(
      yintercept = c(
        -90, -75, -60, -45, -30, -15,
        0,
        15, 30, 45, 60, 75, 90
      ),
      color = "grey90",
      linewidth = 0.3
    ) +

    geom_vline(
      xintercept = c(
        -120, -60, 0, 60, 120
      ),
      color = "grey90",
      linewidth = 0.3
    ) +

    geom_point(
      data = detected_runs,
      aes(
        x = manual_longitude,
        y = manual_latitude,
        size = plot_coverage
      ),
      color = point_color,
      alpha = point_alpha,
      stroke = 0
    ) +

    coord_quickmap(
      xlim = c(-180.5, 180.5),
      ylim = c(-90.5, 90.5),
      expand = TRUE
    ) +

    scale_x_continuous(
      breaks = c(
        -120, -60, 0, 60, 120
      )
    ) +

    scale_y_continuous(
      breaks = c(
        -90, -75, -60, -45, -30, -15,
        0,
        15, 30, 45, 60, 75, 90
      )
    ) +

    scale_size_continuous(
      name = coverage_col,
      limits = size_limits,
      breaks = size_breaks,
      range = size_range
    ) +

    labs(
      title = if (is.null(plot_title)) {
        paste0(
          "Distribution of taxon: ",
          taxon_id
        )
      } else {
        plot_title
      },
      x = "Longitude",
      y = "Latitude"
    ) +

    theme_bw() +

    theme(
      panel.grid = element_blank(),
      axis.text = element_text(
        color = "black"
      ),
      panel.border = element_blank()
    )


  # ----------------------------------------------------------
  # Optional zero/absent locations
  # ----------------------------------------------------------

  if (show_zeroes) {

    p <- p +

      geom_point(
        data = absent_runs,
        aes(
          x = manual_longitude,
          y = manual_latitude
        ),
        shape = 1,
        color = ring_color,
        stroke = ring_stroke,
        size = ring_size
      )
  }


  p
}

# ============================================================
# Get shared size limits for exact taxon IDs
# ============================================================

get_exact_taxa_shared_size_limits <- function(
  derep_df,
  taxon_ids,
  coverage_col = "mean_norm_coverage",
  taxon_col = "sequence",
  run_col = "rep_Run",
  lat_col = "manual_latitude",
  lon_col = "manual_longitude",
  filter_col = NULL,
  filter_value = NULL,
  aggregate_fun = sum
) {

  library(dplyr)


  # ----------------------------------------------------------
  # Optional filtering
  # ----------------------------------------------------------

  derep_df_filt <- derep_df

  if (!is.null(filter_col) && !is.null(filter_value)) {

    derep_df_filt <- derep_df_filt %>%
      dplyr::filter(
        .data[[filter_col]] %in% filter_value
      )
  }


  # ----------------------------------------------------------
  # Check which requested taxa exist
  # ----------------------------------------------------------

  found_ids <- intersect(
    taxon_ids,
    unique(derep_df_filt[[taxon_col]])
  )

  missing_ids <- setdiff(
    taxon_ids,
    found_ids
  )


  if (length(missing_ids) > 0) {

    warning(
      "The following taxon IDs were not found:\n",
      paste(missing_ids, collapse = "\n")
    )
  }


  if (length(found_ids) == 0) {

    stop(
      "None of the requested taxon IDs were found in ",
      taxon_col
    )
  }


  # ----------------------------------------------------------
  # Get requested taxa
  # ----------------------------------------------------------

  matched_df <- derep_df_filt %>%
    dplyr::filter(
      .data[[taxon_col]] %in% found_ids
    ) %>%
    dplyr::filter(
      !is.na(.data[[lat_col]]),
      !is.na(.data[[lon_col]])
    )


  # ----------------------------------------------------------
  # IMPORTANT:
  # calculate coverage separately for each taxon/run
  #
  # This means the maximum represents the maximum point
  # appearing on EITHER individual taxon map.
  # ----------------------------------------------------------

  detected_runs <- matched_df %>%
    dplyr::group_by(
      taxon = .data[[taxon_col]],
      rep_Run = .data[[run_col]],
      manual_latitude = .data[[lat_col]],
      manual_longitude = .data[[lon_col]]
    ) %>%
    dplyr::summarise(
      plot_coverage = aggregate_fun(
        .data[[coverage_col]],
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    dplyr::filter(
      is.finite(plot_coverage),
      plot_coverage > 0
    )


  if (nrow(detected_runs) == 0) {

    return(c(0, 1))
  }


  max_cov <- max(
    detected_runs$plot_coverage,
    na.rm = TRUE
  )


  c(0, max_cov)
}

### Usage:

#taxon_1 <- "TGGTCAAGAAAATCAACTATAATTCCTGAATTTATTGGTGTTAGTTTTTTGATTTATAAT"
#taxon_2 <- "TGGTCTCGTCGTTCTATGATCATCCCAGATATGATTGGGTTGACCATTGCGGTCCATAAC"

taxon_1 <- "1848828|Micromonas_polaris"
taxon_2 <- "1606511|Chloropicon_mariensis"

## Plot for proks
taxon_limits <- get_exact_taxa_shared_size_limits(
  derep_df = derep_prok_sw_epi,
  taxon_col = "sequence",
  taxon_ids = c(
    taxon_1,
    taxon_2
  )
)

taxon_breaks <- pretty(
  taxon_limits,
  n = 4
)

taxon_breaks <- taxon_breaks[
  taxon_breaks > 0
]

p_taxon_1 <- plot_exact_taxon_map(
  derep_df = derep_prok_sw_epi,
  taxon_col = "sequence",
  taxon_id = taxon_1,
  size_limits = taxon_limits,
  size_breaks = taxon_breaks,
  point_color = "olivedrab",
  show_zeroes = TRUE,
  plot_title = "Most abundant pref. high lat. (Pelagibacter)"
)


p_taxon_2 <- plot_exact_taxon_map(
  derep_df = derep_prok_sw_epi,
  taxon_col = "sequence",
  taxon_id = taxon_2,
  size_limits = taxon_limits,
  size_breaks = taxon_breaks,
  point_color = "olivedrab",
  show_zeroes = TRUE,
  plot_title = "Most abundant pref. trop. (Prochlorococcus A)"
)

p_taxon_1/p_taxon_2

## Plot for euks
taxon_limits <- get_exact_taxa_shared_size_limits(
  derep_df = derep_euk_sw_epi,
  taxon_col = "taxon",
  taxon_ids = c(
    taxon_1,
    taxon_2
  )
)

taxon_breaks <- pretty(
  taxon_limits,
  n = 4
)

taxon_breaks <- taxon_breaks[
  taxon_breaks > 0
]

p_taxon_1 <- plot_exact_taxon_map(
  derep_df = derep_euk_sw_epi,
  taxon_col = "taxon",
  taxon_id = taxon_1,
  size_limits = taxon_limits,
  size_breaks = taxon_breaks,
  point_color = "olivedrab",
  show_zeroes = TRUE,
  plot_title = "1848828|Micromonas_polaris"
)


p_taxon_2 <- plot_exact_taxon_map(
  derep_df = derep_euk_sw_epi,
  taxon_col = "taxon",
  taxon_id = taxon_2,
  size_limits = taxon_limits,
  size_breaks = taxon_breaks,
  point_color = "olivedrab",
  show_zeroes = TRUE,
  plot_title = "Taxon2"
)

p_taxon_1/p_taxon_2


##### Building maps with independent point size scales
get_exact_taxon_size_limits <- function(
  derep_df,
  taxon_id,
  coverage_col = "mean_norm_coverage",
  taxon_col = "sequence",
  run_col = "rep_Run",
  lat_col = "manual_latitude",
  lon_col = "manual_longitude",
  aggregate_fun = sum
) {

  library(dplyr)

  detected_runs <- derep_df %>%
    dplyr::filter(.data[[taxon_col]] == taxon_id) %>%
    dplyr::filter(
      !is.na(.data[[lat_col]]),
      !is.na(.data[[lon_col]])
    ) %>%
    dplyr::group_by(
      rep_Run = .data[[run_col]],
      manual_latitude = .data[[lat_col]],
      manual_longitude = .data[[lon_col]]
    ) %>%
    dplyr::summarise(
      plot_coverage = aggregate_fun(
        .data[[coverage_col]],
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    dplyr::filter(
      is.finite(plot_coverage),
      plot_coverage > 0
    )

  if (nrow(detected_runs) == 0) {
    return(c(0, 1))
  }

  c(
    0,
    max(detected_runs$plot_coverage, na.rm = TRUE)
  )
}

taxon_1_limits <- get_exact_taxon_size_limits(
  derep_df = derep_euk_sw_epi,
  taxon_col = "taxon",
  taxon_id = taxon_1
)

taxon_1_breaks <- pretty(
  taxon_1_limits,
  n = 4
)

taxon_1_breaks <- taxon_1_breaks[
  taxon_1_breaks > 0 &
  taxon_1_breaks <= taxon_1_limits[2]
]

taxon_2_limits <- get_exact_taxon_size_limits(
  derep_df = derep_euk_sw_epi,
  taxon_col = "taxon",
  taxon_id = taxon_2
)

taxon_2_breaks <- pretty(
  taxon_2_limits,
  n = 4
)

taxon_2_breaks <- taxon_2_breaks[
  taxon_2_breaks > 0 &
  taxon_2_breaks <= taxon_2_limits[2]
]

p_taxon_1 <- plot_exact_taxon_map(
  derep_df = derep_euk_sw_epi,
  taxon_col = "taxon",
  taxon_id = taxon_1,
  size_limits = taxon_1_limits,
  size_breaks = taxon_1_breaks,
  point_color = "olivedrab",
  show_zeroes = TRUE,
  plot_title = "1848828|Micromonas_polaris"
)


p_taxon_2 <- plot_exact_taxon_map(
  derep_df = derep_euk_sw_epi,
  taxon_col = "taxon",
  taxon_id = taxon_2,
  size_limits = taxon_2_limits,
  size_breaks = taxon_2_breaks,
  point_color = "olivedrab",
  show_zeroes = TRUE,
  plot_title = "1606511|Chloropicon_mariensis"
)


p_taxon_1 / p_taxon_2

######################################################
### To add the GAM curve side by side with the map ###
######################################################

plot_map_with_aligned_gam <- function(
  map_plot,
  gam_df,
  taxon_id_value,
  gam_width = 75,
  gam_gap = 15,
  gam_line_color = "darkgreen",
  gam_ribbon_fill = "darkolivegreen3",
  gam_ribbon_alpha = 0.3,
  gam_line_size = 0.8,
  gam_breaks = NULL,
  backtransform_log1p = TRUE,
  map_breaks = c(-120, -60, 0, 60, 120),
  y_breaks = seq(-90, 90, 15),
  combined_xlab = "Longitude                                      Predicted mean normalized coverage"
) {

  library(ggplot2)
  library(data.table)

  # ----------------------------------------------------------
  # Get GAM predictions for the requested taxon
  # ----------------------------------------------------------

  gam_dat <- gam_df[
    taxon_id == taxon_id_value,
    .(latitude, fit, lo, hi)
  ]

  if (nrow(gam_dat) == 0) {
    stop("No GAM predictions found for taxon: ", taxon_id_value)
  }

  # Sort by latitude so the curve draws correctly
  data.table::setorder(gam_dat, latitude)

  # ----------------------------------------------------------
  # Back-transform from log1p scale if requested
  # ----------------------------------------------------------

  if (backtransform_log1p) {
    gam_dat[, `:=`(
      fit_plot = pmax(expm1(fit), 0),
      lo_plot  = pmax(expm1(lo), 0),
      hi_plot  = pmax(expm1(hi), 0)
    )]
  } else {
    gam_dat[, `:=`(
      fit_plot = pmax(fit, 0),
      lo_plot  = pmax(lo, 0),
      hi_plot  = pmax(hi, 0)
    )]
  }

  # ----------------------------------------------------------
  # Where the GAM strip sits to the right of the world map
  # ----------------------------------------------------------

  gam_start <- 180 + gam_gap
  gam_end   <- gam_start + gam_width

  # ----------------------------------------------------------
  # Choose GAM breaks
  # These are in abundance units (same scale as map abundance)
  # ----------------------------------------------------------

  if (is.null(gam_breaks)) {
    gam_max_tmp <- max(gam_dat$hi_plot, na.rm = TRUE)
    gam_breaks <- pretty(c(0, gam_max_tmp), n = 4)
    gam_breaks <- gam_breaks[gam_breaks >= 0]
  }

  gam_breaks <- sort(unique(gam_breaks))

  # Use the larger of:
  # - actual GAM upper CI
  # - requested breaks
  gam_max <- max(c(gam_dat$hi_plot, gam_breaks), na.rm = TRUE)

  if (gam_max <= 0 || !is.finite(gam_max)) {
    stop("GAM maximum is not valid after transformation.")
  }

  # ----------------------------------------------------------
  # Rescale GAM values into fake longitude coordinates
  # ----------------------------------------------------------

  gam_dat[, `:=`(
    fit_scaled = gam_start + (fit_plot / gam_max) * gam_width,
    lo_scaled  = gam_start + (lo_plot  / gam_max) * gam_width,
    hi_scaled  = gam_start + (hi_plot  / gam_max) * gam_width
  )]

  # Position of GAM tick marks in fake longitude coordinates
  gam_break_positions <- gam_start + (gam_breaks / gam_max) * gam_width

  # ----------------------------------------------------------
  # Build combined x-axis:
  # left side = map longitude
  # right side = GAM abundance ticks
  # ----------------------------------------------------------

  all_x_breaks <- c(map_breaks, gam_break_positions)
  all_x_labels <- c(as.character(map_breaks), as.character(gam_breaks))

  # ----------------------------------------------------------
  # Add GAM to map
  # ----------------------------------------------------------

  p <- map_plot +

    # Ribbon
    geom_ribbon(
      data = gam_dat,
      aes(
        y = latitude,
        xmin = lo_scaled,
        xmax = hi_scaled
      ),
      inherit.aes = FALSE,
      orientation = "y",
      fill = gam_ribbon_fill,
      alpha = gam_ribbon_alpha,
      linewidth = 0
    ) +

    # GAM line
    geom_path(
      data = gam_dat,
      aes(
        x = fit_scaled,
        y = latitude
      ),
      inherit.aes = FALSE,
      color = gam_line_color,
      linewidth = gam_line_size
    ) +

    # Separator between map and GAM strip
    geom_vline(
      xintercept = gam_start,
      color = "grey50",
      linewidth = 0.4
    ) +

    # Replace coord/axes so the map and GAM share latitude perfectly
    coord_quickmap(
      xlim = c(-180, gam_end),
      ylim = c(-90, 90),
      expand = FALSE
    ) +

    scale_y_continuous(
      breaks = y_breaks,
      limits = c(-90, 90),
      expand = c(0, 0)
    ) +

    scale_x_continuous(
      breaks = all_x_breaks,
      labels = all_x_labels,
      limits = c(-180, gam_end),
      expand = c(0, 0)
    ) +

    labs(
      x = combined_xlab,
      y = "Latitude"
    )

  p
}

## For proks
library(data.table)

setDT(proks_gam)
setindex(proks_gam, taxon_id)

p_taxon_1_aligned <- plot_map_with_aligned_gam(
  map_plot = p_taxon_1,
  gam_df = proks_gam,
  taxon_id_value = taxon_1,
  backtransform_log1p = TRUE,
  gam_breaks = c(0, 25, 50, 75, 100)
)

p_taxon_1_aligned

p_taxon_2_aligned <- plot_map_with_aligned_gam(
  map_plot = p_taxon_2,
  gam_df = proks_gam,
  taxon_id_value = taxon_2,
  backtransform_log1p = TRUE,
  gam_breaks = c(0, 0.25, 0.5, 0.75, 1)
)

p_taxon_2_aligned

## For euks
setDT(euks_gam)
setindex(euks_gam, taxon_id)

p_taxon_1_aligned <- plot_map_with_aligned_gam(
  map_plot = p_taxon_1,
  gam_df = euks_gam,
  taxon_id_value = taxon_1,
  backtransform_log1p = TRUE,
  gam_breaks = c(0, 0.25, 0.5)
)

p_taxon_1_aligned

p_taxon_2_aligned <- plot_map_with_aligned_gam(
  map_plot = p_taxon_2,
  gam_df = euks_gam,
  taxon_id_value = taxon_2,
  backtransform_log1p = TRUE,
  gam_breaks = c(0, 0.025, 0.05)
)

p_taxon_2_aligned
