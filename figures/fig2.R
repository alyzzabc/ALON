# Input required: prok_latdist_tax, euk_latdist_tax (ALON Stage 2 final outputs including taxonomy);
# derep_prok_sw_epi, derep_euk_sw_epi (abundance observations)


library(dplyr)
library(stringr)
library(ggplot2)
library(scales)
library(Polychrome)

# All four categories that will always be shown in the plots
geo_plot <- c(
  "tropical",
  "high_latitude",
  "high_latitude_north",
  "high_latitude_south"
)

# ------------------------------------------------------------
# 1. Join prok abundance data to final_geo
# ------------------------------------------------------------

prok_geo <- derep_prok_sw_epi %>%
  left_join(
    prok_latdist_tax %>%
      select(taxon_id, final_geo) %>%
      distinct(),
    by = c("sequence" = "taxon_id"),
    relationship = "many-to-many"
  )


# ------------------------------------------------------------
# 2. Helper to extract an unambiguous prok genus
# ------------------------------------------------------------

extract_prok_genus <- function(x) {
  
  genus <- case_when(
    
    # only accept genus from a SINGLE taxonomy prediction
    !is.na(x) &
      str_count(x, fixed("Root;")) == 1 ~
      str_extract(x, "(?<=^|;\\s)g__[^;]+"),
    
    # multiple predictions / no taxonomy -> unusable
    TRUE ~ NA_character_
  )
  
  str_remove(genus, "^g__")
}


# ------------------------------------------------------------
# 3. Main helper
#
# rank_geos determines WHICH geographic groups are used
# to select the top 10.
#
# The resulting top 10 are ALWAYS plotted across all 4 groups.
# ------------------------------------------------------------

make_prok_rel_abund_plot <- function(rank_geos, plot_title = NULL) {
  
  # ----------------------------------------------------------
  # Find top 10 from the requested geographic union
  # ----------------------------------------------------------
  
  top10 <- prok_geo %>%
    filter(final_geo %in% rank_geos) %>%
    mutate(
      genus_raw = extract_prok_genus(taxonomy)
    ) %>%
    filter(!is.na(genus_raw)) %>%
    group_by(genus_raw) %>%
    summarise(
      total_abundance = sum(mean_norm_coverage, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(total_abundance)) %>%
    slice_head(n = 10)
  
  top10_names <- top10$genus_raw
  
  
  # ----------------------------------------------------------
  # Build plot data across ALL FOUR final_geo categories
  #
  # Everything not in the selected top 10 -> Other
  # Multiple taxonomy predictions -> Other
  # Missing genus -> Other
  # ----------------------------------------------------------
  
  plot_dat <- prok_geo %>%
    filter(final_geo %in% geo_plot) %>%
    mutate(
      genus_raw = extract_prok_genus(taxonomy),
      
      genus_plot = case_when(
        !is.na(genus_raw) &
          genus_raw %in% top10_names ~ genus_raw,
        
        TRUE ~ "Other"
      )
    ) %>%
    group_by(final_geo, genus_plot) %>%
    summarise(
      abundance = sum(mean_norm_coverage, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(final_geo) %>%
    mutate(
      relative_abundance =
        abundance / sum(abundance, na.rm = TRUE)
    ) %>%
    ungroup()
  
  
  # ----------------------------------------------------------
  # Factor order
  # ----------------------------------------------------------
  
  plot_levels <- c(top10_names, "Other")
  
  plot_dat <- plot_dat %>%
    mutate(
      genus_plot = factor(
        genus_plot,
        levels = rev(plot_levels)
      ),
      
      final_geo = factor(
        final_geo,
        levels = c(
          "tropical",
          "high_latitude",
          "high_latitude_north",
          "high_latitude_south"
        ),
        labels = c(
          "Tropical",
          "Bi-high latitude",
          "High latitude north",
          "High latitude south"
        )
      )
    )
  
  
  # ----------------------------------------------------------
  # Polychrome palette
  # Remove low-saturation / greyish colors
  # Other is the only grey
  # ----------------------------------------------------------
  
  poly_cols <- Polychrome::palette36.colors(36)
  hsv_cols <- rgb2hsv(col2rgb(poly_cols))
  
  poly_cols_no_grey <- poly_cols[
    hsv_cols["s", ] > 0.35
  ]
  
  pal_top10 <- poly_cols_no_grey[1:10]
  names(pal_top10) <- top10_names
  
  pal <- c(
    pal_top10,
    Other = "grey70"
  )
  
  
  # ----------------------------------------------------------
  # 4. Plot
  # ----------------------------------------------------------
  
  p <- ggplot(
    plot_dat,
    aes(
      x = final_geo,
      y = relative_abundance,
      fill = genus_plot
    )
  ) +
    geom_col(width = 0.7) +
    scale_fill_manual(
      values = pal,
      breaks = plot_levels
    ) +
    scale_y_continuous(
      labels = percent_format(accuracy = 1),
      expand = c(0, 0)
    ) +
    scale_x_discrete(drop = FALSE) +
    labs(
      title = plot_title,
      x = NULL,
      y = "Relative abundance",
      fill = "Genus"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = 11),
      legend.title = element_text(face = "bold")
    )
  
  
  # Return useful components
  list(
    plot = p,
    top10 = top10,
    plot_data = plot_dat
  )
}

prok_2 <- make_prok_rel_abund_plot(
  rank_geos = c(
    "tropical",
    "high_latitude",
    "high_latitude_north",
    "high_latitude_south"
  ),
  plot_title = "Top 10 from all four final_geo categories"
)

prok_2$plot
prok_2$top10


# ------------------------------------------------------------
# 1. Join euk abundance data to final_geo
# ------------------------------------------------------------

euk_geo <- derep_euk_sw_epi %>%
  left_join(
    euk_latdist_tax %>%
      select(taxon_id, final_geo) %>%
      distinct(),
    by = c("taxon" = "taxon_id"),
    relationship = "many-to-many"
  )


# ------------------------------------------------------------
# 2. Helper to extract an unambiguous genus
# ------------------------------------------------------------

extract_euk_genus <- function(x) {
  genus <- case_when(
    # multiple predictions -> unusable
    str_detect(x, fixed(",")) ~ NA_character_,
    
    # no pipe -> unusable
    !str_detect(x, fixed("|")) ~ NA_character_,
    
    # extract organism name after pipe,
    # then take genus before first underscore
    TRUE ~ x %>%
      str_remove("^.*\\|") %>%
      str_trim() %>%
      word(1) %>%
      str_remove("_.*$")
  )
  
  if_else(
    !is.na(genus) &
      genus != "" &
      str_detect(genus, "^[A-Za-z]"),
    genus,
    NA_character_
  )
}


# ------------------------------------------------------------
# 3. Function:
#    rank top 10 from selected geo union,
#    but plot those top 10 across all four categories
# ------------------------------------------------------------

make_euk_rel_abund_plot <- function(rank_geos, plot_title = NULL) {
  
  # ---- determine top 10 ----
  
  top10 <- euk_geo %>%
    filter(final_geo %in% rank_geos) %>%
    mutate(
      genus_raw = extract_euk_genus(taxon)
    ) %>%
    filter(!is.na(genus_raw)) %>%
    group_by(genus_raw) %>%
    summarise(
      total_abundance = sum(mean_norm_coverage, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(total_abundance)) %>%
    slice_head(n = 10)
  
  top10_names <- top10$genus_raw
  
  
  # ---- calculate relative abundance across all 4 final_geo ----
  
  plot_dat <- euk_geo %>%
    filter(final_geo %in% geo_plot) %>%
    mutate(
      genus_raw = extract_euk_genus(taxon),
      
      genus_plot = case_when(
        !is.na(genus_raw) &
          genus_raw %in% top10_names ~ genus_raw,
        TRUE ~ "Other"
      )
    ) %>%
    group_by(final_geo, genus_plot) %>%
    summarise(
      abundance = sum(mean_norm_coverage, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(final_geo) %>%
    mutate(
      relative_abundance = abundance / sum(abundance, na.rm = TRUE)
    ) %>%
    ungroup()
  
  
  # ---- factor order ----
  
  plot_levels <- c(top10_names, "Other")
  
  plot_dat <- plot_dat %>%
    mutate(
      genus_plot = factor(
        genus_plot,
        levels = rev(plot_levels)
      ),
      final_geo = factor(
        final_geo,
        levels = c(
          "tropical",
          "high_latitude",
          "high_latitude_north",
          "high_latitude_south"
        ),
        labels = c(
          "Tropical",
          "Bi-high latitude",
          "High latitude north",
          "High latitude south"
        )
      )
    )
  
  
  # ---- Polychrome palette, excluding grey-ish colors ----
  
  poly_cols <- Polychrome::palette36.colors(36)
  hsv_cols <- rgb2hsv(col2rgb(poly_cols))
  
  poly_cols_no_grey <- poly_cols[
    hsv_cols["s", ] > 0.35
  ]
  
  pal_top10 <- poly_cols_no_grey[1:10]
  names(pal_top10) <- top10_names
  
  pal <- c(
    pal_top10,
    Other = "grey70"
  )
  
  
  # ---- plot ----
  
  p <- ggplot(
    plot_dat,
    aes(
      x = final_geo,
      y = relative_abundance,
      fill = genus_plot
    )
  ) +
    geom_col(width = 0.7) +
    scale_fill_manual(
      values = pal,
      breaks = plot_levels
    ) +
    scale_y_continuous(
      labels = percent_format(accuracy = 1),
      expand = c(0, 0)
    ) +
    scale_x_discrete(drop = FALSE) +
    labs(
      title = plot_title,
      x = NULL,
      y = "Relative abundance",
      fill = "Genus"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = 11),
      legend.title = element_text(face = "bold")
    )
  
  list(
    plot = p,
    top10 = top10,
    plot_data = plot_dat
  )
}

euk_2 <- make_euk_rel_abund_plot(
  rank_geos = c(
    "tropical",
    "high_latitude",
    "high_latitude_north",
    "high_latitude_south"
  ),
  plot_title = "Top 10 from all four final_geo categories"
)

euk_2$plot
euk_2$top10
