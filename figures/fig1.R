# Input required: prok_latdist, euk_latdist (ALON Stage 2 final outputs)

library(dplyr)
library(ggplot2)


prok_latdist_count <- prok_latdist %>%
  mutate(
    # tropical MUST be monohemispherical
    source = if_else(
      final_geo == "tropical",
      "monohemispherical",
      source
    ),

    source = factor(
      source,
      levels = c(
        "bihemispherical",
        "monohemispherical"
      )
    ),

    final_geo = factor(
      final_geo,
      levels = c(
        # bihemispherical
        "bipolar",
        "subpolar",
        "high_latitude",
        "subtropical",

        # monohemispherical: north -> south
        "polar_north",
        "subpolar_north",
        "high_latitude_north",
        "subtropical_north",
        "tropical",
        "subtropical_south",
        "high_latitude_south",
        "subpolar_south",
        "polar_south"
      )
    ),

    # preferential first so we can put it on top
    excl_status = factor(
      excl_status,
      levels = c(
        "preferential",
        "exclusive"
      )
    )
  ) %>%
  filter(
    source %in% c("bihemispherical", "monohemispherical"),
    excl_status %in% c("preferential", "exclusive"),
    !is.na(final_geo)
  ) %>%
  count(
    source,
    final_geo,
    excl_status,
    name = "n"
  )

p_prok_latdist_count <- ggplot(
  prok_latdist_count,
  aes(
    x = final_geo,
    y = n,
    fill = excl_status
  )
) +
  geom_col(
    width = 0.8,
    position = position_stack(reverse = TRUE)
  ) +

  geom_text(
    aes(label = n),
    position = position_stack(
      vjust = 0.5,
      reverse = TRUE
    ),
    color = "white"
  ) +

  scale_y_continuous(
    labels = scales::comma
  ) +

  facet_grid(
    cols = vars(source),
    scales = "free_x",
    space = "free_x"
  ) +

  scale_fill_manual(
    values = c(
      preferential = "red",
      exclusive = "blue",
      cosmopolitan = "burlywood"
    ),
    breaks = c(
      "preferential",
      "exclusive"
    )
  ) +

  labs(
    x = NULL,
    y = "Number of taxa",
    fill = "Distribution type"
  ) +

  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    legend.position = "right"
  )

p_prok_latdist_count

euk_latdist_count <- euk_latdist %>%
  mutate(
    # tropical MUST be monohemispherical
    source = if_else(
      final_geo == "tropical",
      "monohemispherical",
      source
    ),
    
    source = factor(
      source,
      levels = c(
        "bihemispherical",
        "monohemispherical"
      )
    ),
    
    final_geo = factor(
      final_geo,
      levels = c(
        # bihemispherical
        "bipolar",
        "subpolar",
        "high_latitude",
        "subtropical",
        
        # monohemispherical: north -> south
        "polar_north",
        "subpolar_north",
        "high_latitude_north",
        "subtropical_north",
        "tropical",
        "subtropical_south",
        "high_latitude_south",
        "subpolar_south",
        "polar_south"
      )
    ),
    
    excl_status = factor(
      excl_status,
      levels = c(
        "preferential",
        "exclusive"
      )
    )
  ) %>%
  filter(
    source %in% c("bihemispherical", "monohemispherical"),
    excl_status %in% c("preferential", "exclusive"),
    !is.na(final_geo)
  ) %>%
  count(
    source,
    final_geo,
    excl_status,
    name = "n"
  ) %>%
  bind_rows(
    tibble(
      source = factor(
        c("monohemispherical", "monohemispherical"),
        levels = levels(.$source)
      ),
      final_geo = factor(
        c("subtropical_north", "subtropical_north"),
        levels = levels(.$final_geo)
      ),
      excl_status = factor(
        c("preferential", "exclusive"),
        levels = levels(.$excl_status)
      ),
      n = c(0, 0)
    )
  )

p_euk_latdist_count <- ggplot(
  euk_latdist_count,
  aes(
    x = final_geo,
    y = n,
    fill = excl_status
  )
) +
  geom_col(
    width = 0.8,
    position = position_stack(reverse = TRUE)
  ) +

  geom_text(
    aes(label = n),
    position = position_stack(
      vjust = 0.5,
      reverse = TRUE
    ),
    color = "white"
  ) +

  scale_y_continuous(
    labels = scales::comma
  ) +
  
  facet_grid(
    cols = vars(source),
    scales = "free_x",
    space = "free_x"
  ) +

  scale_fill_manual(
    values = c(
      preferential = "red",
      exclusive = "blue",
      cosmopolitan = "burlywood"
    ),
    breaks = c(
      "preferential",
      "exclusive"
    )
  ) +

  labs(
    x = NULL,
    y = "Number of taxa",
    fill = "Distribution type"
  ) +

  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    legend.position = "right"
  )

p_euk_latdist_count
