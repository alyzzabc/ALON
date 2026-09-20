## ============================================================
## reconcile_size_fractions.R
##
## Reconciles bihemispherical and monohemispherical classifications
## across multiple size fractions.
##
## Expected objects from run_biogeography.R:
##   classification_dirs
##   size_fraction_names
##   output_dir
##
## Expected raw observation object, one of:
##   reconcile_observations_input
##   observations_input
##
## Expected observation columns:
##   taxon_id, latitude, norm_coverage
##
## Optional objects:
##   reconciliation_rule
##     - "min_two_support"  default
##     - "strict_majority"
##
##   map_bihem_polar_to_bipolar
##     - FALSE default
##     - TRUE if you want polar labels renamed to bipolar
##
## Output:
##   <output_dir>/reconciled/
##     reconciled_bihemispherical.rds
##     reconciled_bihemispherical.tsv
##     reconciled_monohemispherical.rds
##     reconciled_monohemispherical.tsv
##     reconciled_biogeography_combined.rds
##     reconciled_biogeography_combined.tsv
##     reconciled_biogeography_combined_split.rds
##     reconciled_biogeography_combined_split.tsv
##     all_size_fraction_classifications.tsv
## ============================================================

suppressPackageStartupMessages({
  library(data.table)
})

## ============================================================
## 1. Required objects
## ============================================================

required_objects <- c(
  "classification_dirs",
  "size_fraction_names",
  "output_dir"
)

missing_objects <- required_objects[!vapply(required_objects, exists, logical(1))]

if (length(missing_objects) > 0L) {
  stop(
    "reconcile_size_fractions.R is missing required object(s): ",
    paste(missing_objects, collapse = ", "),
    call. = FALSE
  )
}

if (length(classification_dirs) != length(size_fraction_names)) {
  stop(
    "classification_dirs and size_fraction_names must have the same length.",
    "\nlength(classification_dirs): ", length(classification_dirs),
    "\nlength(size_fraction_names): ", length(size_fraction_names),
    call. = FALSE
  )
}

if (exists("reconcile_observations_input")) {
  obs_input <- reconcile_observations_input
} else if (exists("observations_input")) {
  obs_input <- observations_input
} else {
  stop(
    "Need either reconcile_observations_input or observations_input ",
    "for recomputing exclusivity.",
    call. = FALSE
  )
}

obs_input <- as.data.table(obs_input)

required_obs_cols <- c("taxon_id", "latitude", "norm_coverage")
missing_obs_cols <- setdiff(required_obs_cols, names(obs_input))

if (length(missing_obs_cols) > 0L) {
  stop(
    "Observation table is missing required column(s): ",
    paste(missing_obs_cols, collapse = ", "),
    call. = FALSE
  )
}

## Default reconciliation rule.
##
## min_two_support:
##   1 real source  -> keep all labels from that source
##   2 real sources -> label must appear in both
##   3 real sources -> label must appear in at least 2
##   >3 real sources -> label must appear in at least 2
##
## strict_majority:
##   1 real source  -> keep all labels from that source
##   2 real sources -> label must appear in both
##   3 real sources -> label must appear in at least 2
##   4 real sources -> label must appear in at least 3
##   5 real sources -> label must appear in at least 3
if (!exists("reconciliation_rule")) {
  reconciliation_rule <- "min_two_support"
}

allowed_reconciliation_rules <- c("min_two_support", "strict_majority")

if (!reconciliation_rule %in% allowed_reconciliation_rules) {
  stop(
    "Unknown reconciliation_rule: ", reconciliation_rule,
    ". Use one of: ",
    paste(allowed_reconciliation_rules, collapse = ", "),
    call. = FALSE
  )
}

if (length(size_fraction_names) > 3L && reconciliation_rule == "min_two_support") {
  warning(
    "More than three size fractions were provided. ",
    "The default reconciliation_rule is 'min_two_support', meaning any label ",
    "supported by at least two real size fractions is retained. ",
    "For >3 size fractions, inspect the reconciled output carefully or set ",
    "reconciliation_rule <- 'strict_majority' before sourcing this script."
  )
}

## Keep classifier-emitted labels as-is.
## Your current bihem output appears to use "polar" rather than "bipolar".
## If you later want polar -> bipolar in reconciled bihem outputs,
## set map_bihem_polar_to_bipolar <- TRUE before sourcing this script.
if (!exists("map_bihem_polar_to_bipolar")) {
  map_bihem_polar_to_bipolar <- FALSE
}

discordant_label <- "discordant"

message("Reconciliation rule: ", reconciliation_rule)

## ============================================================
## 2. Helpers
## ============================================================

clean_fraction_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  tolower(x)
}

strip_excl_prefix <- function(x) {
  sub("^(exclusive_|preferential_)", "", x)
}

collapse_pipe <- function(x) {
  x <- unique(na.omit(as.character(x)))
  x <- x[x != ""]
  if (length(x) == 0L) {
    NA_character_
  } else {
    paste(sort(x), collapse = "|")
  }
}

apply_reconciliation_rule <- function(geo_support, reconciliation_rule) {
  geo_support <- as.data.table(copy(geo_support))
  
  if (reconciliation_rule == "min_two_support") {
    geo_support[
      ,
      keep := fifelse(
        n_real_sources == 1L,
        TRUE,
        n_sources >= 2L
      )
    ]
    
  } else if (reconciliation_rule == "strict_majority") {
    geo_support[
      ,
      keep := fifelse(
        n_real_sources == 1L,
        TRUE,
        n_sources > n_real_sources / 2
      )
    ]
    
  } else {
    stop(
      "Unknown reconciliation_rule: ", reconciliation_rule,
      call. = FALSE
    )
  }
  
  geo_support[]
}

read_one_fraction <- function(classification_dir, size_fraction) {
  temp_dir <- file.path(classification_dir, "temp")
  
  bihem_file <- file.path(temp_dir, "bihemispherical_vOTUs.rds")
  mono_file  <- file.path(temp_dir, "monohemispherical_vOTUs.rds")
  
  out <- list()
  
  if (file.exists(bihem_file)) {
    bihem <- as.data.table(readRDS(bihem_file))
    
    if (!"tag" %in% names(bihem) && "final_tag" %in% names(bihem)) {
      bihem[, tag := final_tag]
    }
    
    required_cols <- c("taxon_id", "tag")
    missing_cols <- setdiff(required_cols, names(bihem))
    
    if (length(missing_cols) > 0L) {
      stop(
        "Bihemispherical file is missing required column(s): ",
        paste(missing_cols, collapse = ", "),
        "\nFile: ", bihem_file,
        call. = FALSE
      )
    }
    
    bihem[, source := "bihemispherical"]
    bihem[, size_fraction := size_fraction]
    out[["bihemispherical"]] <- bihem
  }
  
  if (file.exists(mono_file)) {
    mono <- as.data.table(readRDS(mono_file))
    
    if (!"tag" %in% names(mono) && "final_tag" %in% names(mono)) {
      mono[, tag := final_tag]
    }
    
    required_cols <- c("taxon_id", "tag")
    missing_cols <- setdiff(required_cols, names(mono))
    
    if (length(missing_cols) > 0L) {
      stop(
        "Monohemispherical file is missing required column(s): ",
        paste(missing_cols, collapse = ", "),
        "\nFile: ", mono_file,
        call. = FALSE
      )
    }
    
    mono[, source := "monohemispherical"]
    mono[, size_fraction := size_fraction]
    out[["monohemispherical"]] <- mono
  }
  
  if (length(out) == 0L) {
    warning(
      "No classifier output files found for size fraction: ",
      size_fraction,
      "\nDirectory checked: ",
      temp_dir
    )
    return(NULL)
  }
  
  rbindlist(out, use.names = TRUE, fill = TRUE)
}

## ============================================================
## 3. Generic reconciliation across size fractions
## ============================================================

collapse_size_fraction_tags <- function(dt,
                                        source_name,
                                        discordant_label = "discordant",
                                        map_polar_to_bipolar = FALSE,
                                        reconciliation_rule = "min_two_support") {
  dt <- as.data.table(copy(dt))
  
  if (nrow(dt) == 0L) {
    return(data.table())
  }
  
  required_cols <- c("taxon_id", "tag", "size_fraction")
  missing_cols <- setdiff(required_cols, names(dt))
  
  if (length(missing_cols) > 0L) {
    stop(
      "Input to collapse_size_fraction_tags() is missing column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  
  dt <- dt[!is.na(tag) & tag != ""]
  
  if (nrow(dt) == 0L) {
    return(data.table())
  }
  
  ## Source-level raw tag sets for inspection.
  source_sets <- dt[
    ,
    .(
      tag_set = collapse_pipe(tag)
    ),
    by = .(taxon_id, size_fraction)
  ]
  
  ## Split each pipe-separated tag into individual labels.
  long <- dt[
    ,
    .(
      tag_single = unlist(strsplit(as.character(tag), "\\|"))
    ),
    by = .(taxon_id, size_fraction)
  ]
  
  long <- long[!is.na(tag_single) & tag_single != ""]
  
  ## Strip exclusive/preferential prefix.
  long[, geo := strip_excl_prefix(tag_single)]
  
  if (map_polar_to_bipolar) {
    long[geo == "polar", geo := "bipolar"]
  }
  
  long <- long[!is.na(geo) & geo != ""]
  
  ## Build source-level geo sets for debugging/display.
  geo_source_sets <- long[
    ,
    .(
      geo_set = {
        g <- sort(unique(geo[geo != "not_categorized"]))
        if (length(g) == 0L) "not_categorized" else paste(g, collapse = "|")
      }
    ),
    by = .(taxon_id, size_fraction)
  ]
  
  ## Wide raw tag sets.
  wide_raw <- dcast(
    source_sets,
    taxon_id ~ size_fraction,
    value.var = "tag_set",
    fun.aggregate = function(x) unique(x)[1]
  )
  
  ## Wide stripped geography sets.
  wide_geo <- dcast(
    geo_source_sets,
    taxon_id ~ size_fraction,
    value.var = "geo_set",
    fun.aggregate = function(x) unique(x)[1]
  )
  
  ## Ensure every expected size fraction exists as a column.
  for (sf in size_fraction_names) {
    if (!sf %in% names(wide_raw)) {
      wide_raw[, (sf) := "not_categorized"]
    }
    if (!sf %in% names(wide_geo)) {
      wide_geo[, (sf) := "not_categorized"]
    }
  }
  
  ## Rename wide columns to stable names.
  raw_rename_old <- intersect(size_fraction_names, names(wide_raw))
  geo_rename_old <- intersect(size_fraction_names, names(wide_geo))
  
  raw_rename_new <- paste0(clean_fraction_name(raw_rename_old), "_tag_set")
  geo_rename_new <- paste0(clean_fraction_name(geo_rename_old), "_geo_set")
  
  setnames(wide_raw, raw_rename_old, raw_rename_new)
  setnames(wide_geo, geo_rename_old, geo_rename_new)
  
  wide <- merge(
    wide_geo,
    wide_raw,
    by = "taxon_id",
    all = TRUE,
    sort = FALSE
  )
  
  wide_cols <- setdiff(names(wide), "taxon_id")
  for (col in wide_cols) {
    set(
      wide,
      i = which(is.na(wide[[col]]) | wide[[col]] == ""),
      j = col,
      value = "not_categorized"
    )
  }
  
  ## Voting happens per individual geography label.
  geo_long <- unique(
    long[
      geo != "not_categorized",
      .(taxon_id, size_fraction, geo)
    ]
  )
  
  if (nrow(geo_long) == 0L) {
    final_sets <- data.table(
      taxon_id = unique(dt$taxon_id),
      final_geo = "not_categorized",
      final_geo_support = NA_character_
    )
    
    out <- merge(
      wide,
      final_sets,
      by = "taxon_id",
      all.x = TRUE,
      sort = FALSE
    )
    
    out[, source := source_name]
    setkey(out, taxon_id)
    return(out[])
  }
  
  n_real_sources <- geo_long[
    ,
    .(
      n_real_sources = uniqueN(size_fraction)
    ),
    by = taxon_id
  ]
  
  geo_support <- geo_long[
    ,
    .(
      n_sources = uniqueN(size_fraction),
      sources = paste(sort(unique(size_fraction)), collapse = "|")
    ),
    by = .(taxon_id, geo)
  ]
  
  geo_support <- n_real_sources[
    geo_support,
    on = "taxon_id"
  ]
  
  geo_support <- apply_reconciliation_rule(
    geo_support,
    reconciliation_rule = reconciliation_rule
  )
  
  final_sets <- geo_support[
    keep == TRUE,
    .(
      final_geo = paste(sort(unique(geo)), collapse = "|"),
      final_geo_support = paste(
        paste0(geo, ":", n_sources, "/", n_real_sources),
        collapse = ";"
      )
    ),
    by = taxon_id
  ]
  
  ## Any taxon with no majority-supported geography becomes discordant
  ## if it had at least one real label; otherwise not_categorized.
  all_ids <- unique(dt$taxon_id)
  no_final_ids <- setdiff(all_ids, final_sets$taxon_id)
  
  if (length(no_final_ids) > 0L) {
    ids_with_any_real <- unique(geo_long$taxon_id)
    
    add <- data.table(
      taxon_id = no_final_ids,
      final_geo = fifelse(
        no_final_ids %in% ids_with_any_real,
        discordant_label,
        "not_categorized"
      ),
      final_geo_support = NA_character_
    )
    
    final_sets <- rbindlist(
      list(final_sets, add),
      use.names = TRUE,
      fill = TRUE
    )
  }
  
  out <- merge(
    wide,
    final_sets,
    by = "taxon_id",
    all.x = TRUE,
    sort = FALSE
  )
  
  out[, source := source_name]
  setkey(out, taxon_id)
  out[]
}

## ============================================================
## 4. Recompute bihemispherical exclusivity
## ============================================================

inside_bihem_band <- function(final_geo, lat) {
  a <- abs(lat)
  
  switch(
    final_geo,
    bipolar       = a >= 60,
    polar         = a >= 60,
    subpolar      = a >= 45 & a < 60,
    subtropical   = a >= 15 & a < 45,
    high_latitude = a >= 45,
    tropical      = a < 15,
    rep(FALSE, length(lat))
  )
}

recompute_bihem_exclusivity_from_final_geo <- function(final_tbl,
                                                       observations_dt,
                                                       geo_col = "final_geo",
                                                       id_col = "taxon_id",
                                                       lat_col = "latitude",
                                                       cov_col = "norm_coverage",
                                                       discordant_vals = c(
                                                         "discordant",
                                                         "discordant_classification",
                                                         "not_categorized"
                                                       ),
                                                       add_final_tag = TRUE) {
  final_tbl <- as.data.table(copy(final_tbl))
  obs_dt <- as.data.table(observations_dt)
  
  obs <- obs_dt[
    get(cov_col) > 0 & is.finite(get(lat_col)),
    .(
      taxon_id = get(id_col),
      lat = get(lat_col)
    )
  ]
  
  geo_long <- final_tbl[
    !(get(geo_col) %in% discordant_vals) &
      !is.na(get(geo_col)),
    .(
      final_geo_single = unlist(strsplit(get(geo_col), "\\|"))
    ),
    by = taxon_id
  ]
  
  geo_long <- geo_long[
    !final_geo_single %in% discordant_vals &
      !is.na(final_geo_single) &
      final_geo_single != ""
  ]
  
  if (nrow(geo_long) == 0L) {
    final_tbl[, excl_status := "unknown"]
    if (add_final_tag) {
      final_tbl[, final_tag := get(geo_col)]
    }
    return(final_tbl[])
  }
  
  det <- obs[
    geo_long,
    on = "taxon_id",
    allow.cartesian = TRUE,
    nomatch = 0
  ]
  
  if (nrow(det) == 0L) {
    final_tbl[, excl_status := "unknown"]
    if (add_final_tag) {
      final_tbl[, final_tag := get(geo_col)]
    }
    return(final_tbl[])
  }
  
  det[
    ,
    inside := inside_bihem_band(final_geo_single[1], lat),
    by = .(taxon_id, final_geo_single)
  ]
  
  excl_single <- det[
    ,
    .(
      excl_status_single = fifelse(all(inside), "exclusive", "preferential")
    ),
    by = .(taxon_id, final_geo_single)
  ]
  
  excl_single[
    ,
    final_tag_single := paste0(excl_status_single, "_", final_geo_single)
  ]
  
  excl_collapsed <- excl_single[
    ,
    .(
      excl_status = paste(sort(unique(excl_status_single)), collapse = "|"),
      final_tag = paste(sort(unique(final_tag_single)), collapse = "|")
    ),
    by = taxon_id
  ]
  
  final_tbl[, excl_status := "unknown"]
  final_tbl[, final_tag := get(geo_col)]
  
  final_tbl[
    excl_collapsed,
    `:=`(
      excl_status = i.excl_status,
      final_tag = i.final_tag
    ),
    on = "taxon_id"
  ]
  
  final_tbl[]
}

## ============================================================
## 5. Recompute monohemispherical exclusivity
## ============================================================

parse_mono_geo <- function(tag) {
  hemi <- NA_character_
  base <- tag
  
  if (grepl("_(north|south)$", tag)) {
    hemi <- sub(".*_(north|south)$", "\\1", tag)
    base <- sub("_(north|south)$", "", tag)
  }
  
  list(geo = base, hemi = hemi)
}

inside_mono_band <- function(final_geo, lat) {
  p <- parse_mono_geo(final_geo)
  a <- abs(lat)
  
  hemi_ok <- rep(TRUE, length(lat))
  
  if (!is.na(p$hemi)) {
    hemi_ok <- if (p$hemi == "north") {
      lat >= 0
    } else {
      lat < 0
    }
  }
  
  band_ok <- switch(
    p$geo,
    tropical      = a < 15,
    subtropical   = a >= 15 & a < 45,
    subpolar      = a >= 45 & a < 60,
    polar         = a >= 60,
    high_latitude = a >= 45,
    rep(FALSE, length(lat))
  )
  
  hemi_ok & band_ok
}

recompute_mono_exclusivity_from_final_geo <- function(final_tbl,
                                                      observations_dt,
                                                      geo_col = "final_geo",
                                                      id_col = "taxon_id",
                                                      lat_col = "latitude",
                                                      cov_col = "norm_coverage",
                                                      discordant_vals = c(
                                                        "discordant",
                                                        "discordant_classification",
                                                        "not_categorized"
                                                      ),
                                                      add_final_tag = TRUE) {
  final_tbl <- as.data.table(copy(final_tbl))
  obs_dt <- as.data.table(observations_dt)
  
  obs <- obs_dt[
    get(cov_col) > 0 & is.finite(get(lat_col)),
    .(
      taxon_id = get(id_col),
      lat = get(lat_col)
    )
  ]
  
  geo_long <- final_tbl[
    !(get(geo_col) %in% discordant_vals) &
      !is.na(get(geo_col)),
    .(
      final_geo_single = unlist(strsplit(get(geo_col), "\\|"))
    ),
    by = taxon_id
  ]
  
  geo_long <- geo_long[
    !final_geo_single %in% discordant_vals &
      !is.na(final_geo_single) &
      final_geo_single != ""
  ]
  
  if (nrow(geo_long) == 0L) {
    final_tbl[, excl_status := "unknown"]
    if (add_final_tag) {
      final_tbl[, final_tag := get(geo_col)]
    }
    return(final_tbl[])
  }
  
  det <- obs[
    geo_long,
    on = "taxon_id",
    allow.cartesian = TRUE,
    nomatch = 0
  ]
  
  if (nrow(det) == 0L) {
    final_tbl[, excl_status := "unknown"]
    if (add_final_tag) {
      final_tbl[, final_tag := get(geo_col)]
    }
    return(final_tbl[])
  }
  
  det[
    ,
    inside := inside_mono_band(final_geo_single[1], lat),
    by = .(taxon_id, final_geo_single)
  ]
  
  excl_single <- det[
    ,
    .(
      excl_status_single = fifelse(all(inside), "exclusive", "preferential")
    ),
    by = .(taxon_id, final_geo_single)
  ]
  
  excl_single[
    ,
    final_tag_single := paste0(excl_status_single, "_", final_geo_single)
  ]
  
  excl_collapsed <- excl_single[
    ,
    .(
      excl_status = paste(sort(unique(excl_status_single)), collapse = "|"),
      final_tag = paste(sort(unique(final_tag_single)), collapse = "|")
    ),
    by = taxon_id
  ]
  
  final_tbl[, excl_status := "unknown"]
  final_tbl[, final_tag := get(geo_col)]
  
  final_tbl[
    excl_collapsed,
    `:=`(
      excl_status = i.excl_status,
      final_tag = i.final_tag
    ),
    on = "taxon_id"
  ]
  
  final_tbl[]
}

## ============================================================
## 6. Read all classifier outputs
## ============================================================

message("Reading size-fraction classifier outputs...")

all_classifications <- rbindlist(
  Map(
    read_one_fraction,
    classification_dirs,
    size_fraction_names
  ),
  use.names = TRUE,
  fill = TRUE
)

if (nrow(all_classifications) == 0L) {
  stop("No classification outputs found to reconcile.", call. = FALSE)
}

required_class_cols <- c("taxon_id", "tag", "source", "size_fraction")
missing_class_cols <- setdiff(required_class_cols, names(all_classifications))

if (length(missing_class_cols) > 0L) {
  stop(
    "Combined classifier output is missing required column(s): ",
    paste(missing_class_cols, collapse = ", "),
    call. = FALSE
  )
}

## ============================================================
## 7. Reconcile bihem and mono separately
## ============================================================

message("Reconciling bihemispherical classifications...")

bihem_input <- all_classifications[
  source == "bihemispherical"
]

if (nrow(bihem_input) > 0L) {
  bihem_reconciled <- collapse_size_fraction_tags(
    bihem_input,
    source_name = "bihemispherical",
    discordant_label = discordant_label,
    map_polar_to_bipolar = map_bihem_polar_to_bipolar,
    reconciliation_rule = reconciliation_rule
  )
  
  bihem_reconciled <- recompute_bihem_exclusivity_from_final_geo(
    bihem_reconciled,
    observations_dt = obs_input,
    geo_col = "final_geo",
    id_col = "taxon_id",
    lat_col = "latitude",
    cov_col = "norm_coverage"
  )
} else {
  bihem_reconciled <- data.table()
  warning("No bihemispherical classifications found.")
}

message("Reconciling monohemispherical classifications...")

mono_input <- all_classifications[
  source == "monohemispherical"
]

if (nrow(mono_input) > 0L) {
  mono_reconciled <- collapse_size_fraction_tags(
    mono_input,
    source_name = "monohemispherical",
    discordant_label = discordant_label,
    map_polar_to_bipolar = FALSE,
    reconciliation_rule = reconciliation_rule
  )
  
  mono_reconciled <- recompute_mono_exclusivity_from_final_geo(
    mono_reconciled,
    observations_dt = obs_input,
    geo_col = "final_geo",
    id_col = "taxon_id",
    lat_col = "latitude",
    cov_col = "norm_coverage"
  )
} else {
  mono_reconciled <- data.table()
  warning("No monohemispherical classifications found.")
}

## ============================================================
## 8. Combined reconciled table
## ============================================================

## Bihemispherical classifications take priority in the final
## non-overlapping combined output.
##
## Rationale:
##   The pipeline classifies bihemispherical taxa first.
##   Monohemispherical classification is only meant to resolve taxa
##   that remain unresolved after the bihemispherical stage.
##
## We still save the full mono_reconciled table separately below,
## so no monohemispherical reconciliation information is lost.

## ============================================================
## 8. Combined reconciled table
## ============================================================

resolved_bihem_ids <- character()

if (nrow(bihem_reconciled) > 0L) {
  resolved_bihem_ids <- bihem_reconciled[
    !is.na(final_geo) &
      !final_geo %in% c(
        "not_categorized",
        "discordant",
        "discordant_classification"
      ),
    unique(taxon_id)
  ]
}

mono_reconciled_final <- mono_reconciled

if (length(resolved_bihem_ids) > 0L && nrow(mono_reconciled_final) > 0L) {
  mono_reconciled_final <- mono_reconciled_final[
    !taxon_id %in% resolved_bihem_ids
  ]
}

combined_reconciled <- rbindlist(
  list(
    bihem_reconciled,
    mono_reconciled_final
  ),
  use.names = TRUE,
  fill = TRUE
)

combined_reconciled <- unique(combined_reconciled)

## Optional split table: one row per taxon_id x final_tag component.
split_final_tags <- function(dt) {
  dt <- as.data.table(copy(dt))
  
  if (nrow(dt) == 0L) {
    return(dt)
  }
  
  dt[, .row_id_tmp := .I]
  
  split_dt <- dt[
    ,
    {
      fg0 <- as.character(.SD[["final_geo"]][1])
      ft0 <- as.character(.SD[["final_tag"]][1])
      es0 <- as.character(.SD[["excl_status"]][1])
      
      fg <- if (is.na(fg0) || fg0 == "") {
        character()
      } else {
        unlist(strsplit(fg0, "\\|"))
      }
      
      ft <- if (is.na(ft0) || ft0 == "") {
        character()
      } else {
        unlist(strsplit(ft0, "\\|"))
      }
      
      es <- if (is.na(es0) || es0 == "") {
        character()
      } else {
        unlist(strsplit(es0, "\\|"))
      }
      
      if (length(ft) > 0L && any(grepl("^(exclusive|preferential)_", ft))) {
        es_from_tag <- fifelse(
          grepl("^exclusive_", ft),
          "exclusive",
          fifelse(grepl("^preferential_", ft), "preferential", NA_character_)
        )
        
        fg_from_tag <- sub("^(exclusive_|preferential_)", "", ft)
        
        data.table(
          final_geo = fg_from_tag,
          excl_status = es_from_tag,
          final_tag = ft
        )
        
      } else if (length(fg) > 0L) {
        data.table(
          final_geo = fg,
          excl_status = if (length(es) == length(fg)) es else NA_character_,
          final_tag = if (length(ft) == length(fg)) ft else fg
        )
        
      } else {
        data.table(
          final_geo = NA_character_,
          excl_status = NA_character_,
          final_tag = NA_character_
        )
      }
    },
    by = .row_id_tmp,
    .SDcols = names(dt)
  ]
  
  keep_cols <- setdiff(names(dt), c("final_geo", "excl_status", "final_tag"))
  
  out <- merge(
    dt[, ..keep_cols],
    split_dt,
    by = ".row_id_tmp",
    allow.cartesian = TRUE,
    sort = FALSE
  )
  
  out[, .row_id_tmp := NULL]
  unique(out[])
}

combined_reconciled_split <- split_final_tags(combined_reconciled)

## ============================================================
## 9. Save outputs
## ============================================================

reconciled_dir <- file.path(output_dir, "reconciled")
dir.create(reconciled_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(
  bihem_reconciled,
  file.path(reconciled_dir, "reconciled_bihemispherical.rds")
)

fwrite(
  bihem_reconciled,
  file.path(reconciled_dir, "reconciled_bihemispherical.tsv"),
  sep = "\t"
)

saveRDS(
  mono_reconciled,
  file.path(reconciled_dir, "reconciled_monohemispherical.rds")
)

fwrite(
  mono_reconciled,
  file.path(reconciled_dir, "reconciled_monohemispherical.tsv"),
  sep = "\t"
)

saveRDS(
  combined_reconciled,
  file.path(reconciled_dir, "reconciled_biogeography_combined.rds")
)

fwrite(
  combined_reconciled,
  file.path(reconciled_dir, "reconciled_biogeography_combined.tsv"),
  sep = "\t"
)

saveRDS(
  combined_reconciled_split,
  file.path(reconciled_dir, "reconciled_biogeography_combined_split.rds")
)

fwrite(
  combined_reconciled_split,
  file.path(reconciled_dir, "reconciled_biogeography_combined_split.tsv"),
  sep = "\t"
)

fwrite(
  all_classifications,
  file.path(reconciled_dir, "all_size_fraction_classifications.tsv"),
  sep = "\t"
)

message("Saved reconciled outputs to: ", reconciled_dir)

if (nrow(combined_reconciled) > 0L) {
  message("Final tag counts:")
  print(
    combined_reconciled[
      ,
      .N,
      by = .(source, final_tag)
    ][order(source, -N)]
  )
}