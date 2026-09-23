## ============================================================
## classify_monohemispherical.R
## Classifies remaining taxa as monohemispherical / unimodal
##
## Expected objects from previous scripts:
##   peaks_merged_global
##   observations_input
##   unclassified_ids
##   bihem_dt
##   output_dir
##   temp_dir
## ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(readr)
  library(fst)
})

## ============================================================
## 0. Set latzone cutoffs
## ============================================================

if (!exists("latzone_cutoffs")) {
  latzone_cutoffs <- list(
    tropical = 15,
    subpolar = 45,
    polar = 60
  )
}

tropical_cutoff <- latzone_cutoffs$tropical
subpolar_cutoff <- latzone_cutoffs$subpolar
polar_cutoff <- latzone_cutoffs$polar

## ============================================================
## 1. Check required objects
## ============================================================

required_objects <- c(
  "peaks_merged_global",
  "observations_input",
  "unclassified_ids",
  "bihem_dt",
  "output_dir",
  "temp_dir"
)

missing_objects <- required_objects[!vapply(required_objects, exists, logical(1))]

if (length(missing_objects) > 0) {
  stop(
    "classify_unimodal.R is missing required objects: ",
    paste(missing_objects, collapse = ", "),
    call. = FALSE
  )
}

if (!dir.exists(output_dir)) {
  stop("output_dir does not exist: ", output_dir, call. = FALSE)
}

if (!dir.exists(temp_dir)) {
  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
}

## ============================================================
## 2. Validate input
## ============================================================

required_peak_cols <- c("taxon_id", "latitude", "peak_fit")
missing_peak_cols <- setdiff(required_peak_cols, names(peaks_merged_global))

if (length(missing_peak_cols) > 0) {
  stop(
    "peaks_merged_global is missing required columns: ",
    paste(missing_peak_cols, collapse = ", "),
    call. = FALSE
  )
}

required_obs_cols <- c("taxon_id", "latitude", "norm_coverage")
missing_obs_cols <- setdiff(required_obs_cols, names(observations_input))

if (length(missing_obs_cols) > 0) {
  stop(
    "observations_input is missing required columns: ",
    paste(missing_obs_cols, collapse = ", "),
    call. = FALSE
  )
}

required_bihem_cols <- c("taxon_id", "final_tag")
missing_bihem_cols <- setdiff(required_bihem_cols, names(bihem_dt))

if (length(missing_bihem_cols) > 0) {
  stop(
    "bihem_dt is missing required columns: ",
    paste(missing_bihem_cols, collapse = ", "),
    call. = FALSE
  )
}

## Validate unclassified_ids
if (!exists("unclassified_ids")) {
  stop(
    "Missing required object: unclassified_ids. ",
    "Run classify_bimodal.R before classify_unimodal.R.",
    call. = FALSE
  )
}

if (!is.atomic(unclassified_ids)) {
  stop(
    "unclassified_ids must be an atomic vector of taxon_id values.",
    call. = FALSE
  )
}

unclassified_ids <- as.character(unclassified_ids)
unclassified_ids <- unique(unclassified_ids)

if (length(unclassified_ids) == 0L) {
  warning(
    "unclassified_ids is empty. No taxa remain for unimodal classification."
  )
}

if (anyNA(unclassified_ids)) {
  stop(
    "unclassified_ids contains NA values.",
    call. = FALSE
  )
}

if (any(unclassified_ids == "")) {
  stop(
    "unclassified_ids contains empty string values.",
    call. = FALSE
  )
}

missing_unclassified_ids <- setdiff(
  unclassified_ids,
  unique(as.character(peaks_merged_global$taxon_id))
)

if (length(missing_unclassified_ids) > 0L) {
  stop(
    "Some unclassified_ids are not present in peaks_merged_global$taxon_id. ",
    "Example missing IDs: ",
    paste(head(missing_unclassified_ids, 10), collapse = ", "),
    call. = FALSE
  )
}

## ============================================================
## 3. Prepare input
## ============================================================

library(data.table)

peaks_merged_global <- as.data.table(peaks_merged_global)
observations_input  <- as.data.table(observations_input)
bihem_dt            <- as.data.table(bihem_dt)

## ============================================================
## 4. Load functions
## ============================================================

## --------------------------------------------------------------------- ##
##  make_mono_tags()                                                     ##
##  • Classifies each taxon based on the latitude of its detected peaks. ##
##  • Assigns every peak to one fine latitude zone using these bounds:   ##
##        tropical:          |lat| < 15                                  ##
##        subtropical:       15 ≤ |lat| < 45                             ##
##        subpolar:          45 ≤ |lat| < 60                             ##
##        polar:             |lat| ≥ 60                                  ##
##    with separate north/south labels outside the tropics.              ##
##  • Then summarizes the peak zones per taxon into one broad category:  ##
##        tropical, polar N or S, subpolar N or S, subtropical N or S,   ##
##        high latitude N or S, or unclassified.                         ##
##  • Classification is based on peak locations only.                    ##
##  • Exclusivity/preferential status is added later using observed      ##
##    detections in compute_mono_exclusivity().                          ##
##  • Returns a data.table with one row per taxon_id containing:         ##
##        category          → broad latitude classification              ##
## --------------------------------------------------------------------- ##

make_mono_tags <- function(peaks_dt,
                           unclassified_ids,
                           coverage_dt,
                           id_col  = "taxon_id",
                           lat_col = "latitude") {
  
  library(data.table)
  
  peaks <- as.data.table(peaks_dt)
  unclassified_ids <- as.character(unclassified_ids)
  
  zone_of <- function(lat) {
  a <- abs(lat)
  fcase(
    a < tropical_cutoff, "tropical",
    a < subpolar_cutoff, "subtropical",
    a < polar_cutoff, "subpolar",
    default = "polar"
  )
}
  
  ## 1. Keep only still-unclassified IDs
  peaks <- peaks[get(id_col) %in% unclassified_ids]
  
  if (nrow(peaks) == 0L) {
    return(data.table(
      taxon_id = character(),
      mono_tag = character()
    ))
  }
  
  ## 2. Valid monohemispherical IDs
  ## Tropical exception:
  ##   all peaks tropical, can cross equator.
  ##
  ## Non-tropical:
  ##   all peaks must be in one hemisphere,
  ##   and zones must be either:
  ##     - one zone only
  ##     - polar + subpolar in same hemisphere
  valid_dt <- peaks[
    ,
    {
      lats <- get(lat_col)
      a <- abs(lats)
      
      is_tropical <- all(a < 15)
      
      if (is_tropical) {
        .(valid = TRUE)
      } else {
        hemi_same <- uniqueN(sign(lats[lats != 0])) <= 1L
        zones <- unique(zone_of(lats))
        ok_zone <- length(zones) == 1L || setequal(zones, c("polar", "subpolar"))
        
        .(valid = hemi_same && ok_zone)
      }
    },
    by = id_col
  ]
  
  valid_ids <- valid_dt[valid == TRUE, get(id_col)]
  
  peaks <- peaks[get(id_col) %in% valid_ids]
  
  if (nrow(peaks) == 0L) {
    return(data.table(
      taxon_id = character(),
      mono_tag = character()
    ))
  }
  
  ## 3. Tallest peak per genome determines the specific mono category
  tallest <- peaks[
    ,
    .SD[which.max(peak_fit)],
    by = id_col
  ]
  
  tallest[
    ,
    `:=`(
      tallest_lat = get(lat_col),
      specific_zone = zone_of(get(lat_col))
    )
  ]
  
  ## 4. Assign hemisphere only for non-tropical tags
  tallest[
    ,
    hemi := fifelse(
      specific_zone == "tropical",
      NA_character_,
      fifelse(tallest_lat >= 0, "north", "south")
    )
  ]
  
  ## 5. Specific mono tag from tallest peak
  tallest[
    ,
    specific_tag := fifelse(
      specific_zone == "tropical",
      "tropical",
      paste0(specific_zone, "_", hemi)
    )
  ]
  
  ## 6. Add parent high_latitude tag ONLY when tallest peak is polar/subpolar
  tallest[
    ,
    parent_tag := fifelse(
      specific_zone %in% c("polar", "subpolar"),
      paste0("high_latitude_", hemi),
      NA_character_
    )
  ]
  
  ## 7. Final source-level mono tag set
  ## Never produces high_latitude_north/south alone.
  tallest[
    ,
    mono_tag := mapply(
      function(specific, parent) {
        tags <- c(specific, parent)
        tags <- tags[!is.na(tags)]
        paste(tags, collapse = "|")
      },
      specific_tag,
      parent_tag
    )
  ]
  
  out <- tallest[
    ,
    .(
      taxon_id = get(id_col),
      mono_tag
    )
  ]
  
  unique(out)
}

## --------------------------------------------------------------------- ##
##  compute_mono_exclusivity()                                           ##
##  • A taxon receives one exclusivity / preferential tag based on       ##
##    observations_input.                                                ##
##  • Windows can overlap: polar and subpolar are both subsets of        ##
##    high-latitude.                                                     ##
##  • Returns:                                                           ##
##        $final_tbl   → classification + exclusive / preferential tag   ##
## --------------------------------------------------------------------- ##

compute_mono_exclusivity <- function(final_tbl,
                                                      coverage_dt,
                                                      geo_col = "final_geo",
                                                      id_col = "virus_genome",
                                                      lat_col = "manual_latitude",
                                                      cov_col = "mean_norm_coverage",
                                                      discordant_vals = c(
                                                        "discordant",
                                                        "discordant_classification",
                                                        "not_categorized"
                                                      ),
                                                      add_final_tag = TRUE) {
  
  library(data.table)
  
  final_tbl <- as.data.table(copy(final_tbl))
  cov <- as.data.table(coverage_dt)
  
  ## standardize internally so this works for virus_genome or taxon_id
  final_tbl[, id_tmp := get(id_col)]
  
  parse_geo <- function(tag) {
    hemi <- NA_character_
    base <- tag
    
    if (grepl("_(north|south)$", tag)) {
      hemi <- sub(".*_(north|south)$", "\\1", tag)
      base <- sub("_(north|south)$", "", tag)
    }
    
    list(geo = base, hemi = hemi)
  }
  
  inside_band <- function(final_geo, lat) {
    p <- parse_geo(final_geo)
    a <- abs(lat)
    
    hemi_ok <- rep(TRUE, length(lat))
    
    if (!is.na(p$hemi)) {
      hemi_ok <- if (p$hemi == "north") lat >= 0 else lat < 0
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
  
  obs <- cov[
    get(cov_col) > 0 & is.finite(get(lat_col)),
    .(
      id_tmp = get(id_col),
      lat = get(lat_col)
    )
  ]
  
  geo_long <- final_tbl[
    !(get(geo_col) %in% discordant_vals) &
      !is.na(get(geo_col)),
    .(
      final_geo_single = unlist(strsplit(get(geo_col), "\\|"))
    ),
    by = id_tmp
  ]
  
  geo_long <- geo_long[
    !final_geo_single %in% discordant_vals &
      !is.na(final_geo_single) &
      final_geo_single != ""
  ]
  
  if (nrow(geo_long) == 0L) {
    final_tbl[, excl_status := "unknown"]
    if (add_final_tag) final_tbl[, final_tag := get(geo_col)]
    final_tbl[, id_tmp := NULL]
    return(final_tbl[])
  }
  
  det <- obs[
    geo_long,
    on = "id_tmp",
    allow.cartesian = TRUE,
    nomatch = 0
  ]
  
  det[
    ,
    inside := inside_band(final_geo_single[1], lat),
    by = .(id_tmp, final_geo_single)
  ]
  
  excl_single <- det[
    ,
    .(
      excl_status_single = fifelse(all(inside), "exclusive", "preferential")
    ),
    by = .(id_tmp, final_geo_single)
  ]
  
  excl_single[
    ,
    final_tag_single := paste0(excl_status_single, "_", final_geo_single)
  ]
  
  excl_collapsed <- excl_single[
    ,
    .(
      excl_status = paste(excl_status_single, collapse = "|"),
      final_tag = paste(final_tag_single, collapse = "|")
    ),
    by = id_tmp
  ]
  
  final_tbl[, excl_status := "unknown"]
  final_tbl[, final_tag := get(geo_col)]
  
  final_tbl[
    excl_collapsed,
    `:=`(
      excl_status = i.excl_status,
      final_tag = i.final_tag
    ),
    on = "id_tmp"
  ]
  
  final_tbl[, id_tmp := NULL]
  final_tbl[]
}

## ============================================================
## 4. Execute functions
## ============================================================

# Build monohemispherical–tag table

mono_dt <- make_mono_tags(
  peaks_dt = peaks_merged_global,
  unclassified_ids = unclassified_ids,
  coverage_dt      = observations_input,
  id_col           = "taxon_id"
)

#mono_dt[1:5]

# Add exclusivity tag

mono_dt2 <- compute_mono_exclusivity(
  final_tbl   = mono_dt,
  coverage_dt = observations_input,
  geo_col     = "mono_tag",
  id_col      = "taxon_id",
  lat_col     = "latitude",
  cov_col     = "norm_coverage"
)

monohem_dt <- copy(mono_dt2)
setnames(monohem_dt, "mono_tag", "final_geo")
monohem_dt[, tag := final_tag]
monohem_dt[, source := "monohemispherical"]

monohem_dt[, .N, by = final_tag][order(-N)] # preview

## ============================================================
## 4. Prepare final output
## ============================================================

# Still-unclassified taxa
still_unclassified     <- setdiff(unclassified_ids,
                                      mono_dt2$taxon_id)

## quick counts
#length(still_unclassified)

# Combine bihemispherical and monohemispherical taxa
combined_dt <- rbindlist(
  list(
    monohem_dt[, .(taxon_id, final_geo, excl_status, final_tag, tag, source)],
    bihem_dt[, .(taxon_id, final_tag, tag, source)]
  ),
  use.names = TRUE,
  fill = TRUE
)

combined_dt <- combined_dt[
  ,
  .(taxon_id, final_geo, excl_status, final_tag, source)
]

## ============================================================
## 5. Save temporary and final outputs
## ============================================================

monohem_dt_file <- file.path(temp_dir, "monohemispherical_vOTUs.rds")
still_unclassified_file <- file.path(output_dir, "unclassified_taxid.txt")
combined_out_file <- file.path(output_dir, "biogeog_out.tsv")

saveRDS(monohem_dt, monohem_dt_file)
writeLines(still_unclassified, still_unclassified_file)
write_tsv(combined_dt, combined_out_file)

message("Saved monohemispherical vOTUs to: ", monohem_dt_file)
message("Saved still-unclassified taxon IDs to: ", still_unclassified_file)
message("Saved final biogeography output to: ", combined_out_file)
