## ============================================================
## predict_peaks.R
## Finds latitude peaks from prediction and observation inputs
##
## Expected objects from master script:
##   predictions_input
##   observations_input
##   output_dir
##   temp_dir
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(fst)
})

## ============================================================
## 0. Config
## ============================================================

default_config <- list(
  min_support    = 2,
  min_rel_height = 0.5,
  min_abs_fit    = 0.005,
  min_sep_deg    = 5,
  tropical_cutoff  = 15,
  subpolar_cutoff  = 45,
  polar_cutoff     = 60
)

config <- default_config

if (exists("config_file") && !is.null(config_file) && file.exists(config_file)) {
  config_tbl <- readr::read_tsv(config_file, show_col_types = FALSE)

  required_config_cols <- c("parameter", "value")
  missing_config_cols <- setdiff(required_config_cols, names(config_tbl))

  if (length(missing_config_cols) > 0) {
    stop(
      "Config file is missing required columns: ",
      paste(missing_config_cols, collapse = ", "),
      call. = FALSE
    )
  }

  config_from_file <- as.list(config_tbl$value)
  names(config_from_file) <- config_tbl$parameter

  config_from_file <- lapply(config_from_file, as.numeric)

  config[names(config_from_file)] <- config_from_file
}

# Validate latitude cutoffs after defaults/config are finalized
if (!(config$tropical_cutoff < config$subpolar_cutoff &&
      config$subpolar_cutoff < config$polar_cutoff)) {
  stop(
    "Latitudinal cutoffs must satisfy: ",
    "tropical_cutoff < subpolar_cutoff < polar_cutoff.",
    call. = FALSE
  )
}

latzone_cutoffs <- list(
  tropical = config$tropical_cutoff,
  subpolar = config$subpolar_cutoff,
  polar = config$polar_cutoff
)

message("Using config:")
print(config)


## ============================================================
## 1. Check required objects from master script
## ============================================================

required_objects <- c(
  "predictions_input",
  "observations_input",
  "output_dir",
  "temp_dir"
)

missing_objects <- required_objects[!vapply(required_objects, exists, logical(1))]

if (length(missing_objects) > 0) {
  stop(
    "predict_peaks.R is missing required objects from the master script: ",
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
## 2. Validate input columns
## ============================================================

required_pred_cols <- c("taxon_id", "latitude", "fit")
missing_pred_cols <- setdiff(required_pred_cols, names(predictions_input))

if (length(missing_pred_cols) > 0) {
  stop(
    "predictions_input is missing required columns: ",
    paste(missing_pred_cols, collapse = ", "),
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

## ============================================================
## 3. Prepare prediction input
## ============================================================

pred <- predictions_input %>%
  transmute(
    taxon_id,
    latitude,
    fit = pmax(fit, 0)
  ) %>%
  arrange(taxon_id, latitude) %>%
  filter(is.finite(latitude), is.finite(fit))

setDT(pred)
pred[, lat_bin := round(latitude)]


## ============================================================
## 4. Build support table from merged coverage table
## ============================================================

# Counts how many samples have norm_coverage > 0 per vOTU × 1° bin
support_by_lat <- observations_input %>%
  mutate(lat_bin = round(latitude)) %>%
  group_by(taxon_id, lat_bin) %>%
  summarise(n_pos = sum(norm_coverage > 0, na.rm = TRUE), .groups = "drop") %>%
  as.data.table()


## ============================================================
## 5. Merge support into predictions
## ============================================================

pred_with_support <- merge(
  pred[, .(taxon_id, latitude, lat_bin, fit)],
  support_by_lat,
  by.x = c("taxon_id", "lat_bin"),
  by.y = c("taxon_id", "lat_bin"),
  all.x = TRUE
)

pred_with_support[is.na(n_pos), n_pos := 0L]

## ============================================================
## 6. Assign zones: north, south, equatorial
## ============================================================
equatorial_cutoff <- config$tropical_cutoff

pred_with_support[, zone := fifelse(
  latitude >= equatorial_cutoff, "north",
  fifelse(latitude <= -equatorial_cutoff, "south", "equatorial")
)]

## ==================================================================
## 7. Define peak detection function — returns *all* supported peaks
## ==================================================================

detect_all_supported_peaks <- function(df, min_support = 1) {
  # df: latitude, fit, n_pos
  setorder(df, latitude)
  
  # Keep bins with enough support
  supported <- df[n_pos >= min_support]
  
  if (nrow(supported) == 0) {
    return(data.table(
      latitude = numeric(0),
      peak_fit = numeric(0),
      n_pos = integer(0)
    ))
  }
  
  # Identify contiguous supported regions
  rle_ids <- rleid(supported$latitude - seq_len(nrow(supported)))
  supported[, region_id := rle_ids]
  
  # For each supported region, take the GAM fit maximum as the "peak"
  peaks <- supported[, .SD[which.max(fit)], by = region_id]
  peaks[, region_id := NULL]
  setnames(peaks, "fit", "peak_fit")
  
  return(peaks)
}

## ============================================================
## 8. Detect peaks for every virus × zone
## ============================================================

setorder(pred_with_support, taxon_id, zone, latitude)

peaks_clean <- pred_with_support[
  ,
  detect_all_supported_peaks(.SD, min_support = 1),
  by = .(taxon_id, zone),
  .SDcols = c("latitude", "fit", "n_pos")
]

## ============================================================
## 9. Add local support
## ============================================================

add_local_support <- function(peaks, support_by_lat, radius_deg = 5) {
  peaks[, peak_id := .I]
  
  local_support <- support_by_lat[
    peaks,
    on = .(taxon_id),
    allow.cartesian = TRUE
  ][
    abs(lat_bin - latitude) <= radius_deg,
    .(
      n_pos_local = sum(n_pos, na.rm = TRUE),
      n_bins_local = sum(n_pos > 0, na.rm = TRUE)
    ),
    by = .(peak_id)
  ]
  
  peaks <- merge(peaks, local_support, by = "peak_id", all.x = TRUE)
  peaks[is.na(n_pos_local), n_pos_local := 0L]
  peaks[is.na(n_bins_local), n_bins_local := 0L]
  peaks[, peak_id := NULL]
  
  peaks[]
}

# 9a) add local support
peaks_clean <- add_local_support(
  peaks_clean,
  support_by_lat,
  radius_deg = config$min_sep_deg
)

# 9b) keep only "reliable" peaks first (support + absolute height)
peaks_tmp <- peaks_clean[
  peak_fit >= config$min_abs_fit &
    (
      n_pos >= config$min_support |
        (
          n_pos_local >= config$min_support &
            n_bins_local >= 2
        )
    )
]

## ============================================================
## 10. Apply adaptive filtering
## ============================================================

# 10a) compute relative height only among those reliable peaks (per vOTU x zone)
peaks_tmp[, fit_ratio := peak_fit / max(peak_fit, na.rm = TRUE),
          by = .(taxon_id, zone)]

# 10b) now apply the relative-height threshold
peaks_filtered <- peaks_tmp[
  fit_ratio >= config$min_rel_height
]

## ============================================================
## 11. Merge close peaks
## ============================================================
 
# Merge close peaks ACROSS zones (per vOTU)
merge_close_peaks_global <- function(dt, min_sep_deg = 5) {
  dt[, {
    peaks <- copy(.SD)
    setorder(peaks, -peak_fit)           # tallest first
    kept <- peaks[0]
    
    for (i in seq_len(nrow(peaks))) {
      if (nrow(kept) == 0 || all(abs(peaks$latitude[i] - kept$latitude) > min_sep_deg)) {
        kept <- rbind(kept, peaks[i])
      }
    }
    
    setorder(kept, latitude)
    kept
  }, by = taxon_id]
}

peaks_merged_global <- merge_close_peaks_global(peaks_filtered, min_sep_deg = config$min_sep_deg)

## ============================================================
## 12. Summarize peak properties
## ============================================================
 
# Max peak heights per zone per vOTU
peak_summary_merged <- peaks_merged_global[, .(
  north_peak      = max(peak_fit[zone == "north"], na.rm = TRUE),
  south_peak      = max(peak_fit[zone == "south"], na.rm = TRUE),
  equatorial_peak = max(peak_fit[zone == "equatorial"], na.rm = TRUE),
  n_north         = sum(zone == "north"),
  n_south         = sum(zone == "south"),
  n_equatorial    = sum(zone == "equatorial")
), by = taxon_id]

# Replace missing maxima (no peaks in zone) with 0
cols_to_fix <- c("north_peak", "south_peak", "equatorial_peak")
for (col in cols_to_fix) {
  peak_summary_merged[!is.finite(get(col)), (col) := 0]
}

## ============================================================
## OPTIONAL:  Summarize peak properties
## ============================================================
# Count how many vOTUs remain
message("peaks_merged_global dimensions: ", paste(dim(peaks_merged_global), collapse = " x "))

message("Unique taxa before filtering: ", length(unique(peaks_clean$taxon_id)))
message("Unique taxa after filtering: ", length(unique(peaks_filtered$taxon_id)))
message("Unique taxa after merging: ", length(unique(peaks_merged_global$taxon_id)))

message(
  nrow(peaks_merged_global),
  " peaks remaining out of ",
  nrow(peaks_filtered),
  " after merging close peaks."
)
## ============================================================
## 13.  Save temp files
## ============================================================

if (!dir.exists(temp_dir)) {
  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
}

required_outputs <- c("peaks_merged_global", "peak_summary_merged")
missing_outputs <- required_outputs[!vapply(required_outputs, exists, logical(1))]

if (length(missing_outputs) > 0) {
  stop(
    "predict_peaks.R did not create expected output object(s): ",
    paste(missing_outputs, collapse = ", "),
    call. = FALSE
  )
}

peaks_merged_file <- file.path(temp_dir, "peaks_merged_global.fst")
peak_summary_file <- file.path(temp_dir, "peak_summary_merged.fst")

write_fst(peaks_merged_global, peaks_merged_file)
write_fst(peak_summary_merged, peak_summary_file)

message("Saved peaks to: ", peaks_merged_file)
message("Saved peak summary to: ", peak_summary_file)
