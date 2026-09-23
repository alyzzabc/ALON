## ============================================================
## classify_bihemispherical.R
## Classifies bihemispherical / broad latitudinal tags
##
## Expected objects from master script or predict_peaks.R:
##   peaks_merged_global
##   observations_input
##   output_dir
##   temp_dir
##
## Expected input columns:
##   peaks_merged_global:
##     taxon_id, latitude, peak_fit
##
##   observations_input:
##     taxon_id, latitude, norm_coverage
## ============================================================

suppressPackageStartupMessages({
  library(data.table)
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
  "output_dir",
  "temp_dir"
)

missing_objects <- required_objects[!vapply(required_objects, exists, logical(1))]

if (length(missing_objects) > 0) {
  stop(
    "classify_bihemispherical.R is missing required objects: ",
    paste(missing_objects, collapse = ", "),
    call. = FALSE
  )
}

## ============================================================
## 2. Validate input columns
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

## ============================================================
## 3. Prepare input
## ============================================================

peaks_merged_global <- as.data.table(peaks_merged_global)
observations_input  <- as.data.table(observations_input)

## ============================================================
## 4. Load functions
## ============================================================

## --------------------------------------------------------------------- ##
##  classify_latitude()                                                  ##
##  • Classifies each taxon based on the latitude of its detected peaks. ##
##  • Assigns every peak to one fine latitude zone using these bounds:   ##
##        tropical:          |lat| < 15                                  ##
##        subtropical:       15 ≤ |lat| < 45                             ##
##        subpolar:          45 ≤ |lat| < 60                             ##
##        polar:             |lat| ≥ 60                                  ##
##    with separate north/south labels outside the tropics.              ##
##  • Then summarizes the peak zones per taxon into one broad category:  ##
##        tropical, bipolar, bi_subpolar, bi_subtropical,                ##
##        bi_high_latitudes, or unclassified.                            ##
##  • Classification is based on peak locations only.                    ##
##  • Exclusivity/preferential status is added later using observed      ##
##    detections in tag_exclusivity_multi().                             ##
##  • Returns a data.table with 1-2 rows per taxon_id containing:        ##
##        category          → broad latitude classification              ##
##        top_peaks_summary → readable summary of peak latitudes/fits    ##
## --------------------------------------------------------------------- ##

classify_latitude <- function(peaks_dt) {
  library(data.table)
  peaks <- copy(as.data.table(peaks_dt))
  
  # Reassign coarse zone from Latitude (on *peaks*, not peaks_merged_global)
  peaks[, zone := fifelse(
    latitude >= 15, "north",
    fifelse(latitude <= -15, "south", "equatorial")
  )]
  
  # Reassign fine zone (big_zone) from Latitude (on *peaks*)
  peaks[, big_zone := fcase(
    abs(latitude) < 15,                        "tropical",
    
    latitude >= 15 & latitude < 45,            "north_subtropical",
    latitude <= -15 & latitude > -45,          "south_subtropical",
    
    latitude >= 45 & latitude < 60,            "north_subpolar",
    latitude <= -45 & latitude > -60,          "south_subpolar",
    
    latitude >= 60,                            "north_polar",
    latitude <= -60,                           "south_polar",
    
    default = NA_character_
  )]
  
  classify_one <- function(zones) {
    if (length(zones) > 0 && all(zones == "tropical")) return("tropical")
    if (all(c("north_polar", "south_polar") %in% zones)) return("bipolar")
    if (all(c("north_subpolar", "south_subpolar") %in% zones)) return("bi_subpolar")
    if (all(c("north_subtropical", "south_subtropical") %in% zones)) return("bi_subtropical")
    north_high <- any(zones %in% c("north_polar", "north_subpolar"))
    south_high <- any(zones %in% c("south_polar", "south_subpolar"))
    if (north_high && south_high) return("bi_high_latitudes")
    return("unclassified")
  }
  
  out <- peaks[, .(
    category = classify_one(big_zone),
    top_peaks_summary = paste0(
      sprintf("lat=%.1f° (fit=%.3g)", latitude, peak_fit),
      collapse = "; "
    )
  ), by = taxon_id]
  
  return(out[])
}


## --------------------------------------------------------------------------------- ##
##  tag_exclusivity_votus_multi()                                                    ##
##  • A taxon receives one exclusivity / preferential tag based on                   ##
##    observations_input.                                                            ##
##  • Windows can overlap: polar and subpolar are both subsets of high-latitude.     ##
##  • Returns:                                                                       ##
##        $classified   → original table + ‘exclusivity_tags’ list-col               ##
##        $tag_counts   → data-table with counts per tag                             ##
## --------------------------------------------------------------------------------- ##
tag_exclusivity_multi <- function(classified, obs_input) {
  library(data.table)
  
  ## ---- 0. prep ---------------------------------------------------------- ##
  classified <- as.data.table(classified)
  obs <- as.data.table(observations_input)[
    norm_coverage > 0 & is.finite(latitude),
    .(taxon_id, lat = latitude)
  ]
  
  ## ---- 1. define tag windows & eligible categories --------------------- ##
  tag_spec <- list(
    polar = list(
      lower = 60, upper = Inf,
      cats  = c("bipolar"),
      label = "polar"
    ),
    
    subpolar = list(
      lower = 45, upper = 60,
      cats  = c("bi_subpolar"),
      label = "subpolar"
    ),
    
    subtropical = list(
      lower = 15, upper = 45,
      cats  = c("bi_subtropical"),
      label = "subtropical"
    ),
    
    high_latitude = list(
      lower = 45, upper = Inf,
      cats  = c("bipolar", "bi_subpolar", "bi_high_latitudes"),
      label = "high_latitude"
    ),
    
    tropical = list(            
      lower = 0,  upper = 15,
      cats  = c("tropical"),
      label = "tropical"
    )
  )
  
  ## ---- 2. build a long table of tags ----------------------------------- ##
  tag_rows <- rbindlist(lapply(names(tag_spec), function(tag) {
    win  <- tag_spec[[tag]]
    ids  <- classified[category %chin% win$cats, taxon_id]
    if (!length(ids)) return(NULL)
    
    tmp  <- obs[taxon_id %in% ids][ , .(
      exclusive = all(abs(lat) >= win$lower & abs(lat) < win$upper)
    ), by = taxon_id]
    
    tmp[ , tag := ifelse(exclusive,
                         paste0("exclusive_",    tag),
                         paste0("preferential_", tag))]
    tmp[ , .(taxon_id, tag)]
  }))
  
  ## ---- 3. attach tags back (list-col) ---------------------------------- ##
  tag_list <- tag_rows[ , .(exclusivity_tags = list(tag)), by = taxon_id]
  out      <- merge(classified, tag_list,
                    by = "taxon_id", all.x = TRUE, sort = FALSE)
  
  ## ---- 4. handy count table -------------------------------------------- ##
  counts <- tag_rows[ , .N, by = tag][order(tag)]
  
  return(list(classified = out, tag_counts = counts))
}


## ============================================================
## 5. Find tallest peak in each hemisphere
## ============================================================

top2_hemi <- peaks_merged_global[ 
  , .SD[which.max(peak_fit)],
  by = .(taxon_id, hemi = latitude > 0),
  .SDcols = c("latitude",          
              "zone",
              "peak_fit",
              "n_pos",
              "fit_ratio")
]#[!is.na(hemi)]                           # drop genomes missing N or S

# quick check
stopifnot("latitude" %chin% names(top2_hemi))

## ============================================================
## 6. Assign zones
## ============================================================

top2_hemi[, big_zone := fcase(
  abs(latitude) < tropical_cutoff,
  "tropical",
  
  latitude >= tropical_cutoff & latitude < subpolar_cutoff,
  "north_subtropical",
  
  latitude <= -tropical_cutoff & latitude > -subpolar_cutoff,
  "south_subtropical",
  
  latitude >= subpolar_cutoff & latitude < polar_cutoff,
  "north_subpolar",
  
  latitude <= -subpolar_cutoff & latitude > -polar_cutoff,
  "south_subpolar",
  
  latitude >= polar_cutoff,
  "north_polar",
  
  latitude <= -polar_cutoff,
  "south_polar",
  
  default = NA_character_
)]

## ============================================================
## 7. Execute functions
## ============================================================

tax_classification <- classify_latitude(top2_hemi)
res             <- tag_exclusivity_multi(tax_classification, obs_input)
classification <- res$classified
tag_counts      <- res$tag_counts
#tag_counts

## ============================================================
## 8. Prepare output
## ============================================================

# Pull all unclassified genomes -------------------------------------------
unclassified_ids <- classification[category == "unclassified",
                                    taxon_id]

#length(unclassified_ids)   # how many unclassified?
#head(unclassified_ids)     # preview a few


# Explode only the rows whose tag list length > 0
tag_map <- classification[
  lengths(exclusivity_tags) > 0,          # <- substitute for is.na()
  .(tag = unlist(exclusivity_tags)),      # unnest
  by = taxon_id
]

#head(tag_map)
##          taxon_id                  tag
## 1: ...   exclusive_polar
## 2: ...   exclusive_high_latitude
## ...

exclusive_bipolar_ids <- tag_map[tag == "exclusive_polar", taxon_id]
exclusive_high_lat_ids <- tag_map[tag == "exclusive_high_latitude", taxon_id]
exclusive_subpolar_ids <- tag_map[tag == "exclusive_subpolar", taxon_id]
exclusive_subtropical_ids <- tag_map[tag == "exclusive_subtropical", taxon_id]

preferential_bipolar_ids <- tag_map[tag == "preferential_polar", taxon_id]
preferential_high_lat_ids <- tag_map[tag == "preferential_high_latitude", taxon_id]
preferential_subpolar_ids <- tag_map[tag == "preferential_subpolar", taxon_id]
preferential_subtropical_ids <- tag_map[tag == "preferential_subtropical", taxon_id]

exclusive_tropical_ids <- tag_map[tag == "exclusive_tropical", taxon_id]
preferential_tropical_ids <- tag_map[tag == "preferential_tropical", taxon_id]

# Put every vector in a named list
tag_sets <- list(
  exclusive_bipolar        = exclusive_bipolar_ids,
  exclusive_high_latitude  = exclusive_high_lat_ids,
  exclusive_subpolar       = exclusive_subpolar_ids,
  exclusive_subtropical    = exclusive_subtropical_ids,
  
  preferential_bipolar     = preferential_bipolar_ids,
  preferential_high_latitude = preferential_high_lat_ids,
  preferential_subpolar    = preferential_subpolar_ids,
  preferential_subtropical = preferential_subtropical_ids,
  
  exclusive_tropical       = exclusive_tropical_ids,
  preferential_tropical    = preferential_tropical_ids
)

# Deduplicate dual classifications, i.e. applicable only to polar + high latitude and subpolar + high latitude
bihem_dt <- tag_map[
  ,
  .(
    final_tag = paste(sort(unique(tag)), collapse = "|")
  ),
  by = taxon_id
]

bihem_dt[, tag := final_tag]
bihem_dt[, source := "bihemispherical"]

## ============================================================
## 9. Save temp files
## ============================================================

if (!dir.exists(temp_dir)) {
  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
}

tag_map_file <- file.path(temp_dir, "tag_map.fst")
tag_sets_file <- file.path(temp_dir, "tag_sets.rds")
unclassified_ids_file <- file.path(temp_dir, "unclassified_ids.rds")
bihem_dt_file <- file.path(temp_dir, "bihemispherical_vOTUs.rds")

write_fst(tag_map, tag_map_file)
saveRDS(tag_sets, tag_sets_file)
saveRDS(unclassified_ids, unclassified_ids_file)
saveRDS(bihem_dt, bihem_dt_file)

message("Saved bimodal tag map to: ", tag_map_file)
message("Saved tag sets to: ", tag_sets_file)
message("Saved unclassified IDs to: ", unclassified_ids_file)
message("Saved bihemispherical vOTUs to: ", bihem_dt_file)

