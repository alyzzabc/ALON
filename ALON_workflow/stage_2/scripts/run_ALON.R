#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(data.table)
})

option_list <- list(
  make_option(
    "--pred_input",
    type = "character",
    help = paste(
      "Prediction input TSV(s) with columns: taxon_id, latitude, fit.",
      "For multiple size fractions, provide comma-separated paths."
    )
  ),
  make_option(
    "--obs_input",
    type = "character",
    help = paste(
      "Observation input TSV(s) with columns: taxon_id, latitude, norm_coverage.",
      "For multiple size fractions, provide comma-separated paths."
    )
  ),
  make_option(
    "--size_fraction",
    type = "character",
    default = NULL,
    help = paste(
      "Optional comma-separated size-fraction names.",
      "Example: all,cellular,virome.",
      "If omitted, names are inferred as input_1, input_2, ..."
    )
  ),
  make_option(
    "--reference_obs_input",
    type = "character",
    default = NULL,
    help = "Observation file used to recalculate final exclusive/preferential labels during size-fraction reconciliation. Required when multiple size fractions are provided."
  ),
  make_option(
    "--output_dir",
    type = "character",
    help = "Output directory"
  ),
  make_option(
    "--config_file",
    type = "character",
    default = NULL,
    help = "Optional config TSV with columns: parameter, value"
  ),
  make_option(
    "--script_dir",
    type = "character",
    default = ".",
    help = paste(
      "Directory containing predict_peaks.R, classify_bihemispherical.R,",
      "classify_monohemispherical.R, and reconcile_size_fractions.R"
    )
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

## required arguments
if (is.null(opt$pred_input)) {
  stop("--pred_input is required", call. = FALSE)
}

if (is.null(opt$obs_input)) {
  stop("--obs_input is required", call. = FALSE)
}

if (is.null(opt$output_dir)) {
  stop("--output_dir is required", call. = FALSE)
}

## parse possibly multiple inputs
prediction_files <- trimws(strsplit(opt$pred_input, ",")[[1]])
observation_files <- trimws(strsplit(opt$obs_input, ",")[[1]])

if (length(prediction_files) != length(observation_files)) {
  stop(
    "--pred_input and --obs_input must contain the same number of files.",
    "\nNumber of prediction files: ", length(prediction_files),
    "\nNumber of observation files: ", length(observation_files),
    call. = FALSE
  )
}

if (is.null(opt$size_fraction)) {
  size_fraction_names <- paste0("input_", seq_along(prediction_files))
} else {
  size_fraction_names <- trimws(strsplit(opt$size_fraction, ",")[[1]])
  
  if (length(size_fraction_names) != length(prediction_files)) {
    stop(
      "--size_fraction must contain the same number of names as input files.",
      "\nNumber of size-fraction names: ", length(size_fraction_names),
      "\nNumber of input files: ", length(prediction_files),
      call. = FALSE
    )
  }
}

if (length(prediction_files) > 1L) {
  if (is.null(opt$reference_obs_input)) {
    stop(
      "When multiple size fractions are provided, --reference_obs_input is required. ",
      "This file is used to recalculate final exclusive/preferential labels after reconciliation."
    )
  }
  
  if (!file.exists(opt$reference_obs_input)) {
    stop("--reference_obs_input file not found: ", opt$reference_obs_input)
  }
  
  reconcile_observations_input <- readr::read_tsv(
    opt$reference_obs_input,
    show_col_types = FALSE
  )
}

## base paths
base_output_dir <- opt$output_dir
config_file <- opt$config_file
script_dir <- opt$script_dir

## default latitudinal cutoffs
latzone_cutoffs <- list(
  tropical = 15,
  subpolar = 45,
  polar = 60
)

## allow config file to override defaults
if (!is.null(config_file)) {
  config_dt <- data.table::fread(config_file)
  
  if (!all(c("parameter", "value") %in% names(config_dt))) {
    stop(
      "Config file must contain columns: parameter, value",
      call. = FALSE
    )
  }
  
  get_config_value <- function(parameter_name, default_value) {
    value <- config_dt[parameter == parameter_name, value]
    
    if (length(value) == 0L) {
      return(default_value)
    }
    
    if (length(value) > 1L) {
      stop(
        "Config parameter appears more than once: ",
        parameter_name,
        call. = FALSE
      )
    }
    
    as.numeric(value)
  }
  
  latzone_cutoffs$tropical <- get_config_value(
    "tropical_cutoff",
    latzone_cutoffs$tropical
  )
  
  latzone_cutoffs$subpolar <- get_config_value(
    "subpolar_cutoff",
    latzone_cutoffs$subpolar
  )
  
  latzone_cutoffs$polar <- get_config_value(
    "polar_cutoff",
    latzone_cutoffs$polar
  )
}

## validate cutoffs
if (any(is.na(unlist(latzone_cutoffs)))) {
  stop(
    "Latitudinal cutoffs must be numeric.",
    call. = FALSE
  )
}

if (!(latzone_cutoffs$tropical < latzone_cutoffs$subpolar &&
      latzone_cutoffs$subpolar < latzone_cutoffs$polar)) {
  stop(
    "Latitudinal cutoffs must satisfy: ",
    "tropical_cutoff < subpolar_cutoff < polar_cutoff.",
    call. = FALSE
  )
}

message(
  "Latitudinal cutoffs: tropical/equatorial < ",
  latzone_cutoffs$tropical,
  ", subpolar/high-latitude >= ",
  latzone_cutoffs$subpolar,
  ", polar >= ",
  latzone_cutoffs$polar
)

if (!is.null(config_file) && !file.exists(config_file)) {
  stop("Config file does not exist: ", config_file, call. = FALSE)
}

dir.create(base_output_dir, recursive = TRUE, showWarnings = FALSE)

## source worker scripts
predict_script <- file.path(script_dir, "predict_peaks.R")
bihemispherical_script <- file.path(script_dir, "classify_bihemispherical.R")
monohemispherical_script <- file.path(script_dir, "classify_monohemispherical.R")
reconcile_script <- file.path(script_dir, "reconcile_size_fractions.R")

required_scripts <- c(
  predict_script,
  bihemispherical_script,
  monohemispherical_script
)

for (script in required_scripts) {
  if (!file.exists(script)) {
    stop("Required script does not exist: ", script, call. = FALSE)
  }
}

if (length(prediction_files) > 1L && !file.exists(reconcile_script)) {
  stop(
    "Multiple size fractions were provided, but reconcile_size_fractions.R does not exist: ",
    reconcile_script,
    call. = FALSE
  )
}

classification_dirs <- character()

for (i in seq_along(prediction_files)) {
  predictions_file <- prediction_files[i]
  observations_file <- observation_files[i]
  size_fraction <- size_fraction_names[i]
  
  if (!file.exists(predictions_file)) {
    stop("Prediction input file does not exist: ", predictions_file, call. = FALSE)
  }
  
  if (!file.exists(observations_file)) {
    stop("Observation input file does not exist: ", observations_file, call. = FALSE)
  }
  
  message("============================================================")
  message("Running size fraction: ", size_fraction)
  message("Prediction input: ", predictions_file)
  message("Observation input: ", observations_file)
  
  ## Preserve old behavior for single-input runs.
  ## For multiple inputs, write each fraction to its own subfolder.
  if (length(prediction_files) == 1L) {
    output_dir <- base_output_dir
  } else {
    output_dir <- file.path(base_output_dir, size_fraction)
  }
  
  temp_dir <- file.path(output_dir, "temp")
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
  
  message("Output directory: ", output_dir)
  message("Temporary directory: ", temp_dir)
  
  if (!is.null(config_file)) {
    message("Config file: ", config_file)
  } else {
    message("No config file provided. Using defaults.")
  }
  
  message("Reading prediction input: ", predictions_file)
  predictions_input <- read_tsv(
    predictions_file,
    show_col_types = FALSE
  )
  
  message("Reading observation input: ", observations_file)
  observations_input <- read_tsv(
    observations_file,
    show_col_types = FALSE
  )
  
  message("Running predict_peaks.R")
  source(predict_script)
  
  message("Running classify_bihemispherical.R")
  source(bihemispherical_script)
  
  message("Running classify_monohemispherical.R")
  source(monohemispherical_script)
  
  classification_dirs <- c(classification_dirs, output_dir)
  
  message("Finished size fraction: ", size_fraction)
}

if (length(classification_dirs) > 1L && !exists("reconcile_observations_input")) {
  stop("Reference observations were not loaded. Check --reference_obs_input")
}

if (length(classification_dirs) > 1L) {
  message("============================================================")
  message("Running reconciliation across size fractions")
  
  output_dir <- base_output_dir
  
  map_bihem_polar_to_bipolar <- TRUE
  
  source(reconcile_script)
  
  message("Reconciliation complete.")
} else {
  message("Only one size fraction provided. Skipping reconciliation.")
}

message("Pipeline complete.")
message("Base output directory: ", base_output_dir)

if (length(classification_dirs) == 1L) {
  message("Final output: ", file.path(base_output_dir, "biogeog_out.tsv"))
  message("Unclassified IDs: ", file.path(base_output_dir, "unclassified_taxid.txt"))
} else {
  message("Reconciled output: ", file.path(base_output_dir, "reconciled"))
}