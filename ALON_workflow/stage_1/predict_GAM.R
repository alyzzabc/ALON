#!/usr/bin/env Rscript
options(warn = 1)

suppressPackageStartupMessages({
  library(data.table)
  library(mgcv)
  library(fst)
})

# ---------- CLI ----------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop(paste(
    "Usage:",
    "general_final_predict_gam.R <input.tsv> <outdir>",
    "[chunk_size] [engine] [k_spline] [lat_min_train]",
    "[lat_min_pred] [gamma] [cap_quant] [family_mode]",
    "[min_detections]",
    sep = "\n  "
  ))
}

infile         <- args[[1]]
outdir         <- args[[2]]
chunk_size     <- if (length(args) >= 3  && nzchar(args[[3 ]])) as.integer(args[[3 ]]) else 2000L
engine         <- if (length(args) >= 4  && nzchar(args[[4 ]])) args[[4]] else "bam"       # "bam" or "gam"
k_spline_in    <- if (length(args) >= 5  && nzchar(args[[5 ]])) as.integer(args[[5 ]]) else 12

# refit-tuning knobs
lat_min_train    <- if (length(args) >= 6 && nzchar(args[[6]])) as.numeric(args[[6]]) else -90
lat_min_pred <- if (length(args) >= 7 && nzchar(args[[7]])) {
  as.numeric(args[[7]])
} else {
  -77
}

gamma_in         <- if (length(args) >= 8 && nzchar(args[[8]])) as.numeric(args[[8]]) else 1
cap_quant_in     <- if (length(args) >= 9 && nzchar(args[[9]])) as.numeric(args[[9]]) else 0.995
family_mode <- if (length(args) >= 10 && nzchar(args[[10]])) args[[10]] else "gaussian_log1p"
min_detections <- if (length(args) >= 11 && nzchar(args[[11]])) {
  as.integer(args[[11]])
} else {
  1L
}

family_mode <- tolower(family_mode)
if (!family_mode %in% c("gaussian_log1p", "tweedie")) {
  stop("family_mode must be 'gaussian_log1p' or 'tweedie'")
}

k_spline      <- if (is.na(k_spline_in) || k_spline_in < 3) 12L else as.integer(k_spline_in)
gamma_val     <- if (is.na(gamma_in)) 1.0 else gamma_in
cap_quant     <- if (is.na(cap_quant_in)) 0.995 else cap_quant_in

# ---------- execution mode ----------

slurm_mode <- nzchar(Sys.getenv("SLURM_ARRAY_TASK_ID"))

if (slurm_mode) {

  array_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  cpus     <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))

} else {

  array_id <- NA_integer_
  cpus     <- 1L

}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message(sprintf(
  "Minimum detections required: %d",
  min_detections
))

message(sprintf(
  "Settings: family_mode=%s • engine=%s • k_spline=%d • chunk_size=%d",
  family_mode, engine, k_spline, chunk_size
))

message(sprintf(
  "Refit tweaks: lat_min_train=%.1f • gamma=%.2f • cap_quant=%.3f",
  lat_min_train, gamma_val, cap_quant
))

# ---------- READ tsv file ----------

message(sprintf("Reading %s ...", infile))

epi <- fread(infile, sep = "\t")

required_cols <- c("sample_name", "latitude", "taxon_id", "norm_coverage")

if (!all(required_cols %in% names(epi))) {
  missing <- setdiff(required_cols, names(epi))
  stop(
    "Missing required input columns: ",
    paste(missing, collapse = ", ")
  )
}

epi <- epi[, ..required_cols]

## 
samples_lat <- unique(epi[, .(sample_name, latitude)])
setkey(samples_lat, sample_name)

# per (sample_name,taxon_id) norm_coverage — keep max
wrk_tbl <- epi[, .(norm_coverage = max(norm_coverage, na.rm = TRUE)), by = .(sample_name, taxon_id)]
setkey(wrk_tbl, sample_name, taxon_id)

# discard taxa with less than 5 detections
det_counts <- wrk_tbl[
  norm_coverage > 0,
  .(n_detections = uniqueN(sample_name)),
  by = taxon_id
]

keep_taxa <- det_counts[
  n_detections >= min_detections,
  taxon_id
]

wrk_tbl <- wrk_tbl[taxon_id %in% keep_taxa]

uniq_taxa <- sort(unique(wrk_tbl$taxon_id))

uniq_taxa <- sort(unique(wrk_tbl$taxon_id))
n_total   <- length(uniq_taxa)
n_chunks  <- ceiling(n_total / chunk_size)

# Decide which chunks this R process should handle

if (slurm_mode) {

  if (array_id > n_chunks) {
    message("Array index exceeds number of chunks. Exiting cleanly.")
    quit(save = "no", status = 0)
  }

  chunk_ids <- array_id

  message(sprintf(
    "Running in SLURM mode: array task %d of %d",
    array_id,
    n_chunks
  ))

} else {

  chunk_ids <- seq_len(n_chunks)

  message(sprintf(
    "Running in local mode: processing all %d chunks sequentially",
    n_chunks
  ))

}

message(sprintf("Total taxa: %d | chunk_size: %d | n_chunks: %d ",
                n_total, chunk_size, n_chunks, array_id))

# ---------- MODEL ----------
predict_one <- function(taxon,
                        samples_lat, wrk_tbl,
                        cap_quant   = 0.995,
                        k_spline    = 12L,
                        grid_by     = 1,
                        engine      = "gam",
                        discrete    = TRUE,
                        nthreads    = 1,
                        lat_min_train = -90,
                        lat_min_pred = -90,
                        gamma_val   = 1,
                        family_mode = "gaussian_log1p") {

  # Number of real positive detections for this taxon
  n_detect <- wrk_tbl[
    taxon_id == taxon & norm_coverage > 0,
    uniqueN(sample_name)
  ]
  
  if (n_detect < min_detections) {
    return(data.table(
      taxon_id = taxon,
      latitude = seq(lat_min_pred, 90, by = grid_by),
      fit = 0, lo = 0, hi = 0,
      status = "insufficient_detections",
      engine = engine,
      family_mode = family_mode,
      k_spline = k_spline,
      gamma = gamma_val,
      lat_min_train = lat_min_train
    ))
  }
  
  # Join norm_coverage for this taxon_id to every sample_name/latitude (fills zeros for non-detections)
  xdt <- merge(samples_lat, wrk_tbl[.(unique(samples_lat$sample_name), taxon)],
               by = "sample_name", all.x = TRUE, allow.cartesian = TRUE)
  xdt[is.na(norm_coverage), norm_coverage := 0]
  xdt[, lat := pmin(pmax(latitude, -90), 90)]
  xdt[, x := lat]

  # Cap value from positives only
  cap_val <- 0
  if (sum(xdt$norm_coverage > 0) > 0) {
    cap_val <- as.numeric(quantile(xdt[norm_coverage > 0, norm_coverage],
                                   probs = cap_quant, na.rm = TRUE))
    if (!is.finite(cap_val) || cap_val < 0) cap_val <- 0
  }

  # Always compute capped response
  xdt[, y_cap := pmin(pmax(norm_coverage, 0), cap_val)]

  # Gaussian uses log1p(y_cap); Tweedie uses y_cap directly
  xdt[, log_cov := log1p(y_cap)]

  if (family_mode == "tweedie") {
    xdt <- xdt[is.finite(x) & is.finite(y_cap)]
  } else {
    xdt <- xdt[is.finite(x) & is.finite(log_cov)]
  }

  if (nrow(xdt) < 5L) {
    return(data.table(
      taxon_id = taxon,
      latitude = seq(lat_min_pred, 90, by = grid_by),
      fit = 0, lo = 0, hi = 0,
      status = "insufficient_data",
      engine = engine,
      family_mode = family_mode,
      k_spline = k_spline,
      gamma = gamma_val, lat_min_train = lat_min_train
    ))
  }

  # ---- Trim southern tail in TRAINING ONLY (still predict full grid) ----
  train <- xdt[lat >= lat_min_train]

  # Choose response column for modeling, then store as 'y'
  ycol <- if (family_mode == "tweedie") "y_cap" else "log_cov"

  train <- train[, .(
    x,
    y = get(ycol),
    w = 1
    )]

  n_ux <- uniqueN(train$x)
  if (n_ux < 3L) {
    return(data.table(
      taxon_id = taxon,
      latitude = seq(lat_min_pred, 90, by = grid_by),
      fit = 0, lo = 0, hi = 0,
      status = "insufficient_latitudes",
      engine = engine,
      family_mode = family_mode,
      k_spline = k_spline,
      gamma = gamma_val,
      lat_min_train = lat_min_train
    ))
  }

  # Prediction grid
  newx <- seq(lat_min_pred, 90, by = grid_by)
  fit_ok <- TRUE

  # Helper: fit model with chosen family; keep your existing bam/gam option for gaussian,
  # but for tweedie we use gam (simpler / safer).
  fit_model <- function(form, train) {
    if (family_mode == "tweedie") {
      mgcv::gam(form, data = train, weights = train$w,
                method = "REML", gamma = gamma_val,
                family = mgcv::tw(link = "log"))
    } else {
      if (engine == "bam") {
        mgcv::bam(form, data = train, weights = train$w, method = "fREML",
                  discrete = TRUE, nthreads = nthreads, gamma = gamma_val)
      } else {
        mgcv::gam(form, data = train, weights = train$w, method = "REML",
                  gamma = gamma_val)
      }
    }
  }

  # ---- Fit ----
  k_use <- max(3L, min(k_spline, n_ux - 1L))
  form  <- as.formula(sprintf("y ~ s(x, k = %d)", k_use))
  m <- try(fit_model(form, train), silent = TRUE)
  
  newd <- data.table(latitude = newx, x = newx)
  
  if (inherits(m, "try-error")) {
    fit_msg <- as.character(m)
    
    return(data.table(
      taxon_id = taxon,
      latitude = seq(lat_min_pred, 90, by = grid_by),
      fit = 0,
      lo = 0,
      hi = 0,
      status = paste0("fit_failed: ", fit_msg),
      engine = engine,
      family_mode = family_mode,
      k_spline = k_spline,
      gamma = gamma_val,
      lat_min_train = lat_min_train
    ))
  }

  if (!fit_ok) {
    return(data.table(
      taxon_id = taxon,
      latitude = newx,
      fit = 0, lo = 0, hi = 0,
      status = "fit_failed",
      engine = engine,
      family_mode = family_mode,
      k_spline = k_spline,
      gamma = gamma_val,
      lat_min_train = lat_min_train
    ))
  }

  # ---- Predict + back-transform to "norm_coverage-like" scale ----
  if (family_mode == "tweedie") {
    pp <- predict(m, newdata = newd, se.fit = TRUE, type = "link")
    fit <- pmax(exp(pp$fit), 0)
    lo  <- pmax(exp(pp$fit - 1.96 * pp$se.fit), 0)
    hi  <- pmax(exp(pp$fit + 1.96 * pp$se.fit), 0)
  } else {
    pp <- predict(m, newdata = newd, se.fit = TRUE)
    fit <- pmax(expm1(pp$fit), 0)
    lo  <- pmax(expm1(pp$fit - 1.96 * pp$se.fit), 0)
    hi  <- pmax(expm1(pp$fit + 1.96 * pp$se.fit), 0)
  }

  data.table(
    taxon_id    = taxon,
    latitude = newx,
    fit = fit, lo = lo, hi = hi,
    status = "ok",
    engine = engine,
    family_mode = family_mode,
    k_spline = k_spline,
    gamma = gamma_val, lat_min_train = lat_min_train
  )
}

# ---------- threading ----------
Sys.setenv(OMP_NUM_THREADS = cpus, OPENBLAS_NUM_THREADS = 1, MKL_NUM_THREADS = 1)
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(cpus)
}

# ---------- run this chunk ----------
for (chunk_id in chunk_ids) {

  i_start <- (chunk_id - 1L) * chunk_size + 1L
  i_end   <- min(chunk_id * chunk_size, n_total)

  taxa_chunk <- uniq_taxa[i_start:i_end]

  message(sprintf(
    "Processing chunk %d/%d: taxa %d..%d (n=%d)",
    chunk_id,
    n_chunks,
    i_start,
    i_end,
    length(taxa_chunk)
  ))

  outfile <- file.path(
    outdir,
    sprintf("pred_chunk_%04d.fst", chunk_id)
  )

  if (file.exists(outfile)) {
    message(sprintf(
      "Chunk %d already exists; skipping.",
      chunk_id
    ))
    next
  }

    res_list <- vector("list", length(taxa_chunk))

    for (i in seq_along(taxa_chunk)) {
    taxon <- taxa_chunk[[i]]
    if (i %% 50 == 0) message(sprintf(" ... %d/%d", i, length(taxa_chunk)))
    res_list[[i]] <- tryCatch(
        predict_one(taxon, samples_lat, wrk_tbl,
                    cap_quant     = cap_quant,
                    k_spline      = k_spline,
                    grid_by       = 1,
                    engine        = engine,
                    discrete      = TRUE,
                    nthreads      = cpus,
                    lat_min_train = lat_min_train,
                    lat_min_pred  = lat_min_pred,
                    gamma_val     = gamma_val,
                    family_mode   = family_mode),
        error = function(e) {
        message(sprintf("ERROR on %s: %s", taxon, conditionMessage(e)))
        data.table(
            taxon_id = taxon,
            latitude = seq(lat_min_pred, 90, 1),
            fit = 0, lo = 0, hi = 0,
            status = paste0("error: ", conditionMessage(e)),
            engine = engine,
            family_mode = family_mode,
            k_spline = k_spline,
            gamma = gamma_val, lat_min_train = lat_min_train
        )
        }
    )
}

pred_chunk <- rbindlist(
  res_list,
  use.names = TRUE,
  fill = TRUE
)

tmpfile <- file.path(
  outdir,
  sprintf("pred_chunk_%04d.fst.tmp", chunk_id)
)

outfile <- file.path(
  outdir,
  sprintf("pred_chunk_%04d.fst", chunk_id)
)

write_fst(
  as.data.frame(pred_chunk),
  tmpfile,
  compress = 50
)

if (!file.rename(tmpfile, outfile)) {
  stop(sprintf(
    "Failed to finalize chunk %d",
    chunk_id
  ))
}

message(sprintf(
  "Wrote chunk %d/%d: %d rows",
  chunk_id,
  n_chunks,
  nrow(pred_chunk)
))

}

# ---------- MERGE ALL CHUNKS WHEN COMPLETE ----------

chunk_files <- list.files(
  outdir,
  pattern = "^pred_chunk_[0-9]+\\.fst$",
  full.names = TRUE
)

message(sprintf(
  "Chunks complete: %d/%d",
  length(chunk_files),
  n_chunks
))

if (length(chunk_files) == n_chunks) {

  lockdir <- file.path(outdir, ".merge_lock")

  got_lock <- dir.create(
    lockdir,
    showWarnings = FALSE
  )

  if (got_lock) {

    tryCatch({

      message("All chunks detected. Starting final merge...")

      # Re-check after acquiring the lock
      chunk_files <- list.files(
        outdir,
        pattern = "^pred_chunk_[0-9]+\\.fst$",
        full.names = TRUE
      )

      if (length(chunk_files) != n_chunks) {
        stop(sprintf(
          "Chunk count changed after acquiring merge lock: found %d, expected %d",
          length(chunk_files),
          n_chunks
        ))
      }

      chunk_files <- sort(chunk_files)

      predictions <- rbindlist(
        lapply(chunk_files, read_fst),
        use.names = TRUE,
        fill = TRUE
      )

      # ---------- sanity checks ----------

      required_output_cols <- c(
        "taxon_id",
        "latitude",
        "fit",
        "lo",
        "hi",
        "status"
      )

      missing_output_cols <- setdiff(
        required_output_cols,
        names(predictions)
      )

      if (length(missing_output_cols) > 0) {
        stop(
          "Merged predictions are missing required columns: ",
          paste(missing_output_cols, collapse = ", ")
        )
      }

      n_tax <- uniqueN(predictions$taxon_id)
      n_lats <- uniqueN(predictions$latitude)

      message(sprintf(
        "Merged output: %d taxa | %d latitude grid points | %d rows",
        n_tax,
        n_lats,
        nrow(predictions)
      ))

      # ---------- status summary ----------

      status_summary <- unique(
        predictions[, .(
          taxon_id,
          status
        )]
      )[, .N, by = status][order(status)]

      message("Taxon fit status:")

      for (j in seq_len(nrow(status_summary))) {
        message(sprintf(
          "  %s: %d taxa",
          status_summary$status[j],
          status_summary$N[j]
        ))
      }

      # ---------- write failed taxa report ----------

      failed_taxa <- unique(
        predictions[
          status != "ok",
          .(
            taxon_id,
            status
          )
        ]
      )

      if (nrow(failed_taxa) > 0) {

        failed_tsv <- file.path(
          outdir,
          "failed_taxa.tsv"
        )

        fwrite(
          failed_taxa,
          failed_tsv,
          sep = "\t"
        )

        message(sprintf(
          "Wrote failed_taxa.tsv with %d taxa",
          nrow(failed_taxa)
        ))
      }

      # ---------- write final TSV ----------

      final_tsv <- file.path(
        outdir,
        "predictions_all.tsv"
      )

      tmp_tsv <- paste0(
        final_tsv,
        ".tmp"
      )

      fwrite(
        predictions,
        tmp_tsv,
        sep = "\t"
      )

      if (!file.rename(tmp_tsv, final_tsv)) {
        stop(sprintf(
          "Failed to rename temporary merged TSV: %s -> %s",
          tmp_tsv,
          final_tsv
        ))
      }

      message(sprintf(
        "Final merged TSV written: %s",
        final_tsv
      ))

    }, finally = {

      # Always release the lock, even if anything above fails
      if (dir.exists(lockdir)) {
        unlink(
          lockdir,
          recursive = TRUE,
          force = TRUE
        )
      }

      message("Merge lock released.")
    })
  }
}