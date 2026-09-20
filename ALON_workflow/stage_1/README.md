Stage 1 uses generalized additive models (GAMs) to estimate the latitudinal abundance distribution of individual taxa.

Note: Stage 1 can be skipped if users prefer to model per-latitude abundance independently. To run Stage 2, provide a compatible prediction table. Stage 2 was developed and tested using GAM-derived predictions.

For each `taxon_id`, normalized abundance (`norm_coverage`) across samples is modeled as a function of latitude. Predictions are generated from -90° to 90° latitude at 1° intervals.

Recommended pre-processing: Taxa should be filtered to the desired minimum number of detections before running this script. In our analysis, taxa were required to have at least five non-zero detections before GAM fitting.

The script can run either:
- on a SLURM high-performance computing cluster using an array job, or
- locally without SLURM, in which case all chunks are processed sequentially.

Required R packages:
- data.table
- mgcv
- fst


Optional:
- RhpcBLASctl
`RhpcBLASctl` is used for thread control when available but is not required.

Required input and columns (TSV file)
1. observations_input
	- `sample_name` - unique sample identifier
	- `taxon_id` - taxon identifier
	- `latitude` - sample latitude in decimal degrees
	- `norm_coverage` - normalized abundance/coverage

Additional columns are allowed and will be ignored.

For each taxon:

1. Sample-level normalized coverage is used without latitude binning.
2. Positive coverage values are used to calculate an abundance cap.
3. Coverage is capped at the specified quantile (`cap_quant`; default 0.995).
4. For `gaussian_log1p`, capped coverage is transformed using `log1p`.
5. A single smooth GAM is fitted across latitude.
6. Predictions are generated from -90° to 90° at 1° intervals.
7. Approximate 95% pointwise confidence intervals are calculated for the
   fitted curve and back-transformed to the original abundance scale.

Parameters:
Parameter 		Default 			Description
`chunk_size`	2000				Number of unique `taxon_id`s processed per chunk
`engine`		`gam`				GAM fitting engine (`gam` or `bam` for Gaussian models)
`k_spline` 		12 					Maximum basis dimension for the latitude smooth
`lat_min_train` -90					Minimum latitude included when fitting the model
`gamma` 		1					GAM smoothing penalty multiplier
`cap_quant` 	0.995				Quantile used to cap positive coverage values
`family`  		`gaussian_log1p` 	Model family: `gaussian_log1p` or `tweedie`

Default configuration:
CHUNK=2000
ENGINE=gam
K=12
LATTRIM=-90
GAMMA=1
CAPQ=0.995
FAMILY=gaussian_log1p

Output files:
1. Intermediate prediction chunks: `pred_chunk_0001.fst`, `pred_chunk_0002.fst`
2. Final merged predictions, `predictions_all.tsv` containing the following columns:
	- taxon_id
	- latitude
	- fit - predicted abundance
	- `lo` - lower approximate 95% confidence limit
	- `high` - upper approximate 95% confidence limit
	- `status` - model fitting status
3. `failed_taxa.tsv` containing a list of `taxa_id` for which a succesful GAM prediction was not obtained.

Status values:
	- `ok` - model fitted succesfully
	- `insufficient_data` - too few observations available
	- `insufficient latitudes` - fewer than three unique training latitudes
	- `fit_failed` - GAM fitting failed
	- `error` - another unexpected error occurred during processing

Usage:

IN=path/to/input.tsv
OUT=path/to/output/directory
CHUNK=2000 
ENGINE=gam
K=12
LATTRIM=-90
GAMMA=1
CAPQ=0.995
FAMILY=gaussian_log1p

Rscript predict_gam.R \
  "path/to/input.tsv" \
  "path/to/output/directory" \
  "2000" \
  "gam" \
  "12" \
  "-90" \
  "1" \
  "0.995" \
  "gaussian_log1p"

Local execution
If `SLURM_ARRAY_TAX_ID` is not detected, the script automatically runs in local mode and processes all chunks sequentially. Completed `.fst` are chunks are skipped, allowing an interrupted run to be resumed by running the same command again.

SLURM execution
When `SLURM_ARRAY_TAX_ID` is present, each SLURM array task processes one chunk. The script also reads `SLURM_CPUS_PER_TASK` to determine the number of available CPUs.

After each array task finishes its chunk, it checks whether all expected chunks have been generated. The final task to detect all completed chunks merges them into `predictions_all.tsv`.

A merge lock prevents multiple array tasks from performing the final merge simultaneously.

**Test data**
The eukaryotic all-fraction dataset `stage_1/test_input/input_euks.tsv` is used to test the full Stage 1 to Stage 2 workflow on a smaller input dataset. 
