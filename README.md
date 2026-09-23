# ALON (Abundance-based Latitudinal distributiON)

_Alon means "wave" in Tagalog, one of the over 170 languages spoken by Filipinos. Their ancestors were seafarers, who mastered navigating turbulent waters. But they did not tame the waves, for nature cannot be subdued. They learned the language of the waves, and the waves learned theirs._

#

ALON is an R-based framework designed to assign microbial biogeography using global scale -omics data. 
We applied this framework to microbial hosts and viruses using metagenomes from the ocean. 

The pipeline consists of two stages: Per-latitude abundance prediction using generalized additive models (GAMs) and classification of latitudinal distributions. The initial input requires a sample identifier, feature ID, latitude, and abundance. The final output is a tab-separated file with the feature ID, biogeographic classification, and the exclusivity status. The exclusivity status indicates whether the feature ID is found preferentially or exclusively in the assigned region. ALON can also take multiple size fractions in which case each fraction is first classified using its paired GAM prediction and raw observation files, and then integrated at the final reconciliation step.

#

## Stage 1
Stage 1 uses generalized additive models (GAMs) to estimate the latitudinal abundance distribution of individual taxa.

Note: Stage 1 can be skipped if users prefer to model per-latitude abundance independently. To run Stage 2, provide a compatible prediction table. Stage 2 was developed and tested using GAM-derived predictions.

For each `taxon_id`, normalized abundance (`norm_coverage`) across samples is modeled as a function of latitude. Predictions are generated from -90° to 90° latitude at 1° intervals.

The script can run either:
- on a SLURM high-performance computing cluster using an array job, or
- locally without SLURM, in which case all chunks are processed sequentially.

### Required R packages:
- data.table
- mgcv
- fst

Optional:
-  `RhpcBLASctl` - used for thread control when available but is not required

### Required input and columns (TSV file)
1. observations_input
	- `sample_name` - unique sample identifier
	- `taxon_id` - taxon identifier
	- `latitude` - sample latitude in decimal degrees
	- `norm_coverage` - normalized abundance/coverage

Additional columns are allowed and will be ignored.

For each taxon:

1. Sample-level normalized coverage is used without latitude binning.
2. Positive coverage values are used to calculate an abundance cap.
3. Coverage is capped at the specified quantile (`cap_quant`; default: 0.995).
4. For `gaussian_log1p`, capped coverage is transformed using `log1p`.
5. A single smooth GAM is fitted across latitude.
6. Predictions are generated from -90° to 90° at 1° intervals.
7. Approximate 95% pointwise confidence intervals are calculated for the
   fitted curve and back-transformed to the original abundance scale.

### Parameters:

| Parameter | Default | Description |
|---|---:|---|
| `chunk_size` | 2000 | Number of unique `taxon_id`s processed per chunk |
| `engine` | `gam` | GAM fitting engine (`gam` or `bam` for Gaussian models) |
| `k_spline` | 12 | Maximum basis dimension for the latitude smooth |
| `lat_min_train` | -90 | Southern latitude cutoff below which observed samples are excluded from fitting |
| `lat_min_pred` | -90 | Southern latitude cutoff below which model predictions are not generated |
| `gamma` | 1 | GAM smoothing penalty multiplier |
| `cap_quant` | 0.995 | Quantile used to cap positive coverage values |
| `family` | `gaussian_log1p` | Model family: `gaussian_log1p` or `tweedie` |
|`min_detections`| 1 | Minimum number of positive detections per taxon required for GAM fitting |

### Output files:
1. Intermediate prediction chunks: `pred_chunk_0001.fst`, `pred_chunk_0002.fst`
2. Final merged predictions, `predictions_all.tsv` containing the following columns:
	- `taxon_id`
	- `latitude`
	- `fit` - predicted abundance
	- `lo` - lower approximate 95% confidence limit
	- `high` - upper approximate 95% confidence limit
	- `status` - model fitting status
3. `failed_taxa.tsv` containing a list of `taxa_id` for which a succesful GAM prediction was not obtained.

Status values:
- `ok` - model fitted succesfully
- `insufficient_data` - too few observations available
- `insufficient_latitudes` - fewer than three unique training latitudes
- `fit_failed` - GAM fitting failed
- `error` - another unexpected error occurred during processing

### Usage:
```bash
Rscript predict_gam.R \
  "path/to/input.tsv" \
  "path/to/output/directory" \
  "2000" \
  "gam" \
  "12" \
  "-90" \
  "-90" \
  "1" \
  "0.995" \
  "gaussian_log1p" \
  "1"
```

#### Local execution
If `SLURM_ARRAY_TAX_ID` is not detected, the script automatically runs in local mode and processes all chunks sequentially. Completed `.fst` are chunks are skipped, allowing an interrupted run to be resumed by running the same command again.

#### SLURM execution
When `SLURM_ARRAY_TAX_ID` is present, each SLURM array task processes one chunk. The script also reads `SLURM_CPUS_PER_TASK` to determine the number of available CPUs.

After each array task finishes its chunk, it checks whether all expected chunks have been generated. The final task to detect all completed chunks merges them into `predictions_all.tsv`.

A merge lock prevents multiple array tasks from performing the final merge simultaneously.

### Test data:
The eukaryotic all-fraction dataset `stage_1/test_input/input_euks_small.tsv` is used to test the full Stage 1 to Stage 2 workflow on a smaller input dataset. 

#

## Stage 2
Stage 2 takes taxon_id, latitude, and norm_coverage (observed abundance) and fit (predicted abundance) at given latitude and classifies the latitudinal distribution of taxa. Latitudinal distribution is assigned based on the location of the tallest predicted peak/s. Classification could either be bihemispherical (polar, subpolar, subtropical, high latitude) or monohemispherical (tropical, polar N/S, subpolar N/S, subtropical N/S). Exclusivity is based on actual detection (norm_coverage > 0), e.g. a taxon classified as exclusive_polar peaks in abundance in the polar regions and is never detected outside the polar regions, and a taxon classified as preferential_polar peaks in abundance in the polar regions, but is detected in relatively lower abundance outside of the polar regions. ALON can also take multiple size fractions. See below for usage.

### Required R packages:
- oftparse
- fst
- readr
- data.table
- dplyr
- tidyr

### Required input and columns:
1. observations_input (TSV file)
	- `taxon_id` - taxon identifier
	- `latitude` - sample latitude in decimal degrees
	- `norm_coverage`
  
2. predictions_input (TSV file)
	- `taxon_id` - taxon identifier
	- `latitude` - sample latitude in decimal degrees
	- `fit` (predicted norm_coverage)

3. reference_observation_input (TSV file, required only when providing multiple size fractions)
   - `taxon_id` - taxon identifier
   - `latitude` - sample latitude in decimal degrees
   - `norm_coverage`

**Note in analyzing multiple size fractions:**
During the final reconciliation step, consensus geographic assignments are derived across fractions, and exclusivity status is recalculated from the `reference_obs_input`. This is typically the richest dataset, i.e. when analyses included both all-fraction combined and cellular-fraction inputs, exclusivity was recalculated using the all-fraction observations. Reconciliation was validated for up to three size fractions which used a minimum two-fraction support rule by default For analyses with more than three fractions, users may manually adjust the rule in `reconcile_size_fractions.R`.

***Optional:*** Config can be provided by the user. Default config has been validated for abundance prediction using GAM.
| Parameter | Default | Description |
|---|---|---|
|`min_support`|2|minimum number of actual detections (norm_coverage > 0) required to support a peak |
|`min_rel_height`|0.5|minimum relative height of a secondary peak compared with the tallest peak for the same taxon |
|`min_abs_fit`|0.005|minimum absolute fitted abundance value required for a peak to be considered real |
|`min_sep_deg`|5|minimum distance in latitudinal degrees used when deciding whether nearby peaks should be merged |
|`tropical_cutoff`|15|absolute latitude below which observations are considered tropical/equatorial |
|`subpolar_cutoff`|45|absolute latitude below at which subpolar/high latitude zones begin |
|`polar_cutoff`|60|absolute latitude at which polar zones begin  |

Latitudinal cutoffs must satisfy `tropical_cutoff` < `subpolar_cutoff` < `polar_cutoff`. By default, ALON uses 15°, 45°, and 60° absolute latitude.


Example of optional tab-separated file with columns `parameter` and `value`:

```text
parameter	value
min_support	2
min_rel_height	0.5
min_abs_fit	0.005
min_sep_deg	5
tropical_cutoff 15
subpolar_cutoff 45
polar_cutoff 60
```

### Output files
1. `biogeog_out.tsv` containing the following columns:
  - `taxon_id`
  - `final_geo` - assigned latitudinal category
  - `excl_status`- `exclusive` or `preferential` 
  - `final_tag` -  `excl_status` + `final_geo`

2. `unclassified_taxid.txt` containing a list of `taxa_id` for which a laittudinal distribution was not assigned.

3. `reconciled/reconciled_biogeography_combined.tsv` for multi-fraction analyses and contains the following columns:
  - `taxon_id`
  - `<fraction>_geo`- assigned latitudinal category for each size fraction
  - `<fraction>_tag` - `exclusive` or `preferential` for each size fraction, based on its paired observation input
  - `final_geo` - reconciled latitudinal category
  - `excl_status`- `exclusive` or `preferential`, recalculated based on `reference_obs_input`
  - `final_tag` -  `excl_status` + `final_geo`

### Usage:

```bash
Rscript run_ALON.R \
  --pred_input predictions_all.tsv \
  --obs_input input_euks_small.tsv \
  --output_dir output \
  --script_dir /path/to/scripts
```

Or provide own config file:

```bash
Rscript run_ALON.R \
  --pred_input predictions_all.tsv \
  --obs_input input_euks_small.tsv \
  --output_dir output \
  --config_file config_file.tsv \
  --script_dir /path/to/scripts
```

Or run multiple fractions:
```bash
Rscript run_ALON.R \
  --pred_input gam_predictions_ALL_first100.tsv,gam_predictions_CELLULAR_first100.tsv,gam_predictions_VIRUS_first100.tsv \
  --obs_input observations_ALL_first100.tsv,observations_CELLULAR_first100.tsv,observations_VIRUS_first100.tsv \
  --size_fraction ALL,CELLULAR,VIRUS \
  --reference_obs_input observations_ALL_first100.tsv \
  --output_dir output \
  --script_dir /path/to/scripts
```

### Test data:
The eukaryotic all-fraction inputs in `stage_2/test_input/all_fractions_combined/` are the same as the Stage 1 input and output. A separate multiple-fraction dataset, including all fractions combined, cellular fraction, and viral fraction, is provided in `stage_2/test_input/size_fraction_reconciliation/` to test size fraction reconciliation. 
