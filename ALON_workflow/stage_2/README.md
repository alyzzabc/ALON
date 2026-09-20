Stage 2 takes taxon_id, latitude, and norm_coverage (observed abundance) and fit (predicted abundance) at given latitude and classifies the latitudinal distribution of taxa. Latitudinal distribution is assigned based on the location of the tallest predicted peak/s. Classification could either be bihemispherical (polar, subpolar, subtropical, high latitude) or monohemispherical (tropical, polar N/S, subpolar N/S, subtropical N/S). Exclusivity is based on actual detection (norm_coverage > 0), e.g. a taxon classified as exclusive_polar peaks in abundance in the polar regions and is never detected outside the polar regions, and a taxon classified as preferential_polar peaks in abundance in the polar regions, but is detected in relatively lower abundance outside of the polar regions. ALON can also take multiple size fractions. See below for usage.

Required R packages:
oftparse
fst
readr
data.table
dplyr
tidyr

Required input and columns:
1. observations_input (TSV file)
	- `taxon_id` - taxon identifier
	- `latitude` - sample latitude in decimal degrees
	- norm_coverage
  
2. predictions_input (TSV file)
	- `taxon_id` - taxon identifier
	- `latitude` - sample latitude in decimal degrees
	- `fit` (predicted norm_coverage)

3. reference_observation_input (TSV file, required only when providing multiple size fractions)
  - `taxon_id` - taxon identifier
	- `latitude` - sample latitude in decimal degrees
	- norm_coverage

Note in analyzing multiple size fractions:
During the final reconciliation step, consensus geographic assignments are derived across fractions, and exclusivity status is recalculated from the `reference_obs_input`. This is typically the richest dataset, i.e. when analyses included both all-fraction combined and cellular-fraction inputs, exclusivity was recalculated using the all-fraction observations.

Reconciliation was validated for up to three size fractions which used a minimum two-fraction support rule by default For analyses with more than three fractions, users may manually adjust the rule in `reconcile_size_fractions.R`.

Optional: Config can be provided by the user. Default config has been validated for abundance prediction using GAM.

`min_support`: minimum number of actual detections (norm_coverage > 0) required to support a peak
`min_rel_height`: minimum relative height of a secondary peak compared with the tallest peak for the same taxon
`min_abs_fit`: minimum absolute fitted abundance value required for a peak to be considered real
`min_sep_deg`: minimum distance in latitudinal degrees used when deciding whether nearby peaks should be merged

Optional config file format (tab-separated file with columns parameter and value)
Example:
parameter	        value
`min_support`	    2
`min_rel_height`	0.5
`min_abs_fit`	    0.005
`min_sep_deg`	    5

Output files
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

Usage:

Rscript run_biogeography.R \
  --pred_input gam_predictions_ALL_first100.tsv \
  --obs_input observations_ALL.tsv \
  --output_dir output \
  --script_dir /path/to/scripts

Or provide own config file:

Rscript run_biogeography.R \
  --pred_input gam_predictions_ALL_first100.tsv \
  --obs_input observations_ALL.tsv \
  --output_dir output \
  --config_file config_file.tsv \
  --script_dir /path/to/scripts

Or run multiple fractions:
  --pred_input gam_predictions_ALL_first100.tsv,gam_predictions_CELLULAR_first100.tsv,gam_predictions_VIR_first100.tsv \
  --obs_input observations_ALL.tsv,observations_CELLULAR.tsv,observations_VIR.tsv \
  --reference_obs_input observations_ALL.tsv \
  --output_dir output \
  --script_dir /path/to/scripts

**Test data**
The eukaryotic all-fraction inputs in `stage_2/test_input/all_fractions_combined/` are the same as the Stage 1 input and output. A separate multiple-fraction dataset, including all fractions combined, cellular fraction, and viral fraction, is provided in `stage_2/test_input/size_fraction_reconciliation/` to test size fraction reconciliation. 