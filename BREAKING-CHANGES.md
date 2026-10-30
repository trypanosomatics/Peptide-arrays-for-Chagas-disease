# Breaking changes

## Array-aware pipeline and new output file names

**Book chapter readers:** the commands printed in the chapter (and in `EXAMPLES.md`) belong to the
code *before* this change. To reproduce the chapter exactly, check out commit `60f44bf`:

```bash
git checkout 60f44bf
```

### Why

- Outputs did not say which array design produced them. CHAGASTOPE-v1 (pooled sera) tiles proteins
  with 16-mers at a 4-residue offset, and CHAGASTOPE-v2 (individual sera) tiles with 16-mers at a
  1-residue offset. Code that walks along a protein needs the right offset.
- Parameters (smoothing windows, cutoffs) were either missing from file names or given as free-text
  suffixes. This let runs overwrite each other (03) or get mixed silently (08/09 matched file names with
  regular expressions).
- Peptide length/offset, design files, sample types and default sources were hardcoded in several
  scripts.

### What changed

1. **One array config.** `code/functions/00_common.R` defines `ARRAYS` (`v1`, `v2`) and is sourced by
   every script. Each script has a new `--array` flag that defaults to the array it handles. Steps
   01–04 and 08 support v1, steps 05–07 and 09 support v2, and step 10 supports both. Any other value
   stops with an error.
2. **New output file names**, built and parsed by one function pair (`output_file_name()` /
   `parse_output_file_name()`):

   ```
   <array>.<step>.<field>...<ext>        e.g.  v1.smoothed.AR.PO.med3_mean0_zeros.tsv
   ```

   Fields are separated by dots. Parameter labels are built from the actual arguments:
   `med<median>_mean<mean>_<borders>` for smoothing and `cutoff<k>SD_<n>pep` for peaks. Decimals are
   written with `p` (4.5 SD → `cutoff4p5SD_2pep`). Protein IDs keep their dots
   (`v2.alascan_matrix.TcCLB.511671.60.1.tsv`), because the parser reads the other fields from both
   ends of the name.
3. **`array` column.** Every tabular output now starts with an `array` column (`v1`/`v2`). Readers check
   it and stop if it is missing or doesn't match. The alanine-scan matrices are the exception: they
   stay plain matrices with serum row names.
4. **Bug fix in 04, results change.** Antigenic region sequences (`sequence`, `no_border_sequence` in
   the region summary) were built by appending one residue per peptide. With v1's 4-residue offset,
   every region sequence was truncated (on the test data, 21/21 regions, e.g. 49 instead of 148 aa).
   Region sequences now equal the protein substring of the stated length. Coordinates (`start`,
   `last_start`, `aa_amount`) were already correct. **Discard region tables made with older code.**
5. Other behaviour fixes:
   - 07 accepts several `--selected_protein` and `--sequence_logo_source` values, and skips proteins
     without enough data instead of erroring.
   - 10 works for individual sera with default arguments.
   - 10 groups v2 smoothed data by country, consistent with raw and normalized data.
   - 07's long-format change table no longer has an unnamed row-number column.

All other values are unchanged. On the test data, every normalized, smoothed, statistics, peaks,
region-peptide and plot-index table is identical to the previous version, apart from the new `array`
column.

### Output file names, old → new

Paths are relative to `data/<project>/outputs/`. `<SRC>` is a source (`AR`, `AR_P1`, ...),
`<TYPE>` is `PO`/`NE`, `<smoothing>` is e.g. `med3_mean0_zeros`, `<cutoff>` is e.g.
`cutoff4SD_2pep`, `<scope>` is `all` or the protein given with `--protein`.

| Step | Old name | New name |
|---|---|---|
| 01 | `01_pools_normalized_data/<SRC>_<TYPE>_processed.tsv` | `01_pools_normalized_data/v1.normalized.<SRC>.<TYPE>.tsv` |
| 01 | `01_pools_normalized_data/global_statistics.tsv` | `01_pools_normalized_data/v1.global_statistics.tsv` |
| 02 | `02_pools_smoothed_data/<SRC>_<TYPE>_[<suffix>_]smoothed.tsv` | `02_pools_smoothed_data/v1.smoothed.<SRC>.<TYPE>.<smoothing>.tsv` |
| 03 | `03_pools_antigenic_peaks/pools_peaks_cutoff<k>SD_<n>pep.tsv` | `03_pools_antigenic_peaks/v1.peaks.<smoothing>.<cutoff>.tsv` |
| 04 | `04_antigenic_regions/<k>SD_<n>pep_peptide_region_data.tsv` | `04_antigenic_regions/v1.regions_peptides.<smoothing>.<cutoff>.tsv` |
| 04 | `04_antigenic_regions/<k>SD_<n>pep_summary_region_data.tsv` | `04_antigenic_regions/v1.regions_summary.<smoothing>.<cutoff>.tsv` |
| 05 | `05_individual_serums_normalized_data/<SRC>_PO_processed.tsv` | `05_individual_serums_normalized_data/v2.normalized.<SRC>.PO.tsv` |
| 05 | `05_individual_serums_normalized_data/global_statistics.tsv` | `05_individual_serums_normalized_data/v2.global_statistics.tsv` |
| 06 | `06_individual_serums_smoothed_data/<SRC>_PO_[<suffix>_]smoothed_signals.tsv` | `06_individual_serums_smoothed_data/v2.smoothed.<SRC>.PO.<smoothing>.tsv` |
| 07 | `07_alanine_scan_raw_data/Raw_data_signal_change_alanine_scan_<protein>.tsv` | `07_alanine_scan_raw_data/v2.alascan_change.<protein>.tsv` |
| 07 | `…/Raw_data_signal_change_matrix_alanine_scan_<protein>_group_<i>.tsv` | `…/v2.alascan_matrix.<protein>.<i>.tsv` |
| 07 | `…/Heatmap_alanine_scan_<protein>_group_<i>.pdf` | `…/v2.alascan_heatmap.<protein>.<i>.pdf` |
| 07 | `…/sequence_logo_alanine_scan_<protein>_group_<i>_<sera>.pdf` | `…/v2.alascan_logo.<protein>.<i>.<sera>.pdf` (`all` or `AR_P1+BO_P1`) |
| 08 | `08_09_10_plots/protein_profiles_cruzi_<k>SD[_selected_proteins][_<suffix>].pdf` | `08_09_10_plots/v1.profiles.<smoothing>.<k>SD.<scope>.pdf` |
| 08 | `…/protein_profiles_cruzi_<k>SD[…]_index.tsv` | `…/v1.profiles_index.<smoothing>.<k>SD.<scope>.tsv` |
| 09 | `08_09_10_plots/individual_serums_protein_ridgeplots_<k>SD[…].pdf` | `08_09_10_plots/v2.ridgeplots.<smoothing>.<scope>.pdf` |
| 09 | `…/individual_serums_protein_profiles_cruzi_<k>SD[…]_index.tsv` | `…/v2.ridgeplots_index.<smoothing>.<scope>.tsv` |
| 10 | `08_09_10_plots/Signal_{positives,negatives}_{pools,individual_serums}.pdf` | `08_09_10_plots/<array>.signal_distribution.<smoothing>.{positives,negatives}.pdf` |

Input files (raw data from ArrayExpress, design files from figshare, `data/extra_files/`) keep their
names.

### Command-line flags

All scripts gain `--array`. Removed and renamed flags:

| Script | Removed | Use instead |
|---|---|---|
| 02, 06 | `--output_suffix` | nothing: the label is built from `--smoothing_median_window_size`, `--smoothing_mean_window_size`, `--smooth_borders_option` and printed as "Smoothing label" |
| 03 | `--profile_data_suffix` | `--smoothing_label` (default `med3_mean0_zeros`) |
| 04 | `--peaks_tag` | `--smoothing_label`, `--sd_multiplier_for_cutoff`, `--min_num_of_peptides_in_peak` (same values as given to 03) |
| 08 | `--profile_data_suffix`, `--output_suffix` | `--smoothing_label`; `--fixed_scale` is now parsed as a logical |
| 09 | `--profile_data_suffix`, `--sd_multiplier_for_cutoff` (it was only used in the file name) | `--smoothing_label`; `--only_proteins_above/below` now work from the command line |
| 10 | `--type_data`, `--raw_data_suffix`, `--normalized_data_suffix`, `--smoothed_data_suffix` | `--array v1\|v2`, `--smoothing_label` |

Examples:

```bash
# old
./02_pools_smooth_data.R --smoothing_median_window_size 3 --smoothing_mean_window_size 5 --output_suffix "mean_5_median_3"
./03_calculate_peaks.R --profile_data_suffix "mean_5_median_3" --sd_multiplier_for_cutoff 5
./04_calculate_regions.R --peaks_tag "5SD_2pep"
./08_pools_plot_protein_profiles.R --profile_data_suffix "mean_5_median_3"
./10_plot_signal_variations.R --type_data "individual_serums" --smoothed_data_suffix "mean_5_median_3_smoothed_signals.tsv"

# new
./02_pools_smooth_data.R --smoothing_median_window_size 3 --smoothing_mean_window_size 5   # prints: med3_mean5_zeros
./03_calculate_peaks.R --smoothing_label med3_mean5_zeros --sd_multiplier_for_cutoff 5
./04_calculate_regions.R --smoothing_label med3_mean5_zeros --sd_multiplier_for_cutoff 5
./08_pools_plot_protein_profiles.R --smoothing_label med3_mean5_zeros
./10_plot_signal_variations.R --array v1 --smoothing_label med3_mean5_zeros
```

### Migration

Outputs from older versions cannot be read by the new code. The new names aren't found, and a renamed
old file is rejected because it has no `array` column. Delete old outputs (keep the
`OUTPUT WILL BE GENERATED HERE` placeholders) and rerun steps 01–10:

```bash
find data/chagastope_data/outputs -type f ! -name 'OUTPUT WILL BE GENERATED HERE' -delete
```

For code of your own that reads pipeline outputs, use `output_file_name()` / `find_output_files()`
from `code/functions/00_common.R` instead of hand-written paths.

### Not changed yet

Steps 03 and 04 (peaks and regions) only support v1. Running them on v2 needs peak detection within
each v2 region, not along the whole protein.
