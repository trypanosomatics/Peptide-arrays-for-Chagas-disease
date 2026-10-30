# Review notes (untracked, for discussion)

Code review of the pipeline as of commit `9bb97c6` (2026-10-30). I ran the whole test pipeline
(`--testing TRUE`, steps 01–10) in a scratch copy. All ten steps exit 0 with default arguments, in about
4.5 minutes. Items marked **[confirmed]** were reproduced; **[by reading]** ones were not executed.

## Status on branch `array-config-filenames` (commits 37fd8ea, c6c0f30, 196c2e0)

**Fixed** (see `BREAKING-CHANGES.md`):
- §1.1 region sequences. They now append `peptide_length - offset + 1 .. peptide_length` from the
  array config, and 21/21 test regions equal the protein substring.
- §1.2 run mixing in 08/09 and §1.3 peaks overwriting. Smoothing/cutoff labels are now in the names,
  with exact matching through `find_output_files()`.
- §1.4 step 10 with v2 (`--array v2`).
- §1.5 multi-serum logos, §1.6 multi-protein selection, and the `next()` part of §1.7.
- §6 is implemented: array config, `--array`, dot-separated names, `array` column.
- From §3: wrapper-global dependencies are gone, and 08/09 `prepare_plot_data()` are now identical.
  09 parses `--only_proteins_*`, 08 `--fixed_scale` is a logical, and the undefined `plot_data` in 10
  is gone. 03 uses `copy()`, detects combined peaks with `grepl(",")`, and the dead line is removed.
- New finding, fixed: in step 10 for v2, `all_data_smoothed[,"source"]` on a tibble turned the whole
  column into a single deparsed string, so the "Smoothed" panel had one bogus category. It's now
  grouped by country like the raw and normalized panels.

Normalized, smoothed, statistics, peaks, region-peptide and index tables are identical to the old
outputs (176/176 compared, ignoring the new `array` column).

**Still open:** the rest of §1.7 (heatmap breaks, single-column drop, `make_bold_names`), all of §2,
duplicated `calculateMode`/`smoothVector`, the arg parser, vendored ggseqlogo, `http://` in
`00_install_requirements.R`, §4 doc items, and v2 support in 03/04 (per-region peaks).

## 1. Bugs affecting results or common invocations

1. **Antigenic region sequences are truncated** [confirmed]:
   `code/functions/04_calculate_regions_function.R:143` and `:168`.
   The region sequence is built as first peptide + `substring(peptide, 16, 16)` of each following
   peptide, which is only right for a 1-residue offset. The offset depends on the array:
   CHAGASTOPE-v1 (pools) tiles at offset 4, CHAGASTOPE-v2 (individual sera) at offset 1. The code
   uses `chronic_peaks_peptide_offset` (set to 4 in `04_calculate_regions.R:329`) everywhere else:
   peptide tiling, borders, region merging. Only the sequence assembly ignores it.
   On the v1 test data **21/21** regions have `nchar(sequence) != aa_amount` (e.g. region 1: 49 vs
   148 aa; no-border: 17 vs 20). With offset 1 the current code would happen to be correct.
   Fix: use the offset parameter instead of the literal,
   `substring(peptide, peptide_length - offset + 1, peptide_length)` (as 03 does at
   `03_calculate_peaks_function.R:202`). With offset 1 this reduces to the current `substring(peptide, 16, 16)`, so
   one code path serves both arrays. The array should come from a global config, not a per-script
   literal (see §6). Today 04 only reads v1 pool peaks from 03, so a v2 run also needs a v2 peaks
   input. `start`/`last_start`/`aa_amount` are
   correct; only the `sequence` and `no_border_sequence` columns are wrong for v1. **Check whether
   any published table was derived from these columns.**

2. **08/09 silently mix smoothing runs** [confirmed]:
   `08_pools_plot_protein_profiles_function.R:28-29`, `09_…:257-258`.
   `--profile_data_suffix` is used as a regex (default `smoothed.tsv`), so it also matches
   `m3_smoothed.tsv`. With two runs present, every (source, type, start) appears twice, and profiles,
   error bars and protein ordering combine both runs. Sources are matched with `grepl` against the
   **full path** (`"CO"` matches any path containing CO). Fix: build exact filenames
   `sprintf("%s_%s_%ssmoothed.tsv", source, type, suffix)` as 03 does.

3. **03 output name ignores `--profile_data_suffix`** [by reading]: `03_calculate_peaks_function.R:429`.
   Peaks from different smoothing settings overwrite each other, and 04 cannot tell them apart.

4. **10 crashes for individual sera with default args** [confirmed]: `10_plot_signal_variations.R:916`
   default `_smoothed\.tsv$` doesn't match 06's `*_smoothed_signals.tsv`, so the table is empty and
   `undefined columns selected` is raised at `functions/10_…:97`. It works with
   `--smoothed_data_suffix smoothed_signals.tsv`. Fix: choose the default from `--type_data`.

5. **07 with several `--sequence_logo_source` values crashes** [confirmed]:
   `07_…:116` `if (sources != "all")` → "the condition has length > 1". The PDF filename at `:1077` would
   also vectorize. Fix: `if (!identical(sources, "all"))` and `paste(sources, collapse = "-")`.

6. **07 multi-protein selection** [by reading]: `07_…:1039` uses `== selected_protein`, but the wrapper
   splits `--selected_protein` on commas, so several proteins get recycled and silently filtered. Use
   `%in%`.

7. **07 edge cases** [by reading]: `:28` `next()` outside a loop (error if a protein has <2 non-NA
   signals). `:72-73` heatmap breaks are non-monotonic when all changes have one sign (pheatmap
   error). `:69` drops to a vector when one column remains. `:82` calls undefined `make_bold_names`
   (unreachable now).

## 2. Methodological points worth stating explicitly (manuscript / chapter)

- The peak cutoff `mode + k·sd` is computed on **unsmoothed normalized** signals (01/05) but applied to
  **smoothed** profiles (03/08). The SD is taken over all peptides, reactive ones included, so it depends
  on array content. On the antigen-enriched test subset, sd is about 8,700 vs mode about 160 (cutoff
  about 35k at 4 SD). Consider a robust scale (MAD) or computing the SD on NE/LE data.
- 07 computes alanine-scan effects on **raw**, per-array signals with absolute ±1300 thresholds
  ("fundamental"/"suboptimal" residues, column-clustering filter). Thresholds are scanner/array
  specific. Variables are named `Chips_roche_norm` although data aren't normalized.
- 05 statistics come from the region-only normalization. 09 puts `sd_multiplier` in the filename but
  draws no threshold.
- Border handling `"zeros"` pulls signal down at protein termini (comment acknowledges it). Peaks at the
  N/C-terminus are underestimated with median windows ≥3.
- Defaults are median 3 / mean 0 in code. The EXAMPLES commands use other combinations. State which
  were used for any figure.

## 3. Maintainability (relevant to the planned extensions)

- Function files depend on **wrapper globals** (03: `profile_data_folder`, `sequence_length`,
  `sequence_overlap`, decimals; 05: `design_data_file`; 06: `design_data_file`, `design_groups`,
  `types`, `normalized_data_folder`; 07: `design_alanine_file`, `project_folder`). They can't be reused
  from new scripts without recreating that state. Pass everything explicitly.
- Copy-pasted helpers: `calculateMode` (01, 05), `smoothVector` (02, 06), `prepare_plot_data` (08, 09,
  **same name, different behaviour**). Move them to a shared `functions/utils.R`.
- Hand-rolled arg parser in each of 10 wrappers: no validation, unknown flags ignored. 09 ignores
  `--output_suffix`, `--only_proteins_above/below` (used in EXAMPLES). 08 `--fixed_scale` is not
  coerced to logical. 05 `testing` arg unused. 10 wrapper passes an undefined `plot_data` (harmless due
  to lazy evaluation). Consider `optparse`, which is small and already common in Bioc workflows.
- `setwd(main_folder)` in 01/02/10. `main_folder <- dirname(getwd())` forces running from `code/`.
- 03: `rbindlist` inside nested loops (quadratic, slow on full data). `:142-143` uses a "copy hack"
  instead of `copy()`. `:372` is a dead line referencing a stale variable. `:419` `nchar(source) > 2`
  detects combined peaks, whereas 04 uses `grepl(",")`. Use one predicate.
- 02: smoothed values are assigned back positionally (`:85-86`). This is correct only because of the
  preceding sort. Use `:=` by group instead.
- 07 vendors about 600 lines of modified ggseqlogo (`get_font`, `geom_logo`, ...). Check license
  attribution for the copied code, and expect breakage on ggseqlogo upgrades. Upstream `geom_logo()`
  with `method = "custom"` plus a per-position fill might replace most of it.
- Deprecations already visible: `aes_string()` (warns in 07), `guides(fill = F)`.
- `00_install_requirements.R` uses `http://cran.rstudio.com/`. Switch to https (and/or drop the script
  in favour of renv).
- No tests or CI. A cheap first step is a script that runs the test pipeline and compares key TSVs with
  stored checksums or `all.equal` against a committed baseline (small).

## 4. Docs and metadata

- README: heading "Instructions for Ubuntu Linux 22.04" vs text about 24.04. The DOI link
  `https://10.1038/s41467-023-37522-9` is broken and should be `https://doi.org/10.1038/...`. Code comments
  in 05/06 still refer to old input folder names (`11_…`, `12_…`) and a `-test F` flag that doesn't exist.
- EXAMPLES.md:
  - §3 generic example has an en-dash `–-parameter2`.
  - §3.5: `--smoothing_mean_window_size 5 --smoothing_median_window_size 1` labelled
    `mean_1_median_5` (inverted).
  - §3.11 individual: `--sources AR_P1, AR_P2, ...` with spaces breaks the pair parser. Its
    `--output_suffix` is not supported by 09.
  - Windows §3.9 is missing `--` before `smoothing_median_window_size`.
  - §3.10 `--selected_protein "TCSYLVIO_004530"` with `--testing FALSE` is fine. The test subset only
    contains TcCLB.511671.60.
- `code/.RData` and `code/.Rhistory` are untracked. Add them to `.gitignore`.
- `.gitignore` excludes `*.tsv`, but all test inputs are tracked (force-added). That's fine, but worth a
  comment in `.gitignore` so new fixtures aren't silently skipped.

## 5. Suggestions for the next round of work

1. Fix §1.1–1.5 (small, local changes) and regenerate regions if any downstream use exists.
2. Shared `utils.R` + explicit function arguments, which unblocks `code/tools/` one-liners and reuse.
3. Baseline regression check on test data before the dependency upgrade (R 4.4/4.5, Bioc 3.20+).
4. Tidy long-format export (one table: source, type, protein, start, signal, sd, above_cutoff) for
   downstream stats and plotting outside the pipeline.
5. Plots: per-protein interactive or summary views (peak overlap across countries as an UpSet plot,
   signal heatmap of regions × sera, seroprevalence per region with a stated cutoff).

## 6. Proposal: global array config + array tag in output filenames

Goal: every step knows which array design (v1 or v2) it is working on, from one place, and every
output says which array produced it.

**Config.** One file, e.g. `code/functions/array_config.R`, sourced by every wrapper:

```r
ARRAYS <- list(
  v1 = list(peptide_length = 16, peptide_offset = 4, types = c("PO", "NE"),
            design_file = "inputs/01_pools_array_design/Supplementary File S08 - Mapping of CHAGASTOPE-v1 data to T cruzi proteins.tsv",
            design_group = "Trypanosoma_cruzi"),
  v2 = list(peptide_length = 16, peptide_offset = 1, types = "PO",
            design_file = "inputs/03_individual_serums_array_design/Supplementary File S09 - Mapping of CHAGASTOPE-v2 data to T cruzi proteins.tsv",
            design_group = "antigenic_regions_Tcruzi")
)
get_array <- function(id) { if (!id %in% names(ARRAYS)) stop("unknown array: ", id); c(id = id, ARRAYS[[id]]) }
```

Wrappers get `--array v1|v2`. Each defaults to its natural array (01–04 and 08: v1; 05–07 and 09: v2),
so existing commands in the chapter/EXAMPLES keep working. Everything currently hardcoded (16/12/4 in
03/04, `substring(peptide, 16, 16)`, design paths and groups, `types`) is read from `get_array()`
and passed to functions as arguments. This also removes most of the wrapper-global dependencies in §3.

**Filenames.** Tag outputs with the array, e.g. `v1_AR_PO_smoothed.tsv`,
`v1_pools_peaks_cutoff4SD_2pep.tsv`, `v1_4SD_2pep_summary_region_data.tsv`,
`v2_AR_P1_PO_smoothed_signals.tsv`. Also add an `array` column inside each output table, so the
information survives renames and `rbind`s. Readers check that the tag/column matches the `--array`
they were given and stop on a mismatch.

Points to decide:
- **Prefix vs suffix.** A prefix makes `v1_*` globs trivial, but 07 and 10 parse sample/source names by
  splitting filenames on `_`, and 08/09 match sources by regex. A shared filename builder/parser in
  `utils.R` is needed either way. This also fixes §1.2 (regex over filenames).
- **Backward compatibility.** Renamed outputs break the chapter's printed commands and any existing
  output folders. Option: write the tagged name and accept the old untagged name on read (treated as
  the step's default array, with a warning) for one release.
- 03/08/10 statistics: `global_statistics.tsv` exists per array already (01 vs 05 folders). Tag it
  too, so a peak run can't accidentally pick up the other array's cutoff.
