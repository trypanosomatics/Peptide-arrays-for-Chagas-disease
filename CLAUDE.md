# CLAUDE.md

R pipeline for high-density peptide microarray data from the Chagas Antigen and Epitope Atlas
(CHAGASTOPE; Ricci et al. 2023, Nat Commun, doi:10.1038/s41467-023-37522-9). It accompanies the
Methods in Molecular Biology chapter by Romer, Quinteros & Agüero (2025); `EXAMPLES.md` mirrors the
chapter's code boxes. Known issues and open questions are in `REVIEW_NOTES.md` (untracked); read it
before touching the steps it mentions.

## Two arrays, two data flows

| | Pools (CHAGASTOPE-v1) | Individual sera (CHAGASTOPE-v2) |
|---|---|---|
| Scripts | 01, 02, 03, 04, 08, 10 | 05, 06, 07, 09, 10 |
| Design file | `Supplementary File S08 - Mapping of CHAGASTOPE-v1 ...tsv` | `Supplementary File S09 - Mapping of CHAGASTOPE-v2 ...tsv` |
| Raw data | ArrayExpress E-MTAB-11651, `{SRC}_{PO\|NE}_raw.tsv` | E-MTAB-11655, `{SRC}_{P#\|E#}_PO_raw.tsv` |
| Sources | `AR BO BR CO MX US` (+ `LE`) | `AR_P1`…`US_E6` (country_serum) |
| Tiling | 16-mers, offset 4 (overlap 12), whole proteomes | 16-mers, offset 1, antigenic regions (+ alanine scans) |
| Raw columns | `{SRC}_{TYPE}_rawData_r1/_r2` | `allData.replica1/2` |

- Both arrays are defined once in `ARRAYS` in `code/functions/00_common.R`: peptide length/offset, types,
  default sources, design file and groups. Every script sources it and takes `--array v1|v2`, which
  defaults to the array the step is written for. `get_array(id, supported = ...)` stops on arrays a
  step can't handle (01–04, 08: v1; 05–07, 09: v2; 10: both).
- `PO`/`NE` = pools of Chagas-positive and negative sera. `LE` (Leishmaniasis pool) is normalized
  together with the negatives and excluded from global statistics (`negative_like_sources` in the
  config).
- Design `group` values that are used for filtering: `Trypanosoma_cruzi` (v1), `antigenic_regions_Tcruzi`
  and `alascan_best_peaks_Tcruzi` (v2).
- `data/extra_files/`: T. cruzi Sylvio X10 and CL Brener proteomes (TriTrypDB, Dec 2016, `protein`,
  `sequence`) used by 04; alanine-scan designs used by 07 (full and test subset).

## Pipeline

Every step is a thin CLI wrapper `code/NN_name.R` (parses args, builds paths, `source()`s its function
file) plus `code/functions/NN_name_function.R` (the logic). Outputs go to `data/<project>/outputs/NN_*/`.

**Output names** are `<array>.<step>.<fields>.<ext>`, dot-separated, built by `output_file_name()` and
parsed by `parse_output_file_name()` (fields per step in `FILE_STEPS`). Parameter labels are generated
from the arguments: `smoothing_label()` → `med3_mean0_zeros`, `cutoff_label()` → `cutoff4SD_2pep`
(decimals as `p`). Only one field per step may contain dots (protein IDs like `TcCLB.511671.60`).
Every tabular output starts with an `array` column, checked on read with `check_array()`. Never
hand-write output paths or match them with regexes: use `output_file_name()` / `find_output_files()`.
`BREAKING-CHANGES.md` has the old → new mapping.

| Step | Does | Writes |
|---|---|---|
| 01 | Quantile-normalizes (preprocessCore) all PO columns together and all NE+LE columns together; global `mode` and `sd` of normalized signal | `v1.normalized.<SRC>.<TYPE>.tsv`, `v1.global_statistics.tsv` |
| 02 | Maps peptides to proteins, rolling median then rolling mean per protein and replicate, mean/sd over the 2 replicates | `v1.smoothed.<SRC>.<TYPE>.<smoothing>.tsv` |
| 03 | Peak = run of ≥ `min_num_of_peptides_in_peak` peptides with `mean_smoothed_signal >= mode + k*sd`; also overlap peaks across source combinations (`source` = `"AR, BO"`) and NE/PO signal ratios | `v1.peaks.<smoothing>.<cutoff>.tsv` |
| 04 | PO single-source peaks with NE ratio ≤ 0.2, padded ±16 peptides, merged into antigenic regions | `v1.regions_{peptides,summary}.<smoothing>.<cutoff>.tsv` |
| 05 | Quantile normalization of individual sera, both whole array (`allData.*`) and region peptides only (`onlyRegionsData.*`, used downstream and for the stats) | `v2.normalized.<SRC>.PO.tsv`, `v2.global_statistics.tsv` |
| 06 | Same smoothing as 02, per protein and region | `v2.smoothed.<SRC>.PO.<smoothing>.tsv` |
| 07 | Alanine scan: signal change for each Ala-substituted peptide vs the wild type, per serum, computed from **raw** (not normalized) v2 data; heatmaps and sequence logos | `v2.alascan_{change,matrix,heatmap,logo}.<protein>...` |
| 08 | Pools profiles per protein, one panel per source, threshold line | `v1.profiles[_index].<smoothing>.<k>SD.<all\|protein>.{pdf,tsv}` |
| 09 | Individual-sera ridge plots per protein, faceted by country | `v2.ridgeplots[_index].<smoothing>.<all\|protein>.{pdf,tsv}` |
| 10 | Raw vs normalized vs smoothed signal distributions | `<array>.signal_distribution.<smoothing>.{positives,negatives}.pdf` |

## Running

```bash
cd code                      # required: main_folder defaults to dirname(getwd()); .Rprofile activates renv
./01_pools_normalize_data.R --testing TRUE
./02_pools_smooth_data.R --testing TRUE
./03_calculate_peaks.R --testing TRUE
./04_calculate_regions.R --testing TRUE
./05_individual_serums_normalize_data.R --testing TRUE
./06_individual_serums_smooth_data.R --testing TRUE
./07_alanine_scan_analysis.R --testing TRUE
./08_pools_plot_protein_profiles.R --testing TRUE
./09_individual_serums_plot_protein_profiles.R --testing TRUE
./10_plot_signal_variations.R --testing TRUE   # add --array v2 for individual sera
Rscript tests/test_file_names.R                # file-name builder/parser checks
```

- `--testing TRUE` uses `data/test_data` (subset of about 20 proteins, checked in). The full test run
  takes about 5 minutes; 08 and 09 are the slowest. Without it, scripts read `data/chagastope_data`, which
  you populate from ArrayExpress/figshare (see `EXAMPLES.md` §3.2).
- The test subset is enriched for antigens, so its global SD is inflated (cutoff about 35k at 4 SD).
  Peak calls on test data are a smoke test, not biology.
- The arg parser is hand-rolled: strict `--flag value` pairs, comma lists **without spaces**, unknown
  flags silently ignored, logicals via `as.logical` (`TRUE`/`FALSE`/`T`/`F`). Each wrapper prints the
  values it actually used. Check that echo.
- 02/06 print the smoothing label of their run. Later steps (03, 04, 08, 09, 10) select it with
  `--smoothing_label` (default `med3_mean0_zeros`). 04 selects peaks with the same
  `--sd_multiplier_for_cutoff` / `--min_num_of_peptides_in_peak` given to 03. Runs with different
  parameters coexist safely.
- Full data are RAM-heavy (the whole matrix is quantile-normalized in memory).

## Environment and dependencies

- R 4.3.3, Bioconductor 3.18, pinned with renv (`code/renv.lock`, `code/renv/settings.json`). The project
  library is `code/renv/library/` and the renv root is `code/`, not the repo root.
- Direct deps: data.table, zoo, preprocessCore (Bioc), dplyr, reshape2, tidyr, readr, ggplot2, patchwork,
  gridExtra, ggridges, ggh4x, scico, pheatmap, ggseqlogo. `00_install_requirements.R` is the non-renv
  fallback; keep its package list in sync.
- Upgrading: on a branch, start from a baseline test-run output, run `renv::update()` (or bump R/Bioc in
  settings), rerun all 10 steps with `--testing TRUE`, diff the TSVs against the baseline (below), look at
  the PDFs, then `renv::snapshot()` and update the README "tested with" line.
- Upgrade hazards: 07 contains about 600 lines of **vendored, modified ggseqlogo internals** (it calls
  `get_font()`, which reads ggseqlogo's `extdata` font files); `aes_string()` (already warns) and
  `guides(fill = F)` are deprecated; 07's heatmap recolouring indexes pheatmap grobs by position.

## Conventions (match them)

- Wrapper/function split, numbered pairs, `#### SECTION ####` headers, `snake_case` main functions and
  some camelCase helpers. New pipeline steps get the next number, a wrapper with the same arg-parsing
  block and echo, and an output dir with an `OUTPUT WILL BE GENERATED HERE` placeholder in **both**
  `data/test_data/outputs/` and `data/chagastope_data/outputs/`.
- data.table everywhere except 07/10 (dplyr/tidyr). Read with `fread(..., sep = "\t")`. Write with
  `write.table(..., sep = "\t", quote = TRUE, row.names = FALSE)`.
- Key column names downstream code relies on: `Reporter Name` (= design `array_express_id`), `protein`,
  `start`, `source`, `type`, `mean_smoothed_signal`, `sd_smoothed_signal`, `sequence`, `truncated`,
  `region` (v2).
- Peptide geometry is **array-specific**: v1 = 16-mers at offset 4, v2 = 16-mers at offset 1. Any code
  that walks along a protein (tiling, joining peptides into a sequence, borders) must take
  `peptide_length` / `peptide_offset` from `get_array()`, never assume one. A new array design means a
  new `ARRAYS` entry, not new literals. Still hardcoded: ±1300 signal thresholds in 07. Parametrize
  them when you touch that code.
- Main functions take an `array_config` (from `get_array()`) plus explicit paths and parameters. Don't
  rely on variables defined in the wrapper; earlier code did, and it made functions unusable outside
  their script.
- Inside `[.data.table`, a function argument with the same name as a column (`protein`, `source`,
  `type`, `start`...) is shadowed by the column. Copy it to another name before filtering
  (see `prepare_plot_data()` in 08/09).
- 08 and 09 both define `prepare_plot_data()` (identical today). Don't source both in one session.
- Plots: colour-blind-safe palettes already in use (Okabe-Ito in 07, scico `batlow`/`lipari` in 09/10).
  Keep to those, write PDFs, and add an index TSV when a PDF has one page per protein.

## Where new work goes

- **Pipeline changes / new steps**: `code/` + `code/functions/` as above. Register new output kinds in
  `FILE_STEPS` and add them to `tests/test_file_names.R`. Shared helpers belong in
  `functions/00_common.R`. `calculateMode` and `smoothVector` are still copy-pasted in 01/05 and 02/06,
  so move them there rather than adding a third copy.
- **One-off questions / CLI tools**: `code/tools/` (create it), reading pipeline **outputs** through
  `source("functions/00_common.R")` + `find_output_files()`, never re-implementing normalization or
  smoothing. Small `Rscript` with `commandArgs` or a one-liner (below).
- **Manuscript analyses**: `analysis/<short-name>/` with its own README stating inputs (which outputs,
  which output names and labels), so figures can be regenerated.

## Data and git hygiene

- `.gitignore` ignores `*.tsv *.csv *.pdf`. The test inputs are force-tracked. Use `git add -f` only for
  deliberate small fixtures, and never commit full ArrayExpress data or pipeline outputs.
- Don't commit or push unless asked. `EXAMPLES.md` mirrors the published chapter (code at `60f44bf`) and
  stays as printed. Keep `README.md` and `BREAKING-CHANGES.md` accurate when you change flags or output
  names, and add a new section to `BREAKING-CHANGES.md` for each future breaking change.

## Verifying changes

`tests/test_file_names.R` covers only the name builder/parser. For everything else, rerun the test
pipeline and compare outputs with a baseline run made **before** the change (copy
`data/test_data/outputs` aside first):

```bash
Rscript -e 'a <- data.table::fread(commandArgs(TRUE)[1]); b <- data.table::fread(commandArgs(TRUE)[2]); print(all.equal(a, b))' base/03_pools_antigenic_peaks/v1.peaks.med3_mean0_zeros.cutoff4SD_2pep.tsv new/03_pools_antigenic_peaks/v1.peaks.med3_mean0_zeros.cutoff4SD_2pep.tsv
```

For bug fixes that intentionally change output, say which columns change and why.

## One-liner cookbook (run from `code/`, test-data paths)

```bash
O=../data/test_data/outputs
# strongest single-source PO peaks
Rscript -e 'library(data.table); p <- fread("'$O'/03_pools_antigenic_peaks/v1.peaks.med3_mean0_zeros.cutoff4SD_2pep.tsv"); print(p[type == "PO" & !grepl(",", source)][order(-best_peptide_signal)][1:10, .(source, protein, peak_start, peak_last_end, best_peptide, best_peptide_signal)])'
# proteins by number of countries sharing an overlapping PO peak
Rscript -e 'library(data.table); p <- fread("'$O'/03_pools_antigenic_peaks/v1.peaks.med3_mean0_zeros.cutoff4SD_2pep.tsv"); p[, n_src := lengths(strsplit(source, ", "))]; print(p[type == "PO", .(max_countries = max(n_src)), by = protein][order(-max_countries)])'
# fraction of region peptides above mode + 4 SD per individual serum
Rscript -e 'library(data.table); source("functions/00_common.R"); st <- fread(output_file_name("'$O'/05_individual_serums_normalized_data", "v2", "global_statistics")); cut <- st$mode + 4 * st$sd; d <- rbindlist(lapply(find_output_files("'$O'/06_individual_serums_smoothed_data", "v2", "smoothed", smoothing = "med3_mean0_zeros"), fread)); print(d[, .(frac_above = round(mean(mean_smoothed_signal >= cut), 3)), by = source][order(-frac_above)])'
# top peptides of one protein in one pool
Rscript -e 'library(data.table); d <- fread("'$O'/02_pools_smoothed_data/v1.smoothed.AR.PO.med3_mean0_zeros.tsv"); print(d[protein == "TcCLB.511671.60"][order(-mean_smoothed_signal)][1:5, .(start, sequence, mean_smoothed_signal)])'
```
