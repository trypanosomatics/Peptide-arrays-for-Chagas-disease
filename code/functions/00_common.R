#### SHARED CONFIG AND FILE NAMES ####
## Sourced by every numbered script in code/. It holds:
## 1. ARRAYS: the design of each CHAGASTOPE array (peptide geometry, design file, sources...)
## 2. The builder/parser for output file names, so every output says which array and
##    which parameters produced it.
##
## Output file names are dot-separated fields:
##   <array>.<step>.<field 1>.<field 2>...<ext>
##   e.g. v1.smoothed.AR.PO.med3_mean0_zeros.tsv
## The fields of each step are listed in FILE_STEPS. Field values cannot contain dots, except
## for at most one field per step (protein IDs such as TcCLB.511671.60), which the parser
## recovers by reading the other fields from both ends of the name.

#### ARRAY DESIGNS ####
ARRAYS <- list(
  # CHAGASTOPE-v1: whole proteomes, used with pooled sera (ArrayExpress E-MTAB-11651)
  v1 = list(
    name = "CHAGASTOPE-v1",
    peptide_length = 16,
    peptide_offset = 4,
    types = c("PO", "NE"),
    default_sources = c("AR", "BO", "BR", "CO", "MX", "US", "LE"),
    # LE (leishmaniasis pool) is normalized with the negatives and excluded from global statistics
    negative_like_sources = c("LE"),
    design_file = "inputs/01_pools_array_design/Supplementary File S08 - Mapping of CHAGASTOPE-v1 data to T cruzi proteins.tsv",
    design_group = "Trypanosoma_cruzi"
  ),
  # CHAGASTOPE-v2: antigenic regions at single-residue offset plus alanine scans, used with
  # individual sera (ArrayExpress E-MTAB-11655)
  v2 = list(
    name = "CHAGASTOPE-v2",
    peptide_length = 16,
    peptide_offset = 1,
    types = c("PO"),
    default_sources = c("AR_P1", "AR_P2", "AR_P3", "AR_P4", "AR_P5", "AR_P6", "BO_P1", "BO_P2", "BO_P3", "BO_P4", "BO_P5", "BO_P6", "BR_P1", "BR_P2", "BR_P3", "BR_P4", "BR_P5", "CO_P1", "CO_P2",
                        "CO_P3", "CO_P4", "MX_P1", "MX_P2", "MX_P3", "MX_P4", "MX_P5", "MX_P6", "US_P1", "US_P2", "US_P3", "US_P4", "US_P5", "US_P6", "AR_E1", "AR_E2", "AR_E3", "AR_E4", "AR_E5",
                        "AR_E6", "BO_E1", "BO_E2", "BO_E3", "BO_E4", "BO_E5", "BO_E6", "BR_E1", "BR_E2", "BR_E3", "BR_E4", "BR_E5", "BR_E6", "BR_E7", "CO_E1", "CO_E2", "CO_E3", "CO_E4", "CO_E5",
                        "CO_E6", "CO_E7", "MX_E1", "MX_E2", "MX_E3", "MX_E4", "MX_E5", "MX_E6", "US_E1", "US_E2", "US_E3", "US_E4", "US_E5", "US_E6"),
    negative_like_sources = character(0),
    design_file = "inputs/03_individual_serums_array_design/Supplementary File S09 - Mapping of CHAGASTOPE-v2 data to T cruzi proteins.tsv",
    design_group = "antigenic_regions_Tcruzi",
    alascan_group = "alascan_best_peaks_Tcruzi"
  )
)

# Returns the config of one array, plus its id. Stops if the array is unknown or if the
# calling step does not support it.
get_array <- function(array_id, supported = names(ARRAYS)) {
  if (!(array_id %in% names(ARRAYS))) {
    stop(sprintf("Unknown array '%s'. Valid arrays: %s", array_id, paste(names(ARRAYS), collapse = ", ")))
  }
  if (!(array_id %in% supported)) {
    stop(sprintf("This step does not support array '%s' (supported: %s)", array_id, paste(supported, collapse = ", ")))
  }
  c(list(id = array_id), ARRAYS[[array_id]])
}

#### PARAMETER LABELS ####
# Numbers in file names: no dots, decimals written with "p" (4.5 -> "4p5")
fmt_num <- function(x) {
  gsub(".", "p", format(x, scientific = FALSE, trim = TRUE, drop0trailing = TRUE), fixed = TRUE)
}

smoothing_label <- function(median_window_size, mean_window_size, borders) {
  sprintf("med%s_mean%s_%s", fmt_num(median_window_size), fmt_num(mean_window_size), borders)
}

cutoff_label <- function(sd_multiplier, min_num_of_peptides) {
  sprintf("cutoff%sSD_%spep", fmt_num(sd_multiplier), fmt_num(min_num_of_peptides))
}

#### OUTPUT FILE NAMES ####
# Fields of each step, in file-name order. dotted = the one field allowed to contain dots.
FILE_STEPS <- list(
  normalized          = list(fields = c("source", "type")),
  global_statistics   = list(fields = character(0)),
  smoothed            = list(fields = c("source", "type", "smoothing")),
  peaks               = list(fields = c("smoothing", "cutoff")),
  regions_peptides    = list(fields = c("smoothing", "cutoff")),
  regions_summary     = list(fields = c("smoothing", "cutoff")),
  alascan_change      = list(fields = c("protein"), dotted = "protein"),
  alascan_matrix      = list(fields = c("protein", "group"), dotted = "protein"),
  alascan_heatmap     = list(fields = c("protein", "group"), dotted = "protein"),
  alascan_logo        = list(fields = c("protein", "group", "sera"), dotted = "protein"),
  profiles            = list(fields = c("smoothing", "sd", "scope"), dotted = "scope"),
  profiles_index      = list(fields = c("smoothing", "sd", "scope"), dotted = "scope"),
  ridgeplots          = list(fields = c("smoothing", "scope"), dotted = "scope"),
  ridgeplots_index    = list(fields = c("smoothing", "scope"), dotted = "scope"),
  signal_distribution = list(fields = c("smoothing", "sign"))
)

get_file_step <- function(step) {
  if (!(step %in% names(FILE_STEPS))) {
    stop(sprintf("Unknown output step '%s'", step))
  }
  FILE_STEPS[[step]]
}

# Builds <folder>/<array>.<step>.<fields...>.<ext>. Fields are passed by name, e.g.
# output_file_name(folder, "v1", "smoothed", source = "AR", type = "PO", smoothing = "med3_mean0_zeros")
output_file_name <- function(folder, array, step, ..., ext = "tsv") {
  step_def <- get_file_step(step)
  values <- list(...)
  missing_fields <- setdiff(step_def$fields, names(values))
  extra_fields <- setdiff(names(values), step_def$fields)
  if (length(missing_fields) > 0 || length(extra_fields) > 0) {
    stop(sprintf("Output '%s' needs fields [%s], got [%s]", step,
                 paste(step_def$fields, collapse = ", "), paste(names(values), collapse = ", ")))
  }
  values <- vapply(step_def$fields, function(field) as.character(values[[field]]), character(1))

  parts <- c(array = array, values)
  for (field in names(parts)) {
    value <- parts[[field]]
    if (is.na(value) || value == "") {
      stop(sprintf("Output '%s': field '%s' is empty", step, field))
    }
    if (grepl("[/\\\\]", value)) {
      stop(sprintf("Output '%s': field '%s' contains a path separator: %s", step, field, value))
    }
    if (grepl(".", value, fixed = TRUE) && !identical(field, step_def$dotted)) {
      stop(sprintf("Output '%s': field '%s' cannot contain dots: %s", step, field, value))
    }
  }

  file_name <- paste(c(array, step, values, ext), collapse = ".")
  if (is.null(folder) || folder == "") file_name else file.path(folder, file_name)
}

# Inverse of output_file_name(). Returns NULL for names that don't follow the scheme
# (e.g. the "OUTPUT WILL BE GENERATED HERE" placeholders).
parse_output_file_name <- function(path) {
  tokens <- strsplit(basename(path), ".", fixed = TRUE)[[1]]
  if (length(tokens) < 3 || !(tokens[1] %in% names(ARRAYS)) || !(tokens[2] %in% names(FILE_STEPS))) {
    return(NULL)
  }
  step_def <- FILE_STEPS[[tokens[2]]]
  fields <- step_def$fields
  value_tokens <- tokens[3:(length(tokens) - 1)]
  if (length(tokens) == 3) value_tokens <- character(0)

  if (is.null(step_def$dotted)) {
    if (length(value_tokens) != length(fields)) return(NULL)
    values <- value_tokens
  } else {
    if (length(value_tokens) < length(fields)) return(NULL)
    dotted_i <- which(fields == step_def$dotted)
    n_after <- length(fields) - dotted_i
    before <- value_tokens[seq_len(dotted_i - 1)]
    after <- if (n_after > 0) value_tokens[(length(value_tokens) - n_after + 1):length(value_tokens)] else character(0)
    middle <- value_tokens[dotted_i:(length(value_tokens) - n_after)]
    values <- c(before, paste(middle, collapse = "."), after)
  }

  list(array = tokens[1],
       step = tokens[2],
       fields = setNames(as.list(values), fields),
       ext = tokens[length(tokens)])
}

# Lists the outputs of one step in a folder, keeping only exact matches on the given fields.
# Fields can be vectors (e.g. source = c("AR", "BO")). Returns full paths, sorted.
find_output_files <- function(folder, array, step, ..., ext = "tsv") {
  filters <- list(...)
  get_file_step(step)
  candidates <- list.files(folder, full.names = TRUE)
  keep <- vapply(candidates, function(path) {
    parsed <- parse_output_file_name(path)
    if (is.null(parsed) || parsed$array != array || parsed$step != step || parsed$ext != ext) {
      return(FALSE)
    }
    all(vapply(names(filters), function(field) {
      !is.null(parsed$fields[[field]]) && parsed$fields[[field]] %in% filters[[field]]
    }, logical(1)))
  }, logical(1))
  sort(candidates[keep])
}

#### ARRAY COLUMN ####
# Every tabular output starts with an "array" column, so the array survives renames and rbinds
add_array_column <- function(dt, array) {
  dt <- data.table::as.data.table(dt)
  dt[, array := array]
  data.table::setcolorder(dt, c("array", setdiff(colnames(dt), "array")))
  dt
}

check_array <- function(dt, array, file = "") {
  if (!("array" %in% colnames(dt))) {
    stop(sprintf("No 'array' column in %s: was it produced by an older version of the pipeline?", file))
  }
  found <- unique(dt$array)
  if (!identical(found, array)) {
    stop(sprintf("%s contains data for array(s) %s, expected %s", file, paste(found, collapse = ", "), array))
  }
  invisible(TRUE)
}
