#!/usr/bin/env Rscript
# Checks for the output file-name builder/parser in functions/00_common.R
# Run from the code/ folder: Rscript tests/test_file_names.R

source("functions/00_common.R")

expect_error <- function(expr) {
  failed <- tryCatch({ expr; FALSE }, error = function(e) TRUE)
  if (!failed) stop("expected an error: ", deparse(substitute(expr)))
}

#### Labels ####
stopifnot(smoothing_label(3, 0, "zeros") == "med3_mean0_zeros")
stopifnot(cutoff_label(4, 2) == "cutoff4SD_2pep")
stopifnot(cutoff_label(4.5, 2) == "cutoff4p5SD_2pep")
stopifnot(fmt_num(10) == "10", fmt_num(100000) == "100000")

#### Round trip for every step ####
examples <- list(
  list("v1", "normalized", source = "AR", type = "PO"),
  list("v2", "normalized", source = "AR_P1", type = "PO"),
  list("v1", "global_statistics"),
  list("v1", "smoothed", source = "LE", type = "NE", smoothing = "med3_mean0_zeros"),
  list("v1", "peaks", smoothing = "med3_mean0_zeros", cutoff = "cutoff4p5SD_2pep"),
  list("v1", "regions_peptides", smoothing = "med3_mean0_zeros", cutoff = "cutoff4SD_2pep"),
  list("v1", "regions_summary", smoothing = "med3_mean0_zeros", cutoff = "cutoff4SD_2pep"),
  list("v2", "alascan_change", protein = "TcCLB.511671.60"),
  list("v2", "alascan_matrix", protein = "TcCLB.511671.60", group = "1"),
  list("v2", "alascan_heatmap", protein = "TCSYLVIO_004530", group = "12"),
  list("v2", "alascan_logo", protein = "TcCLB.511671.60", group = "3", sera = "AR_P1+BO_P1"),
  list("v1", "profiles", smoothing = "med3_mean0_zeros", sd = "4SD", scope = "all"),
  list("v1", "profiles_index", smoothing = "med3_mean0_zeros", sd = "4SD", scope = "TcCLB.511671.60"),
  list("v2", "ridgeplots", smoothing = "med3_mean0_zeros", scope = "TcCLB.511671.60"),
  list("v2", "ridgeplots_index", smoothing = "med3_mean0_zeros", scope = "all"),
  list("v2", "signal_distribution", smoothing = "med3_mean0_zeros", sign = "positives")
)
for (example in examples) {
  array <- example[[1]]
  step <- example[[2]]
  fields <- example[-(1:2)]
  for (ext in c("tsv", "pdf")) {
    path <- do.call(output_file_name, c(list(folder = "some/folder", array = array, step = step), fields, list(ext = ext)))
    stopifnot(dirname(path) == "some/folder")
    parsed <- parse_output_file_name(path)
    stopifnot(!is.null(parsed), parsed$array == array, parsed$step == step, parsed$ext == ext)
    stopifnot(identical(as.character(names(parsed$fields)), as.character(names(fields))),
              identical(unlist(parsed$fields), unlist(lapply(fields, as.character))))
    rebuilt <- do.call(output_file_name, c(list(folder = "some/folder", array = parsed$array, step = parsed$step), parsed$fields, list(ext = parsed$ext)))
    stopifnot(rebuilt == path)
  }
}

stopifnot(basename(output_file_name("x", "v2", "alascan_matrix", protein = "TcCLB.511671.60", group = 1)) ==
            "v2.alascan_matrix.TcCLB.511671.60.1.tsv")

#### Invalid names ####
expect_error(output_file_name("x", "v1", "smoothed", source = "AR", type = "PO", smoothing = "med3.5"))
expect_error(output_file_name("x", "v1", "smoothed", source = "AR", type = "PO"))
expect_error(output_file_name("x", "v1", "smoothed", source = "AR", type = "PO", smoothing = "s", extra = "y"))
expect_error(output_file_name("x", "v1", "normalized", source = "", type = "PO"))
expect_error(output_file_name("x", "v1", "normalized", source = "a/b", type = "PO"))
expect_error(output_file_name("x", "v1", "no_such_step"))
expect_error(get_array("v3"))
expect_error(get_array("v2", supported = "v1"))
stopifnot(get_array("v1")$peptide_offset == 4, get_array("v2")$peptide_offset == 1)

stopifnot(is.null(parse_output_file_name("OUTPUT WILL BE GENERATED HERE")))
stopifnot(is.null(parse_output_file_name("AR_PO_processed.tsv")))
stopifnot(is.null(parse_output_file_name("v1.normalized.AR.tsv")))

#### find_output_files: exact matching ####
folder <- tempfile("names_test")
dir.create(folder)
files <- c(output_file_name(folder, "v1", "smoothed", source = "AR", type = "PO", smoothing = "med3_mean0_zeros"),
           output_file_name(folder, "v1", "smoothed", source = "AR", type = "NE", smoothing = "med3_mean0_zeros"),
           output_file_name(folder, "v1", "smoothed", source = "AR", type = "PO", smoothing = "med3_mean3_zeros"),
           output_file_name(folder, "v2", "smoothed", source = "AR_P1", type = "PO", smoothing = "med3_mean0_zeros"),
           file.path(folder, "OUTPUT WILL BE GENERATED HERE"))
invisible(file.create(files))

found <- find_output_files(folder, "v1", "smoothed", source = "AR", smoothing = "med3_mean0_zeros")
stopifnot(length(found) == 2)
stopifnot(length(find_output_files(folder, "v1", "smoothed", source = "AR", type = "PO", smoothing = "med3_mean3_zeros")) == 1)
stopifnot(length(find_output_files(folder, "v2", "smoothed", source = "AR")) == 0)
stopifnot(length(find_output_files(folder, "v2", "smoothed", source = c("AR_P1", "BO_P1"))) == 1)
stopifnot(length(find_output_files(folder, "v1", "smoothed", ext = "pdf")) == 0)
unlink(folder, recursive = TRUE)

#### array column ####
dt <- add_array_column(data.table::data.table(a = 1:2), "v1")
stopifnot(colnames(dt)[1] == "array", check_array(dt, "v1"))
expect_error(check_array(dt, "v2"))
expect_error(check_array(data.table::data.table(a = 1), "v1"))

cat("All file-name tests passed\n")
