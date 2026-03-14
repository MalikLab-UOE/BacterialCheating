#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(readr)
  library(tools)
})

# ---------------------------
# helpers
# ---------------------------
stop_if_missing <- function(path, label = NULL) {
  if (!file.exists(path)) {
    msg <- if (!is.null(label)) paste0(label, " not found: ") else "File not found: "
    stop(msg, path, call. = FALSE)
  }
}

make_dir <- function(dir) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}

# convert to numeric and report NA inflation
as_numeric_report <- function(df, cols, label = "") {
  if (length(cols) == 0) return(df)
  before_na <- sum(is.na(as.matrix(df[, cols, drop = FALSE])))
  out <- df %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(cols), ~ suppressWarnings(as.numeric(.x))))
  after_na <- sum(is.na(as.matrix(out[, cols, drop = FALSE])))
  delta <- after_na - before_na
  if (delta > 0) {
    warning(label, ": numeric conversion introduced ", delta, " additional NA values.")
  }
  out
}

# ---------------------------
# args / paths
# ---------------------------
args <- commandArgs(trailingOnly = TRUE)

hmm_csv   <- ifelse(length(args) >= 1, args[1], "Input_Data/hmm_isolates.csv")
gh_xlsx   <- ifelse(length(args) >= 2, args[2], "Input_Data/microtrait_GH.xlsx")
abc_xlsx  <- ifelse(length(args) >= 3, args[3], "Input_Data/microtrait_transporters.xlsx")
abc_sheet <- ifelse(length(args) >= 4, as.integer(args[4]), 2)
meta_csv  <- ifelse(length(args) >= 5, args[5], "Input_Data/isolates_metadata_r226.csv")

out_rds   <- ifelse(length(args) >= 6, args[6], "data/processed/isolates_input.rds")
out_csv   <- ifelse(length(args) >= 7, args[7], "data/processed/isolates_input.csv")

message("=== Build isolates input (GH/ABC totals + ratio) ===")
message("HMM CSV:        ", hmm_csv)
message("GH rules xlsx:  ", gh_xlsx)
message("ABC rules xlsx: ", abc_xlsx, " (sheet=", abc_sheet, ")")
message("Meta CSV:       ", meta_csv)
message("Output RDS:     ", out_rds)
message("Output CSV:     ", out_csv)

stop_if_missing(hmm_csv,  "HMM CSV")
stop_if_missing(gh_xlsx,  "GH xlsx")
stop_if_missing(abc_xlsx, "ABC xlsx")
stop_if_missing(meta_csv, "Metadata CSV")

# ---------------------------
# load
# ---------------------------
hmm_data   <- readr::read_csv(hmm_csv, show_col_types = FALSE)
gh_traits  <- readxl::read_excel(gh_xlsx)
abc_traits <- readxl::read_excel(abc_xlsx, sheet = abc_sheet)
img_meta   <- readr::read_csv(meta_csv, show_col_types = FALSE)

# validate columns
if (!("id" %in% names(hmm_data))) {
  stop("Column 'id' missing in hmm csv.", call. = FALSE)
}
if (!("taxon_oid" %in% names(img_meta))) {
  stop("Column 'taxon_oid' missing in metadata.", call. = FALSE)
}
if (!("GH" %in% names(gh_traits))) {
  stop("Column 'GH' missing in GH rules xlsx.", call. = FALSE)
}
if (!("ABC" %in% names(abc_traits))) {
  stop("Column 'ABC' missing in ABC rules xlsx (selected sheet).", call. = FALSE)
}

# ---------------------------
# gene lists + matched columns
# ---------------------------
gh_genes  <- unique(stats::na.omit(gh_traits$GH))
abc_genes <- unique(stats::na.omit(abc_traits$ABC))

gh_columns  <- intersect(names(hmm_data), gh_genes)
abc_columns <- intersect(names(hmm_data), abc_genes)

message("Matched GH columns:  ", length(gh_columns))
message("Matched ABC columns: ", length(abc_columns))

if (length(gh_columns) == 0) stop("No GH columns matched between hmm_data and GH list.", call. = FALSE)
if (length(abc_columns) == 0) stop("No ABC columns matched between hmm_data and ABC list.", call. = FALSE)

# ---------------------------
# compute totals
# ---------------------------
gh_only <- hmm_data %>%
  dplyr::select(id, dplyr::all_of(gh_columns))
gh_only <- as_numeric_report(gh_only, gh_columns, label = "GH")
gh_only <- gh_only %>%
  dplyr::mutate(GH_total_abundance = rowSums(dplyr::across(dplyr::all_of(gh_columns)), na.rm = TRUE)) %>%
  dplyr::select(id, GH_total_abundance)

abc_only <- hmm_data %>%
  dplyr::select(id, dplyr::all_of(abc_columns))
abc_only <- as_numeric_report(abc_only, abc_columns, label = "ABC")
abc_only <- abc_only %>%
  dplyr::mutate(ABC_total_abundance = rowSums(dplyr::across(dplyr::all_of(abc_columns)), na.rm = TRUE)) %>%
  dplyr::select(id, ABC_total_abundance)

merged_abundance <- dplyr::full_join(gh_only, abc_only, by = "id")

# ---------------------------
# join metadata + compute ratio
# ---------------------------
merged_table <- dplyr::left_join(img_meta, merged_abundance, by = c("taxon_oid" = "id")) %>%
  dplyr::mutate(
    GH_total_abundance  = dplyr::coalesce(GH_total_abundance, 0),
    ABC_total_abundance = dplyr::coalesce(ABC_total_abundance, 0)
  ) %>%
  dplyr::filter(GH_total_abundance > 0, ABC_total_abundance > 0) %>%
  dplyr::mutate(
    ABC_GH_ratio = ABC_total_abundance / GH_total_abundance,
    ABC_GH_ratio = dplyr::if_else(is.finite(ABC_GH_ratio), ABC_GH_ratio, NA_real_)
  )

message("Final rows: ", nrow(merged_table))

# ---------------------------
# save
# ---------------------------
make_dir(dirname(out_rds))
make_dir(dirname(out_csv))

saveRDS(merged_table, out_rds)
readr::write_csv(merged_table, out_csv)

message("Done ✅")