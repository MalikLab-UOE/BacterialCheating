#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tools)
})

# ---------------------------
# Helpers
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

pick_columns_by_gene_list <- function(hmm_df, gene_list, id_col = "id") {
  gene_list <- unique(na.omit(gene_list))
  cols <- intersect(colnames(hmm_df), gene_list)
  
  out <- hmm_df %>%
    select(all_of(id_col), all_of(cols))
  
  list(df = out, matched_cols = cols)
}

add_total_abundance <- function(df, id_col = "id", total_col = "total_abundance") {
  # If there are no gene columns matched, rowSums() would fail; handle gracefully
  gene_cols <- setdiff(colnames(df), id_col)
  if (length(gene_cols) == 0) {
    df[[total_col]] <- 0
    return(df)
  }
  
  df %>%
    mutate(
      "{total_col}" := rowSums(select(., all_of(gene_cols)), na.rm = TRUE)
    )
}

# ---------------------------
# Main
# ---------------------------
args <- commandArgs(trailingOnly = TRUE)

# Default paths (match your current script)
hmm_csv      <- ifelse(length(args) >= 1, args[1], "Input_Data/hmm_MAG_IMG.csv")
gh_xlsx      <- ifelse(length(args) >= 2, args[2], "Input_Data/microtrait_GH.xlsx")
abc_xlsx     <- ifelse(length(args) >= 3, args[3], "Input_Data/microtrait_transporters.xlsx")
abc_sheet    <- ifelse(length(args) >= 4, as.integer(args[4]), 2)
meta_csv     <- ifelse(length(args) >= 5, args[5], "Input_Data/IMG_bindata_withmeta_norestricted1.csv")
out_rds      <- ifelse(length(args) >= 6, args[6], "data/processed/img_input.rds")
out_csv      <- ifelse(length(args) >= 7, args[7], "data/processed/img_input.csv")

message("=== Build img_input (RDS/CSV) ===")
message("HMM CSV:      ", hmm_csv)
message("GH rules:     ", gh_xlsx)
message("ABC rules:    ", abc_xlsx, " (sheet=", abc_sheet, ")")
message("Meta CSV:     ", meta_csv)
message("Output RDS:   ", out_rds)
message("Output CSV:   ", out_csv)

# Check inputs
stop_if_missing(hmm_csv,  "HMM CSV")
stop_if_missing(gh_xlsx,  "GH xlsx")
stop_if_missing(abc_xlsx, "ABC xlsx")
stop_if_missing(meta_csv, "Meta CSV")

# Read inputs
hmm_data  <- read.csv(hmm_csv, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
gh_traits <- read_excel(gh_xlsx)
abc_traits <- read_excel(abc_xlsx, sheet = abc_sheet)
img_meta  <- read.csv(meta_csv, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

# Validate expected columns
if (!("id" %in% colnames(hmm_data))) stop("Column 'id' is missing in HMM CSV.", call. = FALSE)
if (!("Bin.ID" %in% colnames(img_meta))) stop("Column 'Bin.ID' is missing in Meta CSV.", call. = FALSE)
if (!("Domain" %in% colnames(img_meta))) stop("Column 'Domain' is missing in Meta CSV.", call. = FALSE)

# Extract gene lists (be forgiving about column name)
if (!("GH" %in% colnames(gh_traits))) stop("Column 'GH' is missing in GH rules xlsx.", call. = FALSE)
if (!("ABC" %in% colnames(abc_traits))) stop("Column 'ABC' is missing in ABC rules xlsx (selected sheet).", call. = FALSE)

gh_genes  <- gh_traits$GH
abc_genes <- abc_traits$ABC

# Pick columns & compute totals
gh_pick <- pick_columns_by_gene_list(hmm_data, gh_genes, id_col = "id")
gh_only <- add_total_abundance(gh_pick$df, id_col = "id", total_col = "GH_total_abundance")

abc_pick <- pick_columns_by_gene_list(hmm_data, abc_genes, id_col = "id")
abc_only <- add_total_abundance(abc_pick$df, id_col = "id", total_col = "ABC_total_abundance")

message("Matched GH columns:  ", length(gh_pick$matched_cols))
message("Matched ABC columns: ", length(abc_pick$matched_cols))

# Merge totals (and optionally keep all gene columns you selected)
merged_abundance <- full_join(gh_only, abc_only, by = "id")

# Join meta
merged_table <- left_join(img_meta, merged_abundance, by = c("Bin.ID" = "id"))

# Filter non-zero & bacteria
merged_table <- merged_table %>%
  filter(
    !is.na(ABC_total_abundance),
    !is.na(GH_total_abundance),
    ABC_total_abundance != 0,
    GH_total_abundance != 0,
    Domain == "Bacteria"
  ) %>%
  mutate(
    ABC_GH_ratio = ABC_total_abundance / GH_total_abundance
  )

message("Final rows after filters: ", nrow(merged_table))

# Write outputs
make_dir(dirname(out_rds))
make_dir(dirname(out_csv))

saveRDS(merged_table, file = out_rds)
write.csv(merged_table, file = out_csv, row.names = FALSE)

message("Done ✅")