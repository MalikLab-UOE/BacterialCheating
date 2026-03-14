#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
  library(openxlsx)
  library(dplyr)
  library(purrr)
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

# robust numeric conversion for selected columns
to_numeric_cols <- function(df, cols) {
  if (length(cols) == 0) return(df)
  df %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(cols), ~ suppressWarnings(as.numeric(.x))))
}

# safe rowSums for all columns except id
add_total <- function(df, id_col, total_name) {
  gene_cols <- setdiff(names(df), id_col)
  if (length(gene_cols) == 0) {
    df[[total_name]] <- 0
    return(df)
  }
  df <- to_numeric_cols(df, gene_cols)
  df %>% dplyr::mutate("{total_name}" := rowSums(dplyr::across(dplyr::all_of(gene_cols)), na.rm = TRUE))
}

# fill taxonomy "NA"/"" -> NA and then fill down to "unclassified ..."
fill_taxonomy <- function(df) {
  df[df == "NA" | df == ""] <- NA
  
  # only do if these columns exist
  tax_cols <- c("Phylum", "Class", "Order", "Family", "Genus")
  have <- tax_cols[tax_cols %in% names(df)]
  if (length(have) == 0) return(df)
  
  # fill in hierarchical way (your original logic)
  df %>%
    dplyr::mutate(
      Class  = if ("Class"  %in% names(.)) ifelse(is.na(Class),
                                                  ifelse(!is.na(Phylum), paste("unclassified", Phylum), "unclassified"),
                                                  Class) else Class,
      Order  = if ("Order"  %in% names(.)) ifelse(is.na(Order),
                                                  ifelse(!is.na(Class), paste("unclassified", Class), "unclassified"),
                                                  Order) else Order,
      Family = if ("Family" %in% names(.)) ifelse(is.na(Family),
                                                  ifelse(!is.na(Order), paste("unclassified", Order), "unclassified"),
                                                  Family) else Family,
      Genus  = if ("Genus"  %in% names(.)) ifelse(is.na(Genus),
                                                  ifelse(!is.na(Family), paste("unclassified", Family), "unclassified"),
                                                  Genus) else Genus
    )
}

# ---------------------------
# args / default paths
# ---------------------------
args <- commandArgs(trailingOnly = TRUE)

hmm_csv   <- ifelse(length(args) >= 1, args[1], "Input_Data/hmm_Loma.csv")
gh_xlsx   <- ifelse(length(args) >= 2, args[2], "Input_Data/microtrait_GH.xlsx")
abc_xlsx  <- ifelse(length(args) >= 3, args[3], "Input_Data/microtrait_transporters.xlsx")
abc_sheet <- ifelse(length(args) >= 4, as.integer(args[4]), 2)

gtdb_csv  <- ifelse(length(args) >= 5, args[5], "Input_Data/loma_gtdbtk.bac120.summary_r226.csv")

bins_xlsx <- ifelse(length(args) >= 6, args[6],
                    "Input_Data/loma_output_seqtk_MegaHITcoa_MW121_BinsStats_45_10.xlsx")

min_comp  <- ifelse(length(args) >= 7, as.numeric(args[7]), 50)
max_cont  <- ifelse(length(args) >= 8, as.numeric(args[8]), 10)

out_rds   <- ifelse(length(args) >= 9,  args[9],  "data/processed/loma_input.rds")
out_csv   <- ifelse(length(args) >= 10, args[10], "data/processed/loma_input.csv")

message("=== Build LOMA MAG input (GH/ABC totals + taxonomy + binQC) ===")
message("HMM CSV:     ", hmm_csv)
message("GH rules:    ", gh_xlsx)
message("ABC rules:   ", abc_xlsx, " (sheet=", abc_sheet, ")")
message("GTDB csv:    ", gtdb_csv)
message("Bins xlsx:   ", bins_xlsx)
message("QC:          completeness >= ", min_comp, ", contamination <= ", max_cont)
message("Out RDS:     ", out_rds)
message("Out CSV:     ", out_csv)

stop_if_missing(hmm_csv,  "HMM CSV")
stop_if_missing(gh_xlsx,  "GH rules xlsx")
stop_if_missing(abc_xlsx, "ABC rules xlsx")
stop_if_missing(gtdb_csv, "GTDB summary csv")
stop_if_missing(bins_xlsx,"BinsStats xlsx")

# ---------------------------
# read inputs
# ---------------------------
hmm_data  <- readr::read_csv(hmm_csv, show_col_types = FALSE)
gh_traits <- readxl::read_excel(gh_xlsx)
abc_traits <- readxl::read_excel(abc_xlsx, sheet = abc_sheet)
gtdb <- readr::read_csv(gtdb_csv, show_col_types = FALSE)

# validate required columns
if (!("id" %in% names(hmm_data))) stop("Column 'id' missing in hmm_Loma.csv", call. = FALSE)
if (!("GH" %in% names(gh_traits))) stop("Column 'GH' missing in microtrait_GH.xlsx", call. = FALSE)
if (!("ABC" %in% names(abc_traits))) stop("Column 'ABC' missing in microtrait_transporters.xlsx (selected sheet)", call. = FALSE)
if (!("Bin.ID" %in% names(gtdb))) stop("Column 'Bin.ID' missing in gtdbtk summary", call. = FALSE)

# gene lists
gh_genes  <- unique(stats::na.omit(gh_traits$GH))
abc_genes <- unique(stats::na.omit(abc_traits$ABC))

gh_cols  <- intersect(names(hmm_data), gh_genes)
abc_cols <- intersect(names(hmm_data), abc_genes)

message("Matched GH cols:  ", length(gh_cols))
message("Matched ABC cols: ", length(abc_cols))
if (length(gh_cols) == 0) stop("No GH columns matched in HMM table.", call. = FALSE)
if (length(abc_cols) == 0) stop("No ABC columns matched in HMM table.", call. = FALSE)

# ---------------------------
# compute totals
# ---------------------------
gh_only <- hmm_data %>%
  dplyr::select(id, dplyr::all_of(gh_cols)) %>%
  add_total(id_col = "id", total_name = "GH_total_abundance") %>%
  dplyr::select(id, GH_total_abundance)

abc_only <- hmm_data %>%
  dplyr::select(id, dplyr::all_of(abc_cols)) %>%
  add_total(id_col = "id", total_name = "ABC_total_abundance") %>%
  dplyr::select(id, ABC_total_abundance)

merged_abundance <- dplyr::full_join(gh_only, abc_only, by = "id")

# ---------------------------
# join GTDB taxonomy
# ---------------------------
merged_table <- dplyr::left_join(gtdb, merged_abundance, by = c("Bin.ID" = "id"))

# ---------------------------
# read bins xlsx (all sheets) + QC filter
# ---------------------------
wb <- openxlsx::loadWorkbook(bins_xlsx)
sheets <- names(wb)

read_and_tag <- function(sh) {
  df <- openxlsx::read.xlsx(wb, sheet = sh)
  if (!("bin" %in% names(df))) stop("Sheet ", sh, " is missing column 'bin'.", call. = FALSE)
  
  df %>%
    dplyr::mutate(
      bin = paste0(sh, "_", bin),
      sample = sh
    )
}

merged_bins <- purrr::map_dfr(sheets, read_and_tag)

# validate QC columns
need_qc <- c("completeness", "contamination")
miss_qc <- setdiff(need_qc, names(merged_bins))
if (length(miss_qc) > 0) stop("BinsStats missing columns: ", paste(miss_qc, collapse = ", "), call. = FALSE)

filtered_bins <- merged_bins %>%
  dplyr::filter(completeness >= min_comp, contamination <= max_cont)

# ---------------------------
# join QC-filtered bins to main table
# ---------------------------
merged_table <- merged_table %>%
  dplyr::inner_join(filtered_bins, by = c("Bin.ID" = "bin"))

# ---------------------------
# filter non-zero + compute ratio
# ---------------------------
merged_table <- merged_table %>%
  dplyr::mutate(
    GH_total_abundance  = dplyr::coalesce(GH_total_abundance, 0),
    ABC_total_abundance = dplyr::coalesce(ABC_total_abundance, 0)
  ) %>%
  dplyr::filter(GH_total_abundance != 0, ABC_total_abundance != 0) %>%
  dplyr::mutate(
    ABC_GH_ratio = ABC_total_abundance / GH_total_abundance,
    ABC_GH_ratio = dplyr::if_else(is.finite(ABC_GH_ratio), ABC_GH_ratio, NA_real_)
  )

# taxonomy cleanup (optional)
merged_table <- fill_taxonomy(merged_table)

message("Final rows: ", nrow(merged_table))

# ---------------------------
# save
# ---------------------------
make_dir(dirname(out_rds))
make_dir(dirname(out_csv))

saveRDS(merged_table, out_rds)
readr::write_csv(merged_table, out_csv)

message("Done ✅")