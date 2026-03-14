#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
})

# =========================================================
# FULL SCRIPT
# isolates scatter (GH vs ABC)
# - compute TRUE global rank across ALL databases
# - highlight global top10% (rank_ratio_global >= 0.9) in green
# - BF-enriched species (ISS rule; NO significance/FDR) as red outline squares
#   Rule implemented:
#     log2FC = log2( (BF+eps)/(B+eps) )
#     BF enriched: log2FC > 0
#     |log2FC| >= LOG2FC_CUTOFF
#     mean_abs >= MEAN_ABS_CUTOFF
# =========================================================

# -------------------------
# helpers
# -------------------------
stop_if_missing <- function(path, label = NULL) {
  if (!file.exists(path)) {
    msg <- if (!is.null(label)) paste0(label, " not found: ") else "File not found: "
    stop(msg, path, call. = FALSE)
  }
}

make_dir <- function(dir) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}

# -------------------------
# input / output paths
# -------------------------
ALL_DF_CSV   <- "data/processed/all_db_merged.csv"
GROUPS_CSV   <- "Input_Data/isolates_group.csv"
READS_CSV    <- "Input_Data/isolates_reads.csv"

# merged_table must contain species annotation: taxon_oid + Species
# set to NULL if merged_table is already loaded in memory
MERGED_TABLE_CSV <- "data/processed/isolates_input.csv"
# MERGED_TABLE_CSV <- NULL

OUT_DIR <- "results/figures"
OUT_PDF <- file.path(OUT_DIR, "Fig_S2.pdf")
OUT_PNG <- file.path(OUT_DIR, "Fig_S2.png")

# -------------------------
# params
# -------------------------
RANK_CUTOFF      <- 0.9
LOG2FC_CUTOFF    <- 1
MEAN_ABS_CUTOFF  <- 0.1
EPS              <- 1e-12

# colors / sizes
COL_BG      <- "grey85"
COL_TOP     <- "#1b7837"
COL_BOX     <- "#d73027"
PT_BG_SIZE  <- 2.0
PT_BG_ALPHA <- 0.55
PT_TOP_SIZE <- 2.6
PT_TOP_ALPHA <- 0.95
BOX_SIZE    <- 4.0
BOX_STROKE  <- 1.2

# =========================================================
# 0) checks
# =========================================================
stop_if_missing(ALL_DF_CSV, "all_db_merged.csv")
stop_if_missing(GROUPS_CSV, "isolates_group.csv")
stop_if_missing(READS_CSV,  "isolates_reads.csv")
if (!is.null(MERGED_TABLE_CSV)) stop_if_missing(MERGED_TABLE_CSV, "merged_table.csv")
make_dir(OUT_DIR)

# =========================================================
# 1) read all_df and compute TRUE GLOBAL rank across ALL databases
# =========================================================
all_df_raw <- readr::read_csv(ALL_DF_CSV, show_col_types = FALSE)

need_all <- c("database", "genome_id", "GH_total_abundance", "ABC_total_abundance", "ABC_GH_ratio")
miss_all <- setdiff(need_all, names(all_df_raw))
if (length(miss_all) > 0) {
  stop("Missing columns in all_df: ", paste(miss_all, collapse = ", "), call. = FALSE)
}

all_df <- all_df_raw %>%
  mutate(
    database = as.character(database),
    genome_id = as.character(genome_id),
    GH_total_abundance  = suppressWarnings(as.numeric(GH_total_abundance)),
    ABC_total_abundance = suppressWarnings(as.numeric(ABC_total_abundance)),
    ABC_GH_ratio        = suppressWarnings(as.numeric(ABC_GH_ratio))
  ) %>%
  filter(
    is.finite(GH_total_abundance),
    is.finite(ABC_total_abundance),
    is.finite(ABC_GH_ratio)
  ) %>%
  mutate(
    rank_ratio_global = dplyr::percent_rank(ABC_GH_ratio)
  )

cat("Global pool rows:", nrow(all_df), "\n")
cat("Global top10 rows:", sum(all_df$rank_ratio_global >= RANK_CUTOFF), "\n")
cat("Rows per database:\n")
print(all_df %>% count(database))

# =========================================================
# 2) isolates subset (KEEP global rank)
# =========================================================
df_iso <- all_df %>%
  filter(database == "isolates") %>%
  mutate(is_top10 = rank_ratio_global >= RANK_CUTOFF)

if (nrow(df_iso) == 0) stop("No isolates rows after filtering. Check database labels.", call. = FALSE)

cat("Isolates rows:", nrow(df_iso), "\n")
cat("Isolates top10 (global):", sum(df_iso$is_top10), "\n")

# =========================================================
# 3) build abs_wide from reads + groups (NO significance)
# =========================================================
groups <- read.csv(GROUPS_CSV, header = TRUE, stringsAsFactors = FALSE)
reads  <- read.csv(READS_CSV,  header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

need_groups <- c("sample", "member", "Abundance")
miss_groups <- setdiff(need_groups, names(groups))
if (length(miss_groups) > 0) {
  stop("Missing columns in groups: ", paste(miss_groups, collapse = ", "), call. = FALSE)
}
if (!("Species.names" %in% names(reads))) {
  stop("reads must contain column: Species.names", call. = FALSE)
}

# clean column names
names(reads)  <- names(reads)  %>% str_replace_all(" ", ".")
names(groups) <- names(groups) %>% str_replace_all(" ", ".")

# ensure member labels are exactly B / BF
groups <- groups %>%
  mutate(
    sample = as.character(sample),
    member = as.character(member),
    Abundance = suppressWarnings(as.numeric(Abundance))
  ) %>%
  filter(member %in% c("B", "BF"), is.finite(Abundance))

# match samples
sample_cols <- intersect(names(reads), groups$sample)
if (length(sample_cols) == 0) {
  stop("No overlapping sample names between reads columns and groups$sample.\n",
       "Check sample naming/formatting.", call. = FALSE)
}

# convert sample columns to numeric
reads <- reads %>%
  mutate(across(all_of(sample_cols), ~ suppressWarnings(as.numeric(.x))))

# long table + absolute abundance
dat_rel <- reads %>%
  dplyr::select(Species.names, all_of(sample_cols)) %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "sample",
    values_to = "rel"
  ) %>%
  mutate(
    Species.names = as.character(Species.names),
    sample = as.character(sample),
    rel = suppressWarnings(as.numeric(rel))
  ) %>%
  left_join(groups %>% dplyr::select(sample, member, Abundance), by = "sample") %>%
  filter(is.finite(rel), !is.na(member), is.finite(Abundance)) %>%
  mutate(abs = rel * Abundance)

# mean absolute abundance per species per member
mean_abs_by_member <- dat_rel %>%
  group_by(Species.names, member) %>%
  summarise(
    mean_abs = mean(abs, na.rm = TRUE),
    .groups = "drop"
  )

# wide table: B and BF columns
abs_wide <- mean_abs_by_member %>%
  pivot_wider(
    names_from = member,
    values_from = mean_abs
  )

# =========================================================
# 4) BF enriched species table (ISS rule; NO significance)
# =========================================================
enrich_BF_tbl <- abs_wide %>%
  transmute(
    Species.names = as.character(Species.names),
    B  = coalesce(suppressWarnings(as.numeric(B)),  0),
    BF = coalesce(suppressWarnings(as.numeric(BF)), 0),
    mean_abundance = pmax(B, BF),
    log2FC = log2((BF + EPS) / (B + EPS))
  ) %>%
  filter(
    log2FC > 0,
    abs(log2FC) >= LOG2FC_CUTOFF,
    mean_abundance >= MEAN_ABS_CUTOFF
  ) %>%
  distinct(Species.names, .keep_all = TRUE)

cat("BF enriched species (rule-based) n =", nrow(enrich_BF_tbl), "\n")

# =========================================================
# 5) map isolates genome_id -> Species (via merged_table)
# =========================================================
if (!is.null(MERGED_TABLE_CSV)) {
  merged_table <- readr::read_csv(MERGED_TABLE_CSV, show_col_types = FALSE)
}

if (!exists("merged_table")) {
  stop("merged_table not found.\n",
       "Provide MERGED_TABLE_CSV or load merged_table before running.", call. = FALSE)
}

need_mt <- c("taxon_oid", "Species")
miss_mt <- setdiff(need_mt, names(merged_table))
if (length(miss_mt) > 0) {
  stop("merged_table missing columns: ", paste(miss_mt, collapse = ", "), call. = FALSE)
}

species_map <- merged_table %>%
  transmute(
    genome_id = as.character(taxon_oid),
    Species = as.character(Species),
    Species_clean = sub("^s__", "", Species)
  ) %>%
  distinct(genome_id, .keep_all = TRUE)

df_iso2 <- df_iso %>%
  left_join(species_map, by = "genome_id") %>%
  mutate(
    Species_clean = as.character(Species_clean)
  ) %>%
  left_join(enrich_BF_tbl, by = c("Species_clean" = "Species.names")) %>%
  mutate(
    is_BF_enriched = !is.na(log2FC)
  )

cat("Isolates with Species mapped:", sum(!is.na(df_iso2$Species_clean)), " / ", nrow(df_iso2), "\n")
cat("Isolates BF-enriched (mapped):", sum(df_iso2$is_BF_enriched, na.rm = TRUE), "\n")
cat("Isolates top10 (global) & BF-enriched:", sum(df_iso2$is_top10 & df_iso2$is_BF_enriched, na.rm = TRUE), "\n")

# =========================================================
# 6) plot
# =========================================================
p <- ggplot(df_iso2, aes(x = GH_total_abundance, y = ABC_total_abundance)) +
  # layer 1: all isolates (grey)
  geom_point(
    color = COL_BG,
    size = PT_BG_SIZE,
    alpha = PT_BG_ALPHA
  ) +
  # layer 2: top10% global rank (green)
  geom_point(
    data = df_iso2 %>% filter(is_top10),
    color = COL_TOP,
    size = PT_TOP_SIZE,
    alpha = PT_TOP_ALPHA
  ) +
  # layer 3: BF enriched species (red outline square)
  geom_point(
    data = df_iso2 %>% filter(is_BF_enriched),
    shape = 22,
    fill = NA,
    color = COL_BOX,
    size = BOX_SIZE,
    stroke = BOX_STROKE
  ) +
  labs(
    x = "Complex carbohydrate depolymerization",
    y = "Carbohydrate transporters",
    title = paste0(
      "Isolates: global rank >= ", RANK_CUTOFF, " (green)\n",
      "BF enriched species (|log2FC| >= ", LOG2FC_CUTOFF,
      "; mean abs >= ", MEAN_ABS_CUTOFF, ")"
    )
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(size = 15, face = "bold"),
    axis.text  = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )

print(p)

ggsave(OUT_PDF, p, width = 7.5, height = 6)
ggsave(OUT_PNG, p, width = 7.5, height = 6, dpi = 300)

message("Saved:\n - ", OUT_PDF, "\n - ", OUT_PNG, "\nDone ✅")