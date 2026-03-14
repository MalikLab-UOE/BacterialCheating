#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(scales)
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

# unify columns (keeps Phylum + genome_size if present)
unify_cols <- function(df){
  df %>%
    mutate(
      taxon_oid = if ("taxon_oid" %in% names(.)) as.character(.data$taxon_oid) else NA_character_,
      id        = if ("id"        %in% names(.)) as.character(.data$id)        else NA_character_,
      Bin.ID    = if ("Bin.ID"    %in% names(.)) as.character(.data$`Bin.ID`)  else NA_character_,
      genome_id = coalesce(taxon_oid, id, Bin.ID),
      
      Phylum = if ("Phylum" %in% names(.)) as.character(.data$Phylum) else NA_character_,
      genome_size = if ("genome_size" %in% names(.)) suppressWarnings(as.numeric(.data$genome_size)) else NA_real_,
      
      ABC_total_abundance = suppressWarnings(as.numeric(.data$ABC_total_abundance)),
      GH_total_abundance  = suppressWarnings(as.numeric(.data$GH_total_abundance)),
      ABC_GH_ratio        = suppressWarnings(as.numeric(.data$ABC_GH_ratio))
    ) %>%
    dplyr::select(
      genome_id, database, Phylum, genome_size,
      GH_total_abundance, ABC_total_abundance, ABC_GH_ratio
    )
}

# ---------------------------
# inputs
# ---------------------------
img_rds  <- "data/processed/img_input.rds"
loma_rds <- "data/processed/loma_input.rds"
iso_rds  <- "data/processed/isolates_input.rds"
glo_rds  <- "data/processed/glo_input.rds"

out_dir <- "results/figures"
out_pdf <- file.path(out_dir, "Fig_S5.pdf")
out_png <- file.path(out_dir, "Fig_S5.png")

stop_if_missing(img_rds,  "img_input.rds")
stop_if_missing(loma_rds, "loma_input.rds")
stop_if_missing(iso_rds,  "isolates_input.rds")
stop_if_missing(glo_rds,  "glo_input.rds")
make_dir(out_dir)

# ---------------------------
# load + tag database
# ---------------------------
img  <- readRDS(img_rds)  %>% mutate(database = "IMG")
loma <- readRDS(loma_rds) %>% mutate(database = "LOMA")
iso  <- readRDS(iso_rds)  %>% mutate(database = "isolates")
glo  <- readRDS(glo_rds)  %>% mutate(database = "glo")  # keep lowercase consistent

img  <- unify_cols(img)
loma <- unify_cols(loma)
iso  <- unify_cols(iso)
glo  <- unify_cols(glo)

# ---------------------------
# merge + clean (global)
# ---------------------------
db_levels <- c("LOMA", "IMG", "isolates", "glo")

all_df <- bind_rows(loma, img, iso, glo) %>%
  mutate(database = factor(database, levels = db_levels)) %>%
  filter(
    is.finite(ABC_GH_ratio),
    is.finite(ABC_total_abundance),
    is.finite(GH_total_abundance)
  )

if (nrow(all_df) == 0) stop("No rows after filtering for ratio/abundance.", call. = FALSE)

# ---------------------------
# global ranked cheating index + highlight flag
# ---------------------------
all_df <- all_df %>%
  mutate(
    rank_ratio_global = percent_rank(ABC_GH_ratio),
    highlight = ifelse(rank_ratio_global >= 0.9, "Top10", "Other")
  )

# ---------------------------
# subset to glo + clean Phylum name (remove p__)
# ---------------------------
glo_df <- all_df %>%
  filter(database == "glo") %>%
  filter(!is.na(Phylum), Phylum != "") %>%
  mutate(Phylum = sub("^p__", "", Phylum))

if (nrow(glo_df) == 0) stop("No glo rows with Phylum available.", call. = FALSE)

# ---------------------------
# choose Phylum order by ranked cheating index median (within glo)
# ---------------------------
TOP_N_PHYLUM <- 30
MIN_N_PHYLUM <- 1000

phylum_stat <- glo_df %>%
  group_by(Phylum) %>%
  summarise(
    n = n(),
    median_ranked = median(rank_ratio_global, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n >= MIN_N_PHYLUM) %>%
  arrange(desc(median_ranked)) %>%
  slice_head(n = TOP_N_PHYLUM)

if (nrow(phylum_stat) == 0) {
  stop("No phylum passed MIN_N_PHYLUM = ", MIN_N_PHYLUM, " in glo.", call. = FALSE)
}

# Build y-axis labels: "Phylum (n=xxx)"
phylum_stat <- phylum_stat %>%
  mutate(label = paste0(Phylum, " (n=", n, ")"))

# keep only selected phyla and attach label
glo_df <- glo_df %>%
  semi_join(phylum_stat %>% dplyr::select(Phylum), by = "Phylum") %>%
  left_join(phylum_stat %>% dplyr::select(Phylum, label), by = "Phylum") %>%
  mutate(
    label = factor(label, levels = rev(phylum_stat$label))
  )

# ---------------------------
# pivot long for 4 panels (ORDER: ABC, GH, ranked, genome_size LAST)
# ---------------------------
df_long <- glo_df %>%
  transmute(
    genome_id,
    label,
    highlight,
    `Carbohydrate transporters` = ABC_total_abundance,
    `Complex carbohydrate depolymerization` = GH_total_abundance,
    `Ranked cheating index` = rank_ratio_global,
    `Genome size` = genome_size
  ) %>%
  pivot_longer(
    cols = -c(genome_id, label, highlight),
    names_to = "metric",
    values_to = "value"
  ) %>%
  filter(is.finite(value))

df_long$metric <- factor(df_long$metric, levels = c(
  "Complex carbohydrate depolymerization",
  "Carbohydrate transporters",
  "Ranked cheating index",
  "Genome size"
))

# ---------------------------
# plot: POINTS ONLY + median marker (NO violin)
# ---------------------------
# Optional: jitter width per panel (a bit smaller for ranked)
# We'll keep one jitter width that works generally.
p <- ggplot(df_long, aes(x = value, y = label)) +
  
  # grey points first
  geom_jitter(
    data = df_long %>% filter(highlight == "Other"),
    color = "grey75",
    width = 0.18,
    size = 1.2,
    alpha = 0.55
  ) +
  
  # green points on top (Top 10%)
  geom_jitter(
    data = df_long %>% filter(highlight == "Top10"),
    color = "#1b9e77",
    width = 0.18,
    size = 1.4,
    alpha = 0.55
  ) +
  
  # median marker per phylum within each panel
  stat_summary(
    fun = median,
    geom = "point",
    shape = 124,   # vertical tick
    size = 9,
    color = "black"
  ) +
  
  facet_wrap(~ metric, nrow = 1, scales = "free_x") +
  
  labs(
    x = NULL,
    y = NULL,
    title = "Phylum distributions (Top 10% global ranked cheating index highlighted)"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 12),
    axis.text.y = element_text(face = "bold"),
    panel.spacing.x = unit(0.8, "lines"),
    plot.margin = margin(10, 10, 10, 10),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  coord_cartesian(clip = "off")

print(p)

# ---------------------------
# save
# ---------------------------
ggsave(out_pdf, p, width = 18, height = 9)
ggsave(out_png, p, width = 18, height = 9, dpi = 300)

message("Saved:")
message(" - ", out_pdf)
message(" - ", out_png)
message("Done ✅")
