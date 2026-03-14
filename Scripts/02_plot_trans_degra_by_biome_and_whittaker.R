#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(rstatix)
  library(multcompView)
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

# ---------------------------
# Args / paths
# ---------------------------
args <- commandArgs(trailingOnly = TRUE)

in_rds     <- ifelse(length(args) >= 1, args[1], "data/processed/img_input.rds")
fig_prefix <- ifelse(length(args) >= 2, args[2], "Fig_3")   # e.g., Fig2/Fig3/SFig1
min_n      <- ifelse(length(args) >= 3, as.integer(args[3]), 100)

fig_dir   <- "results/figures"
table_dir <- "results/tables"

out_box_pdf  <- file.path(fig_dir,   paste0(fig_prefix, "_trans_degra_by_biome_log2_box.pdf"))
out_whit_pdf <- file.path(fig_dir,   paste0(fig_prefix, "_trans_degra_by_biome_RAWmedian_whittaker.pdf"))

out_cld_csv  <- file.path(table_dir, paste0(fig_prefix, "_trans_degra_by_biome_CLD_log2.csv"))
out_med_csv  <- file.path(table_dir, paste0(fig_prefix, "_trans_degra_by_biome_medians_RAW.csv"))

message("=== Plot trans/degra by biome + Whittaker annotation ===")
message("Input RDS:      ", in_rds)
message("Fig prefix:     ", fig_prefix)
message("min_n:          ", min_n)
message("Boxplot PDF:    ", out_box_pdf)
message("Whittaker PDF:  ", out_whit_pdf)
message("CLD CSV:        ", out_cld_csv)
message("RAW medians:    ", out_med_csv)

stop_if_missing(in_rds, "Input RDS")
make_dir(fig_dir)
make_dir(table_dir)

# ---------------------------
# Load data
# ---------------------------
merged_table <- readRDS(in_rds)

required_cols <- c("Biome_new", "ABC_GH_ratio")
missing <- setdiff(required_cols, colnames(merged_table))
if (length(missing) > 0) {
  stop("Missing required columns in merged_table: ", paste(missing, collapse = ", "), call. = FALSE)
}

# ---------------------------
# Colors: Whittaker
# ---------------------------
whittaker_cols <- c(
  "Tundra"                           = "#C1E1DD",
  "Boreal forest"                    = "#A5C790",
  "Temperate seasonal forest"        = "#97B669",
  "Temperate rain forest"            = "#75A95E",
  "Tropical rain forest"             = "#317A22",
  "Tropical seasonal forest/savanna" = "#A09700",
  "Subtropical desert"               = "#DCBB50",
  "Temperate grassland/desert"       = "#FCD57A",
  "Woodland/shrubland"               = "#D16E3F"
)

# ---------------------------
# 1) Preprocess + filter groups by n
# ---------------------------
plot_table <- merged_table %>%
  rename(Biome_label_raw = Biome_new) %>%
  filter(!is.na(Biome_label_raw), is.finite(ABC_GH_ratio)) %>%
  mutate(
    # for stats/boxplot
    trans_degra_log2 = log2(pmax(ABC_GH_ratio, .Machine$double.eps)),
    # for Whittaker annotation (raw)
    trans_degra_raw = ABC_GH_ratio
  )

n_df <- plot_table %>% count(Biome_label_raw, name = "n")

plot_table <- plot_table %>%
  left_join(n_df, by = "Biome_label_raw") %>%
  filter(n >= min_n) %>%
  mutate(Biome_label = paste0(Biome_label_raw, " (n=", n, ")")) %>%
  droplevels()

if (nrow(plot_table) == 0) {
  stop("After filtering, plot_table has 0 rows. Try lowering min_n or check Biome_new / ABC_GH_ratio.", call. = FALSE)
}

# ---------------------------
# 2) RAW medians per biome (for Whittaker labels) + export
# ---------------------------
medians_raw <- plot_table %>%
  group_by(Biome_label_raw) %>%
  summarise(
    n = dplyr::n(),
    median_raw = median(trans_degra_raw, na.rm = TRUE),
    q25_raw = quantile(trans_degra_raw, 0.25, na.rm = TRUE),
    q75_raw = quantile(trans_degra_raw, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median_raw))

write.csv(medians_raw, out_med_csv, row.names = FALSE)

# For boxplot ordering we keep log2 medians (so plot looks like before)
medians_for_order <- plot_table %>%
  group_by(Biome_label) %>%
  summarise(med = median(trans_degra_log2, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(med))

plot_table <- plot_table %>%
  mutate(Biome_label = factor(Biome_label, levels = rev(medians_for_order$Biome_label))) %>%
  droplevels()

# ---------------------------
# 3) Stats: Kruskal + Dunn(BH) + CLD export (log2 scale)
# ---------------------------
kw_res <- plot_table %>% kruskal_test(trans_degra_log2 ~ Biome_label)
message("Kruskal-Wallis result:")
print(kw_res)

dunn_res <- plot_table %>%
  dunn_test(trans_degra_log2 ~ Biome_label, p.adjust.method = "BH")

pv <- dunn_res$p.adj
names(pv) <- paste(dunn_res$group1, dunn_res$group2, sep = "-")

cld_letters <- multcompLetters(pv, threshold = 0.05)$Letters

ypos_df <- plot_table %>%
  group_by(Biome_label) %>%
  summarise(
    n = dplyr::n(),
    ypos = quantile(trans_degra_log2, 0.98, na.rm = TRUE) * 1.06,
    .groups = "drop"
  )

cld_table <- data.frame(
  Biome_label = names(cld_letters),
  cld = unname(cld_letters),
  row.names = NULL
) %>%
  left_join(ypos_df, by = "Biome_label")

write.csv(cld_table, out_cld_csv, row.names = FALSE)
message("Exported CLD table: ", out_cld_csv)

# ---------------------------
# 4) Color mapping for boxplot
# ---------------------------
lab_levels  <- levels(plot_table$Biome_label)
base_biomes <- sub(" \\(n=.*\\)$", "", lab_levels)
fill_cols <- setNames(whittaker_cols[base_biomes], lab_levels)

if (any(is.na(fill_cols))) {
  warning(
    "These Biome_label values are not matched in Whittaker palette:\n",
    paste(lab_levels[is.na(fill_cols)], collapse = ", "),
    "\nConsider adding colors for them."
  )
}

# ---------------------------
# 5) Boxplot (log2; no CLD letters on plot)
# ---------------------------
set.seed(1)
p_box <- ggplot(plot_table, aes(x = Biome_label, y = trans_degra_log2)) +
  geom_point(
    aes(color = Biome_label),
    position = position_jitter(width = 0.15, height = 0),
    alpha = 0.20, size = 0.8, show.legend = FALSE
  ) +
  geom_boxplot(
    aes(fill = Biome_label),
    width = 0.6, outlier.shape = NA, alpha = 0.70, color = "black"
  ) +
  stat_summary(fun = median, geom = "point", size = 1.6, color = "black") +
  coord_flip() +
  scale_fill_manual(values = fill_cols, drop = FALSE) +
  scale_color_manual(values = fill_cols, drop = FALSE) +
  labs(
    x = NULL,
    y = expression(log[2]("trans/degra")),
    subtitle = paste0("Ordered by median (log2). Dunn (BH) CLD exported. min_n=", min_n)
  ) +
  theme_classic(base_size = 15) +
  theme(
    axis.text = element_text(color = "black"),
    axis.text.y = element_text(size = 10),
    legend.position = "none",
    plot.margin = margin(10, 60, 10, 10)
  )

print(p_box)

ggsave(
  filename = out_box_pdf,
  plot = p_box,
  width = 6,
  height = 0.8 * length(levels(plot_table$Biome_label))
)
message("Saved boxplot: ", out_box_pdf)

# ---------------------------
# 6) Whittaker base plot + annotate RAW medians
# ---------------------------
if (!requireNamespace("plotbiomes", quietly = TRUE)) {
  warning("Package 'plotbiomes' is not installed. Skipping Whittaker plot. Install with install.packages('plotbiomes').")
} else {
  suppressPackageStartupMessages(library(plotbiomes))
  data("Whittaker_biomes", package = "plotbiomes")
  
  pos_df <- Whittaker_biomes %>%
    group_by(biome) %>%
    summarise(
      x = median(temp_c, na.rm = TRUE),
      y = median(precp_cm, na.rm = TRUE),
      .groups = "drop"
    )
  
  # 你的 Biome_new 名称 ↔ plotbiomes biome 名称映射
  map_vec <- c(
    "Tundra"                           = "Tundra",
    "Boreal forest"                    = "Boreal forest",
    "Temperate seasonal forest"        = "Temperate seasonal forest",
    "Temperate rain forest"            = "Temperate rain forest",
    "Tropical rain forest"             = "Tropical rain forest",
    "Tropical seasonal forest/savanna" = "Tropical seasonal forest/savanna",
    "Subtropical desert"               = "Subtropical desert",
    "Temperate grassland/desert"       = "Temperate grassland/desert",
    "Woodland/shrubland"               = "Woodland/shrubland"
  )
  
  label_df <- medians_raw %>%
    mutate(
      biome = unname(map_vec[Biome_label_raw]),
      label = sprintf("%.2f", median_raw)   # <- RAW median（不 log）
    ) %>%
    filter(!is.na(biome)) %>%
    left_join(pos_df, by = "biome")
  
  missing_map <- setdiff(unique(medians_raw$Biome_label_raw), names(map_vec))
  if (length(missing_map) > 0) {
    warning("These Biome_label_raw are not in map_vec (won't be annotated): ",
            paste(missing_map, collapse = ", "))
  }
  
  p_whit <- whittaker_base_plot() +
    geom_text(
      data = label_df,
      aes(x = x, y = y, label = label),
      fontface = "bold",
      size = 4
    ) +
    labs(
      title = paste0(fig_prefix, ""),
      subtitle = paste0("Numbers are RAW medians (not log-transformed); min_n=", min_n),
      x = "Mean annual temperature (°C)",
      y = "Mean annual precipitation (cm)"
    )
  
  print(p_whit)
  
  ggsave(
    filename = out_whit_pdf,
    plot = p_whit,
    width = 6,
    height = 5
  )
  message("Saved Whittaker plot: ", out_whit_pdf)
}

message("Done ✅")