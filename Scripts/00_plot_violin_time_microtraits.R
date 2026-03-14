#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tools)
})

# ------------------
# inputs (edit)
# ------------------
in_csv  <- "data/processed/loma_input.csv"
out_pdf <- "results/figures/Raw_cheating_index_across_Families_top30.pdf"
out_png <- "results/figures/Raw_cheating_index_across_Families_top30.png"
out_dir <- dirname(out_pdf)

# ------------------
# params (keep your selection rule)
# ------------------
TOP_N_FAMILY  <- 30
MIN_N_FAMILY  <- 5

time_cols <- c(
  "T4" = "#0571b0",
  "T3" = "#92c5de",
  "T2" = "#f4a582",
  "T1" = "#ca0020"
)

# ------------------
# helpers
# ------------------
stop_if_missing <- function(path, label = NULL) {
  if (!file.exists(path)) {
    msg <- if (!is.null(label)) paste0(label, " not found: ") else "File not found: "
    stop(msg, path, call. = FALSE)
  }
}
make_dir <- function(dir) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}

# ------------------
# load
# ------------------
stop_if_missing(in_csv, "Input CSV")
make_dir(out_dir)

merged_table <- readr::read_csv(in_csv, show_col_types = FALSE)

# ------------------
# choose column for Raw cheating index
# ------------------
# 兼容两种情况：
# 1) 你的数据仍是 ABC_GH_ratio
# 2) 你把列名改成 "Raw cheating index"
raw_col <- if ("Raw cheating index" %in% names(merged_table)) {
  "Raw cheating index"
} else if ("ABC_GH_ratio" %in% names(merged_table)) {
  "ABC_GH_ratio"
} else {
  stop("Cannot find Raw cheating index column. Need either 'ABC_GH_ratio' or 'Raw cheating index'.", call. = FALSE)
}

# ------------------
# checks
# ------------------
need_cols <- c("Family", "time", raw_col)
miss <- setdiff(need_cols, names(merged_table))
if (length(miss) > 0) {
  stop("Missing columns in input CSV: ", paste(miss, collapse = ", "), call. = FALSE)
}

# ------------------
# clean
# ------------------
df <- merged_table %>%
  mutate(
    Family = as.character(.data$Family),
    time   = as.character(.data$time),
    Raw_cheating_index = suppressWarnings(as.numeric(.data[[raw_col]]))
  ) %>%
  filter(
    !is.na(Family), Family != "",
    !is.na(time), time != "",
    is.finite(Raw_cheating_index)
  ) %>%
  mutate(
    time = toupper(trimws(time)),
    time = factor(time, levels = c("T1", "T2", "T3", "T4"))
  ) %>%
  filter(!is.na(time))

# ------------------
# Family order (your rule)
# ------------------
family_order <- df %>%
  group_by(Family) %>%
  filter(n() >= MIN_N_FAMILY) %>%   # keep your rule: n() >= 5
  summarise(
    median_index = median(Raw_cheating_index, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median_index)) %>%
  slice_head(n = TOP_N_FAMILY) %>%
  pull(Family)

if (length(family_order) == 0) {
  stop("No families passed filter: n() >= ", MIN_N_FAMILY, " with finite Raw cheating index.", call. = FALSE)
}

df_plot <- df %>%
  filter(Family %in% family_order) %>%
  mutate(Family = factor(Family, levels = rev(family_order)))

# ------------------
# plot
# ------------------
p <- ggplot(df_plot, aes(x = Family, y = Raw_cheating_index)) +
  geom_violin(
    width = 0.8,
    fill  = "grey90",
    color = "white",
    scale = "width",
    trim  = FALSE
  ) +
  geom_jitter(
    aes(color = time),
    width = 0.2, size = 2, alpha = 0.8
  ) +
  stat_summary(
    fun = median, geom = "crossbar",
    width = 0.5, color = "black", linewidth = 0.3
  ) +
  scale_color_manual(values = time_cols, name = "Time", drop = FALSE) +
  labs(
    x = NULL,
    y = "Raw cheating index",
    title = paste0("Raw cheating index across Families (Top ", TOP_N_FAMILY, " by median, n≥", MIN_N_FAMILY, ")")
  ) +
  coord_flip() +
  theme_classic(base_size = 15) +
  theme(
    axis.text = element_text(color = "black"),
    axis.text.y = element_text(size = 12),
    axis.line = element_line(linewidth = 1),
    axis.ticks = element_line(linewidth = 1),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    plot.margin = margin(10, 10, 10, 10)
  )

print(p)

# ------------------
# save
# ------------------
ggsave(out_pdf, plot = p, width = 6, height = 10)
ggsave(out_png, plot = p, width = 6, height = 10, dpi = 300)

message("Saved:")
message(" - ", out_pdf)
message(" - ", out_png)
message("Done ✅")