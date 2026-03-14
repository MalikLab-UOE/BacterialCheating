#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(rstatix)
  library(tibble)
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

# ---------------------------
# args
# ---------------------------
args <- commandArgs(trailingOnly = TRUE)

in_csv  <- ifelse(length(args) >= 1, args[1], "Input_Data/LOMA_merged_mags_50_10.csv")
out_pdf <- ifelse(length(args) >= 2, args[2], "results/figures/GH_violin_taxonomy_ab_medianline.pdf")

# which groups (edit if your file uses different labels)
tax_levels <- c("Bacteria", "euk")

# colors
fill_cols <- c("Bacteria" = "#E76F51", "euk" = "#2A9D8F")

# ---------------------------
# checks + load
# ---------------------------
stop_if_missing(in_csv, "Input CSV")
make_dir(dirname(out_pdf))

df <- readr::read_csv(in_csv, show_col_types = FALSE)

need <- c("tax", "GH")
miss <- setdiff(need, names(df))
if (length(miss) > 0) {
  stop("Missing columns in input: ", paste(miss, collapse = ", "), call. = FALSE)
}

# ---------------------------
# clean
# ---------------------------
df_clean <- df %>%
  mutate(
    tax = as.character(tax),
    GH  = suppressWarnings(as.numeric(GH))
  ) %>%
  filter(
    !is.na(tax),
    !is.na(GH),
    is.finite(GH),
    tax %in% tax_levels
  ) %>%
  mutate(
    tax = factor(tax, levels = tax_levels)
  )

if (nrow(df_clean) == 0) {
  stop("No rows left after filtering. Check tax labels and GH column.", call. = FALSE)
}

# ---------------------------
# stats: Wilcoxon + letters
# ---------------------------
wilc <- rstatix::wilcox_test(df_clean, GH ~ tax)

pval <- wilc$p[1]
letters_df <- if (is.finite(pval) && pval < 0.05) {
  tibble(tax = tax_levels, label = c("a", "b"))
} else {
  tibble(tax = tax_levels, label = c("a", "a"))
}

# label position: max + 5%
y_top <- df_clean %>%
  group_by(tax) %>%
  summarise(y = max(GH, na.rm = TRUE), .groups = "drop") %>%
  left_join(letters_df, by = "tax") %>%
  mutate(y = y * 1.05)

# ---------------------------
# plot
# ---------------------------
p <- ggplot(df_clean, aes(x = tax, y = GH, fill = tax)) +
  geom_violin(
    color = "black",
    linewidth = 0.5,
    trim = FALSE,
    scale = "width",
    alpha = 0.45
  ) +
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.45,
    color = "black",
    fill = NA,
    fatten = 1.3
  ) +
  geom_jitter(
    aes(color = tax),
    width = 0.2,
    size = 1.8,
    alpha = 0.7
  ) +
  geom_text(
    data = y_top,
    aes(x = tax, y = y, label = label),
    inherit.aes = FALSE,
    vjust = 0,
    fontface = "bold",
    size = 5
  ) +
  scale_fill_manual(values = fill_cols) +
  scale_color_manual(values = fill_cols) +
  labs(
    x = NULL,
    y = "GH",
    title = "Violin plot of GH by taxonomy",
    subtitle = paste0("Wilcoxon p = ", signif(pval, 3))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title  = element_text(size = 14, face = "bold"),
    panel.grid.major = element_line(color = "gray80", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )



print(p)
ggsave(out_pdf, p, width = 5, height = 3, dpi = 300)

message("Saved: ", out_pdf)
message("Rows used: ", nrow(df_clean))
message("Done ✅")