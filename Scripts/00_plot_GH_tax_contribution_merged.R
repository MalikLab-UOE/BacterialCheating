#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(rstatix)
})

# ---------------------------
# helpers
# ---------------------------
stop_if_missing <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path, call. = FALSE)
}
make_dir <- function(dir) if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)

# ---------------------------
# args
# ---------------------------
args <- commandArgs(trailingOnly = TRUE)
in_csv  <- ifelse(length(args) >= 1, args[1], "Input_Data/GH_tax_contribuction.csv")
out_dir <- ifelse(length(args) >= 2, args[2], "results/figures")

stop_if_missing(in_csv)
make_dir(out_dir)

out_pdf <- file.path(out_dir, "Fig_1_GH_TPM_Contribution_tax_mergedEco_ab_points.pdf")

# ---------------------------
# 1) read + clean
# ---------------------------
final_df <- readr::read_csv(in_csv, show_col_types = FALSE)

need <- c("Sample", "Group", "tax", "GH_FPKM_Contribution")
miss <- setdiff(need, names(final_df))
if (length(miss) > 0) {
  stop("Missing columns: ", paste(miss, collapse = ", "), call. = FALSE)
}

data_clean <- final_df %>%
  filter(!is.na(tax), tax %in% c("Bacteria", "Fungi")) %>%
  mutate(
    tax = factor(tax, levels = c("Bacteria", "Fungi")),
    GH_FPKM_Contribution = as.numeric(GH_FPKM_Contribution)
  ) %>%
  filter(is.finite(GH_FPKM_Contribution))

# ---------------------------
# 2) t-test (Bacteria vs Fungi)
# ---------------------------
tt <- t_test(data_clean, GH_FPKM_Contribution ~ tax)

letters_df <- if (tt$p[1] < 0.05) {
  tibble(
    tax = factor(c("Bacteria", "Fungi"), levels = c("Bacteria", "Fungi")),
    label = c("a", "b")
  )
} else {
  tibble(
    tax = factor(c("Bacteria", "Fungi"), levels = c("Bacteria", "Fungi")),
    label = c("a", "a")
  )
}

# ---------------------------
# 3) mean ± SE (Eco merged)
# ---------------------------
plot_data <- data_clean %>%
  group_by(tax) %>%
  summarise(
    Mean_FPKM = mean(GH_FPKM_Contribution, na.rm = TRUE),
    SE_FPKM   = sd(GH_FPKM_Contribution, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  left_join(letters_df, by = "tax") %>%
  mutate(
    y_pos = Mean_FPKM + SE_FPKM +
      0.05 * max(Mean_FPKM + SE_FPKM, na.rm = TRUE)
  )

# ---------------------------
# 4) plot
# ---------------------------
cols <- c("Bacteria" = "#ff6e69", "Fungi" = "#00c3c4")

p <- ggplot(plot_data, aes(x = tax, y = Mean_FPKM, fill = tax)) +
  geom_col(color = "black", width = 0.7, alpha = 0.9) +
  geom_jitter(
    data = data_clean,
    aes(x = tax, y = GH_FPKM_Contribution, color = tax),
    width = 0.15, size = 2.2, alpha = 0.6,
    show.legend = FALSE, inherit.aes = FALSE
  ) +
  geom_errorbar(
    aes(ymin = Mean_FPKM - SE_FPKM, ymax = Mean_FPKM + SE_FPKM),
    width = 0.15, color = "black"
  ) +
  geom_text(
    aes(y = y_pos, label = label),
    vjust = 0, fontface = "bold", size = 5
  ) +
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  labs(
    title = "GH_FPKM Contribution (Eco merged)",
    subtitle = paste0("t-test p = ", signif(tt$p[1], 3)),
    x = NULL,
    y = "Mean TPM (± SE)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.text.x  = element_text(size = 12, color = "black"),
    axis.text.y  = element_text(size = 12, color = "black"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(p)

ggsave(out_pdf, p, width = 5, height = 3, dpi = 300)

message("Saved: ", out_pdf)
message("Done ✅")