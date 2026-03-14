#!/usr/bin/env Rscript

# =========================
# Full script: bar(mean)+error(SE)+dots + significance + colors
#  - bars translucent
#  - significance stars capped at ***
# =========================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(rstatix)
  library(ggpubr)
  library(readr)
})

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
# INPUT / OUTPUT
# -------------------------
in_csv  <- "Input_Data/isolates_group.csv"
out_dir <- "results/figures"

out_pdf <- file.path(out_dir, "Fig_S1.pdf")
out_png <- file.path(out_dir, "Fig_S1.png")

stop_if_missing(in_csv, "isolates_group.csv")
make_dir(out_dir)

# -------------------------
# READ INPUT
# -------------------------
groups <- readr::read_csv(in_csv, show_col_types = FALSE)

# -------------------------
# REQUIRED columns in `groups`
# -------------------------
need_cols <- c("member", "Abundance")
miss <- setdiff(need_cols, names(groups))
if (length(miss) > 0) {
  stop("Missing columns in groups: ", paste(miss, collapse = ", "), call. = FALSE)
}

# -------------------------
# Basic cleaning
# -------------------------
groups <- groups %>%
  mutate(
    member = as.factor(member),
    Abundance = suppressWarnings(as.numeric(Abundance))
  ) %>%
  filter(!is.na(member), is.finite(Abundance))

# Optional: set order
# groups$member <- factor(groups$member, levels = c("B", "BF"))

# -------------------------
# Summary stats (mean + SE)
# -------------------------
groups_sum <- groups %>%
  group_by(member) %>%
  summarise(
    n = sum(!is.na(Abundance)),
    mean_abundance = mean(Abundance, na.rm = TRUE),
    sd_abundance   = sd(Abundance, na.rm = TRUE),
    se_abundance   = sd_abundance / sqrt(n),
    .groups = "drop"
  )

# -------------------------
# Significance test
#   2 groups -> t-test
#   >=3     -> ANOVA + Tukey HSD
# -------------------------
k <- dplyr::n_distinct(groups$member)

if (k == 2) {
  stat_res <- groups %>%
    t_test(Abundance ~ member) %>%
    mutate(p = as.numeric(p)) %>%
    mutate(
      group1 = as.character(group1),
      group2 = as.character(group2),
      y.position = max(groups$Abundance, na.rm = TRUE) * 1.08
    )
} else if (k >= 3) {
  tuk <- groups %>%
    tukey_hsd(Abundance ~ member) %>%
    mutate(
      group1 = as.character(group1),
      group2 = as.character(group2),
      p = as.numeric(p.adj)
    )
  
  y0 <- max(groups$Abundance, na.rm = TRUE)
  step <- (y0 - min(groups$Abundance, na.rm = TRUE)) * 0.08
  if (!is.finite(step) || step == 0) step <- y0 * 0.05 + 0.05
  
  stat_res <- tuk %>%
    mutate(y.position = y0 + step * row_number())
} else {
  stop("Need at least 2 groups in `member`.", call. = FALSE)
}

# -------------------------
# Cap stars at *** (no ****)
# -------------------------
stat_res <- stat_res %>%
  mutate(
    p.signif = case_when(
      p <= 0.001 ~ "***",
      p <= 0.01  ~ "**",
      p <= 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )

# -------------------------
# Colors (optional fixed palette)
# -------------------------
my_cols <- NULL
# my_cols <- c("B" = "#F4A7A1", "BF" = "#8DD3C7")

# -------------------------
# Plot
# -------------------------
p <- ggplot() +
  geom_col(
    data = groups_sum,
    aes(x = member, y = mean_abundance, fill = member),
    width = 0.65,
    color = "black",
    alpha = 0.30
  ) +
  geom_errorbar(
    data = groups_sum,
    aes(
      x = member,
      ymin = mean_abundance - se_abundance,
      ymax = mean_abundance + se_abundance
    ),
    width = 0.2,
    linewidth = 0.6
  ) +
  geom_jitter(
    data = groups,
    aes(x = member, y = Abundance, color = member),
    width = 0.15,
    height = 0,
    size = 2,
    alpha = 0.65
  ) +
  stat_pvalue_manual(
    stat_res,
    label = "p.signif",
    tip.length = 0.01,
    hide.ns = TRUE
  ) +
  labs(
    x = "Member",
    y = "Abundance"
  ) +
  theme_classic(base_size = 14) +
  guides(color = "none")

if (!is.null(my_cols)) {
  p <- p +
    scale_fill_manual(values = my_cols) +
    scale_color_manual(values = my_cols)
}

print(p)

# -------------------------
# Save
# -------------------------
ggsave(out_pdf, p, width = 6, height = 4, dpi = 300)
ggsave(out_png, p, width = 6, height = 4, dpi = 300)

message("Saved:")
message(" - ", out_pdf)
message(" - ", out_png)
message("Done ✅")