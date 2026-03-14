#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
  library(tools)
  library(multcompView)  # Tukey letters
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
# args / paths
# ---------------------------
args <- commandArgs(trailingOnly = TRUE)

in_tsv     <- ifelse(length(args) >= 1, args[1], "Input_Data/LOMA_FBratio.txt")
out_prefix <- ifelse(length(args) >= 2, args[2], "Fig_1_FB")
out_dir    <- ifelse(length(args) >= 3, args[3], "results/figures")

out_pdf_bf <- file.path(out_dir, paste0(out_prefix, "_bac_fun_bar_points.pdf"))
out_png_bf <- file.path(out_dir, paste0(out_prefix, "_bac_fun_bar_points.png"))

out_pdf_fb <- file.path(out_dir, paste0(out_prefix, "_ratio_fb_Tukey_boxplot.pdf"))
out_png_fb <- file.path(out_dir, paste0(out_prefix, "_ratio_fb_Tukey_boxplot.png"))

message("=== Plot bac/fun + fb ratio ===")
message("Input: ", in_tsv)
message("Output dir: ", out_dir)

stop_if_missing(in_tsv, "Input table")
make_dir(out_dir)

# ---------------------------
# load + validate
# ---------------------------
dat <- readr::read_delim(in_tsv, delim = "\t", show_col_types = FALSE)
names(dat) <- trimws(names(dat))

need <- c("samp", "bac", "fun")
miss <- setdiff(need, names(dat))
if (length(miss) > 0) stop("Missing columns: ", paste(miss, collapse = ", "), call. = FALSE)

# ---------------------------
# clean + compute ratio
# ---------------------------
dat2 <- dat %>%
  mutate(
    samp = factor(samp, levels = c("T1","T2","T3","T4")),
    bac  = suppressWarnings(as.numeric(bac)),
    fun  = suppressWarnings(as.numeric(fun))
  ) %>%
  filter(!is.na(samp), is.finite(bac), is.finite(fun)) %>%
  mutate(
    fb = fun / bac,
    fb = if_else(is.finite(fb), fb, NA_real_)
  )

if (nrow(dat2) == 0) stop("0 rows after filtering. Check samp/bac/fun.", call. = FALSE)

# =========================================================
# A) Figure 1: bac/fun bar (mean ± SEM) + jitter points
# =========================================================
long <- dat2 %>%
  pivot_longer(cols = c(bac, fun), names_to = "group", values_to = "abundance") %>%
  mutate(
    group = recode(group, bac = "Bacteria", fun = "Fungi"),
    group = factor(group, levels = c("Bacteria","Fungi"))
  )

sum_df <- long %>%
  group_by(samp, group) %>%
  summarise(
    n = dplyr::n(),
    mean = mean(abundance, na.rm = TRUE),
    sd   = sd(abundance, na.rm = TRUE),
    sem  = sd / sqrt(n),
    .groups = "drop"
  )

fill_cols_bf <- c("Bacteria" = "#F4A3A3", "Fungi" = "#59DAD8")

p_bf <- ggplot(sum_df, aes(x = samp, y = mean, fill = group)) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.68,
    color = "black",
    linewidth = 0.6,
    alpha = 0.85
  ) +
  geom_errorbar(
    aes(ymin = mean - sem, ymax = mean + sem),
    position = position_dodge(width = 0.75),
    width = 0.20,
    linewidth = 0.6
  ) +
  geom_point(
    data = long,
    aes(x = samp, y = abundance, color = group),
    position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75),
    size = 1.6,
    alpha = 0.65,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = fill_cols_bf) +
  scale_color_manual(values = fill_cols_bf) +
  labs(x = NULL, y = "Abundance", fill = NULL, color = NULL) +
  theme_classic(base_size = 14) +
  theme(axis.text = element_text(color = "black"))

print(p_bf)
ggsave(out_pdf_bf, p_bf, width = 6.2, height = 4.5, dpi = 600)
ggsave(out_png_bf, p_bf, width = 6.2, height = 4.5, dpi = 600)

# =========================================================
# B) Figure 1: fb ratio boxplot + Tukey letters
# =========================================================
data_clean <- dat2 %>%
  filter(is.finite(fb)) %>%
  droplevels()

if (nrow(data_clean) < 3) stop("Not enough rows for fb ANOVA after filtering.", call. = FALSE)

# 颜色：按时间点
my_colors <- c('T4' = '#0571b0', 'T3' = '#92c5de', 'T2' = '#f4a582', 'T1' = '#ca0020')

aov_mod   <- aov(fb ~ samp, data = data_clean)
tukey_res <- TukeyHSD(aov_mod)
letters   <- multcompLetters4(aov_mod, tukey_res)

letters_df <- tibble::tibble(
  samp  = names(letters$samp$Letters),
  label = letters$samp$Letters
)

ann_df <- data_clean %>%
  group_by(samp) %>%
  summarise(y_pos = max(fb, na.rm = TRUE) * 1.08, .groups = "drop") %>%
  left_join(letters_df, by = "samp") %>%
  mutate(samp = factor(samp, levels = levels(data_clean$samp)))

p_fb <- ggplot(data_clean, aes(x = samp, y = fb, color = samp)) +
  geom_boxplot(
    outlier.shape = 21,
    outlier.fill = "white",
    outlier.size = 2,
    linewidth = 1.1,
    alpha = 0.8
  ) +
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.6) +
  geom_text(
    data = ann_df,
    aes(x = samp, y = y_pos, label = label, color = samp),
    inherit.aes = FALSE,
    vjust = 0,
    fontface = "bold",
    size = 5
  ) +
  scale_color_manual(values = my_colors) +
  labs(
    x = NULL,
    y = "Fungi:Bacteria Ratio (fb = fun/bac)",
    title = "fb across sampling time (Tukey groups)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text = element_text(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2)
  )

print(p_fb)
ggsave(out_pdf_fb, p_fb, width = 5, height = 5, dpi = 600)
ggsave(out_png_fb, p_fb, width = 5, height = 5, dpi = 600)

message("Saved:")
message(" - ", out_pdf_bf)
message(" - ", out_png_bf)
message(" - ", out_pdf_fb)
message(" - ", out_png_fb)
message("Done ✅")