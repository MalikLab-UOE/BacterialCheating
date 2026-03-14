#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(scales)
  library(tools)
})

# =========================================================
# Helpers
# =========================================================
stop_if_missing <- function(path, label = NULL) {
  if (!file.exists(path)) {
    msg <- if (!is.null(label)) paste0(label, " not found: ") else "File not found: "
    stop(msg, path, call. = FALSE)
  }
}
make_dir <- function(dir) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}

# unify columns across datasets (genome_size here, not mgt)
unify_cols <- function(df) {
  ratio <- NA_real_
  if ("ABC_GH_ratio" %in% names(df)) {
    ratio <- suppressWarnings(as.numeric(df$ABC_GH_ratio))
  } else if ("trans_degra" %in% names(df)) {
    ratio <- suppressWarnings(as.numeric(df$trans_degra))
  }
  
  df %>%
    mutate(
      taxon_oid = if ("taxon_oid" %in% names(.)) as.character(taxon_oid) else NA_character_,
      id        = if ("id"        %in% names(.)) as.character(id)        else NA_character_,
      Bin.ID    = if ("Bin.ID"    %in% names(.)) as.character(`Bin.ID`)  else NA_character_,
      genome_id = coalesce(taxon_oid, id, Bin.ID),
      
      genome_size = if ("genome_size" %in% names(.))
        suppressWarnings(as.numeric(.data$genome_size)) else NA_real_,
      
      ABC_total_abundance = suppressWarnings(as.numeric(ABC_total_abundance)),
      GH_total_abundance  = suppressWarnings(as.numeric(GH_total_abundance)),
      ABC_GH_ratio        = ratio
    ) %>%
    dplyr::select(
      genome_id,
      genome_size,
      ABC_total_abundance,
      GH_total_abundance,
      ABC_GH_ratio,
      database
    )
}

# =========================================================
# Inputs / outputs
# =========================================================
img_rds  <- "data/processed/img_input.rds"
loma_rds <- "data/processed/loma_input.rds"
iso_rds  <- "data/processed/isolates_input.rds"
glo_rds  <- "data/processed/glo_input.rds"

out_dir <- "results/figures"
out_pdf <- file.path(out_dir, "Fig_S4.pdf")
out_png <- file.path(out_dir, "Fig_S4.png")

stop_if_missing(img_rds,  "IMG rds")
stop_if_missing(loma_rds, "LOMA rds")
stop_if_missing(iso_rds,  "isolates rds")
stop_if_missing(glo_rds,  "GLO rds")
make_dir(out_dir)

# =========================================================
# 1) Load + tag
# =========================================================
img  <- readRDS(img_rds)  %>% mutate(database = "IMG")
loma <- readRDS(loma_rds) %>% mutate(database = "LOMA")
iso  <- readRDS(iso_rds)  %>% mutate(database = "isolates")
glo  <- readRDS(glo_rds)  %>% mutate(database = "GLO")

img  <- unify_cols(img)
loma <- unify_cols(loma)
iso  <- unify_cols(iso)
glo  <- unify_cols(glo)

# =========================================================
# 2) Merge + clean global pool (for ranking)
# =========================================================
db_levels <- c("LOMA", "IMG", "isolates", "GLO")

all_df <- bind_rows(loma, img, iso, glo) %>%
  mutate(database = factor(database, levels = db_levels)) %>%
  filter(
    is.finite(ABC_GH_ratio),
    is.finite(ABC_total_abundance),
    is.finite(GH_total_abundance)
  )

if (nrow(all_df) == 0) stop("No rows after filtering ratio/abundance.", call. = FALSE)

# =========================================================
# 3) Global rank computed ONCE across all datasets
# =========================================================
all_df <- all_df %>%
  mutate(rank_ratio_global = percent_rank(ABC_GH_ratio))

message("Rows per dataset (global pool):")
print(all_df %>% count(database))

# =========================================================
# 4) Subset to GLO and keep genome_size
# =========================================================
plot_df <- all_df %>%
  filter(database == "GLO") %>%
  filter(is.finite(genome_size), genome_size > 0) %>%
  mutate(
    cheating_class = ifelse(rank_ratio_global >= 0.9,
                            "Ranked cheating index ≥ 0.9 (Top 10%)",
                            "Ranked cheating index < 0.9")
  )

if (nrow(plot_df) == 0) {
  stop("No rows left after filtering GLO and genome_size.", call. = FALSE)
}

plot_df$cheating_class <- factor(
  plot_df$cheating_class,
  levels = c("Ranked cheating index < 0.9", "Ranked cheating index ≥ 0.9 (Top 10%)"),
  ordered = TRUE
)

cat("Counts per class:\n")
print(table(plot_df$cheating_class, useNA = "ifany"))

# =========================================================
# 5) Significance test (Wilcoxon)
# =========================================================
wilcox_res <- wilcox.test(genome_size ~ cheating_class, data = plot_df)
p_value <- wilcox_res$p.value

p_label <- ifelse(
  is.na(p_value), "p = NA",
  ifelse(p_value < 0.001, "p < 0.001", paste0("p = ", signif(p_value, 3)))
)

x_max <- max(plot_df$genome_size, na.rm = TRUE)

# =========================================================
# 6) Plot: violin + points colored by group + median diamond + p label
# =========================================================
FILL_COLS <- c(
  "Ranked cheating index < 0.9" = "grey85",
  "Ranked cheating index ≥ 0.9 (Top 10%)" = "#1b9e77"
)

POINT_COLS <- c(
  "Ranked cheating index < 0.9" = "grey65",
  "Ranked cheating index ≥ 0.9 (Top 10%)" = "#66c2a4"
)

p <- ggplot(plot_df, aes(y = cheating_class, x = genome_size)) +
  
  geom_violin(
    aes(fill = cheating_class),
    trim = TRUE,
    scale = "width",
    color = "grey30",
    linewidth = 0.55,
    alpha = 0.85
  ) +
  
  geom_point(
    aes(color = cheating_class),
    position = position_jitter(height = 0.15, width = 0),
    size = 0.7,
    alpha = 0.18
  ) +
  
  stat_summary(
    fun = median,
    geom = "point",
    shape = 23,
    size = 3.0,
    fill = "white",
    color = "black",
    stroke = 0.6
  ) +
  
  scale_fill_manual(values = FILL_COLS) +
  scale_color_manual(values = POINT_COLS) +
  
  annotate(
    "text",
    x = x_max * 0.95,
    y = 1.5,
    label = p_label,
    size = 5,
    fontface = "bold"
  ) +
  
  scale_x_continuous(labels = scales::label_comma()) +
  
  labs(
    y = NULL,
    x = "Genome size (bp)",
    title = "Genome size by GLOBAL ranked cheating index\n(GLO; Wilcoxon test)"
  ) +
  
  theme_bw(base_size = 15) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 17, margin = margin(b = 8)),
    axis.title.x = element_text(face = "bold", margin = margin(t = 8)),
    axis.text = element_text(color = "black", size = 12),
    legend.position = "none"
  )

print(p)

# =========================================================
# 7) Save
# =========================================================
ggsave(out_pdf, p, width = 7.5, height = 4.8)
ggsave(out_png, p, width = 7.5, height = 4.8, dpi = 300)

message("Saved:")
message(" - ", out_pdf)
message(" - ", out_png)
message("Done ✅")