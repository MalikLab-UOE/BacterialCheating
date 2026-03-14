#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
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

# unify columns
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
      ABC_total_abundance = suppressWarnings(as.numeric(ABC_total_abundance)),
      GH_total_abundance  = suppressWarnings(as.numeric(GH_total_abundance)),
      ABC_GH_ratio        = ratio
    ) %>%
    dplyr::select(
      genome_id,
      ABC_total_abundance,
      GH_total_abundance,
      ABC_GH_ratio,
      database
    )
}

# ---------------------------
# inputs
# ---------------------------
img_rds  <- "data/processed/img_input.rds"
loma_rds <- "data/processed/loma_input.rds"
iso_rds  <- "data/processed/isolates_input.rds"
glo_rds  <- "data/processed/glo_input.rds"

out_png  <- "results/figures/ALL_scatter_rank_global_facet_binary.png"

stop_if_missing(img_rds,  "IMG rds")
stop_if_missing(loma_rds, "LOMA rds")
stop_if_missing(iso_rds,  "isolates rds")
stop_if_missing(glo_rds,  "GLO rds")

# ---------------------------
# load + tag
# ---------------------------
img  <- readRDS(img_rds)  %>% mutate(database = "IMG")
loma <- readRDS(loma_rds) %>% mutate(database = "LOMA")
iso  <- readRDS(iso_rds)  %>% mutate(database = "isolates")
glo  <- readRDS(glo_rds)  %>% mutate(database = "GLO")

img  <- unify_cols(img)
loma <- unify_cols(loma)
iso  <- unify_cols(iso)
glo  <- unify_cols(glo)

# ---------------------------
# merge + clean
# ---------------------------
all_df <- bind_rows(loma, img, iso, glo) %>%
  filter(
    is.finite(ABC_GH_ratio),
    is.finite(ABC_total_abundance),
    is.finite(GH_total_abundance)
  )

if (nrow(all_df) == 0) {
  stop("No rows left after filtering. Check ratio / abundance columns.", call. = FALSE)
}

# ---------------------------
# compute global rank
# ---------------------------
db_levels <- c("LOMA", "IMG", "isolates", "GLO")

all_df <- all_df %>%
  mutate(
    rank_ratio_global = percent_rank(ABC_GH_ratio),
    database = factor(database, levels = db_levels),
    highlight = ifelse(rank_ratio_global >= 0.9,
                       "Top 10% Ranked cheating index",
                       "Other genomes")
  )

message("Rows per dataset:")
print(all_df %>% count(database, highlight))

# ---------------------------
# plot
# ---------------------------
p_rank_facet <- ggplot(
  all_df,
  aes(
    x = GH_total_abundance,
    y = ABC_total_abundance,
    color = highlight
  )
) +
  geom_point(size = 2, alpha = 0.9) +
  scale_color_manual(
    values = c(
      "Other genomes" = "grey80",
      "Top 10% Ranked cheating index" = "#1b9e77"
    ),
    name = NULL
  ) +
  facet_wrap(~ database, ncol = 2, nrow = 2, scales = "free") +
  labs(
    x = "Glycoside hydrolases",
    y = "Carbohydrate transporters",
    title = "Transport vs depolymerization (Top 10% Ranked cheating index highlighted)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 11),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )

print(p_rank_facet)

make_dir(dirname(out_png))
ggsave(out_png, p_rank_facet, width = 10, height = 10, dpi = 300)

message("Saved: ", out_png)
message("Done ✅")