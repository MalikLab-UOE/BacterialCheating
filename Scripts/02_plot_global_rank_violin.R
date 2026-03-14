#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
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
# inputs / outputs
# ---------------------------
img_rds  <- "data/processed/img_input.rds"
loma_rds <- "data/processed/loma_input.rds"
iso_rds  <- "data/processed/isolates_input.rds"
glo_rds  <- "data/processed/glo_input.rds"

out_png <- "results/figures/Fig_2_global_rank_violin.png"

stop_if_missing(img_rds,  "img_input.rds")
stop_if_missing(loma_rds, "loma_input.rds")
stop_if_missing(iso_rds,  "isolates_input.rds")
stop_if_missing(glo_rds,  "glo_input.rds")
make_dir(dirname(out_png))

# ---------------------------
# load + tag dataset
# ---------------------------
img  <- readRDS(img_rds)  %>% mutate(database = "IMG")
loma <- readRDS(loma_rds) %>% mutate(database = "LOMA")
iso  <- readRDS(iso_rds)  %>% mutate(database = "isolates")
glo  <- readRDS(glo_rds)  %>% mutate(database = "glo")

# ---------------------------
# unify columns
# ---------------------------
unify_cols <- function(df){
  df %>%
    mutate(
      taxon_oid = if ("taxon_oid" %in% names(.)) as.character(taxon_oid) else NA_character_,
      id        = if ("id"        %in% names(.)) as.character(id)        else NA_character_,
      Bin.ID    = if ("Bin.ID"    %in% names(.)) as.character(`Bin.ID`)  else NA_character_,
      genome_id = coalesce(taxon_oid, id, Bin.ID),
      ABC_GH_ratio = suppressWarnings(as.numeric(ABC_GH_ratio))
    ) %>%
    dplyr::select(genome_id, ABC_GH_ratio, database)
}

all_df <- bind_rows(
  unify_cols(loma),
  unify_cols(img),
  unify_cols(iso),
  unify_cols(glo)
) %>%
  filter(is.finite(ABC_GH_ratio)) %>%
  mutate(
    database = factor(database, levels = c("LOMA", "IMG", "isolates", "glo")),
    rank_ratio_global = percent_rank(ABC_GH_ratio),
    point_alpha = ifelse(database == "glo", 0.01, 0.1)
  )

stopifnot(nrow(all_df) > 0)

# ---------------------------
# labels + counts
# ---------------------------
new_labels <- c(
  "LOMA"     = "Plant \nlitter MAGs",
  "IMG"      = "Soil MAGs from \nglobal biomes",
  "isolates" = "Soil isolate \ngenomes",
  "glo"      = "Environmental and \nhost-associated MAGs"
)

count_df <- all_df %>%
  count(database, name = "N") %>%
  mutate(
    label = paste0("N = ", N),
    y = -0.06
  )

# point colors
point_cols <- c(
  "LOMA"     = "#66C2A5",
  "IMG"      = "#8DA0CB",
  "isolates" = "#FC8D62",
  "glo"      = "#E78AC3"
)

# ---------------------------
# plot
# ---------------------------
p_fig2_violin <- ggplot(all_df, aes(x = database, y = rank_ratio_global)) +
  geom_violin(
    fill = NA,
    color = "grey35",
    linewidth = 1.0,
    trim = TRUE
  ) +
  geom_jitter(
    aes(color = database, alpha = point_alpha),
    width = 0.15,
    size = 0.7,
    show.legend = FALSE
  ) +
  scale_color_manual(values = point_cols) +
  scale_alpha_identity() +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 23,
    size = 3.0,
    fill = "white",
    color = "black",
    stroke = 0.8
  ) +
  geom_text(
    data = count_df,
    aes(x = database, y = y, label = label),
    inherit.aes = FALSE,
    size = 3.2,
    vjust = 1
  ) +
  scale_x_discrete(labels = new_labels) +
  coord_cartesian(ylim = c(-0.08, 1.05), clip = "off") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x  = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(face = "bold"),
    plot.title   = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "none",
    plot.margin = margin(10, 10, 25, 10)
  ) +
  labs(
    x = NULL,
    y = "Global rank-normalized cheating index",
    title = "Global distribution of microbial cheating index"
  )

print(p_fig2_violin)

# ---------------------------
# save
# ---------------------------
ggsave(
  filename = out_png,
  plot = p_fig2_violin,
  width = 8,
  height = 7,
  dpi = 300
)

message("Saved: ", out_png)
message("Done ✅")