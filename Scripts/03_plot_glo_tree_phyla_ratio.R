#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(dplyr)
  library(ggplot2)
  library(ggtree)
  library(ggtreeExtra)
  library(ggnewscale)
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

make_phyla_palette <- function(levels_vec) {
  cols <- scales::hue_pal()(length(levels_vec))
  names(cols) <- levels_vec
  cols
}

# add one heatmap ring
add_ring <- function(p, df, value_col, legend_title, offset = 0.08, pwidth = 0.10) {
  p + ggnewscale::new_scale_fill() +
    geom_fruit(
      data = df,
      geom = geom_tile,
      mapping = aes(y = ID, x = legend_title, fill = .data[[value_col]]),
      offset = offset,
      pwidth = pwidth
    ) +
    scale_fill_gradientn(
      colours = c("#ffffcc", "#c2e699", "#78c679", "#238443"),
      limits  = c(0, 1),
      oob     = scales::squish,
      na.value= "grey85",
      name    = legend_title
    )
}

# ---------------------------
# inputs (ONLY GLO)
# ---------------------------
tree_file <- "Input_Data/bacteria_globdb_r226.tree"
glo_rds   <- "data/processed/glo_input.rds"

out_png    <- "results/figures/Fig4_tree_GLO_phyla_ratio_ranked.png"
out_pruned <- "results/trees/pruned_tree_GLO.tree"
out_layers_dir <- "results/layers"  # optional exports

stop_if_missing(tree_file, "tree file")
stop_if_missing(glo_rds,   "glo_input.rds")

make_dir(dirname(out_png))
make_dir(dirname(out_pruned))
make_dir(out_layers_dir)

# ---------------------------
# 1) read tree
# ---------------------------
tree <- ape::read.tree(tree_file)

# ---------------------------
# 2) read glo table
# ---------------------------
glo <- readRDS(glo_rds)

# pick an ID column that exists
id_col <- if ("id" %in% names(glo)) "id" else if ("taxon_oid" %in% names(glo)) "taxon_oid" else if ("Bin.ID" %in% names(glo)) "Bin.ID" else NA_character_
if (is.na(id_col)) stop("No ID column found in glo_input.rds (need id/taxon_oid/Bin.ID).", call. = FALSE)

# pick a phylum column that exists
phy_col <- if ("Phylum" %in% names(glo)) "Phylum" else if ("phylum" %in% names(glo)) "phylum" else if ("Phyla" %in% names(glo)) "Phyla" else NA_character_
if (is.na(phy_col)) stop("No phylum column found in glo_input.rds (need Phylum/phylum/Phyla).", call. = FALSE)

# ratio column must exist
if (!("ABC_GH_ratio" %in% names(glo))) stop("Column 'ABC_GH_ratio' missing in glo_input.rds.", call. = FALSE)

glo <- glo %>%
  mutate(
    ID = as.character(.data[[id_col]]),
    Phyla = trimws(as.character(.data[[phy_col]])),
    ABC_GH_ratio = suppressWarnings(as.numeric(ABC_GH_ratio))
  ) %>%
  filter(
    !is.na(ID), ID != "",
    !is.na(Phyla), Phyla != "",
    is.finite(ABC_GH_ratio)
  ) %>%
  distinct(ID, .keep_all = TRUE)

if (nrow(glo) < 3) stop("Too few rows in GLO after cleaning; check ID/Phyla/ratio columns.", call. = FALSE)

# ---------------------------
# 3) build dt1 (ID + Phyla, top15 + Other)
# ---------------------------
dt1 <- glo %>% transmute(ID, Phyla)

top15 <- dt1 %>% count(Phyla, sort = TRUE) %>% slice_head(n = 15) %>% pull(Phyla)

dt1 <- dt1 %>%
  mutate(Phyla = ifelse(Phyla %in% top15, Phyla, "Other"))

phyla_levels <- dt1 %>% count(Phyla, sort = TRUE) %>% pull(Phyla)
dt1 <- dt1 %>% mutate(Phyla = factor(Phyla, levels = phyla_levels))

phyla_cols <- make_phyla_palette(levels(dt1$Phyla))

# ---------------------------
# 4) build dt2 (ID + Ratio_rank, ranked WITHIN GLO)
# ---------------------------
dt2 <- glo %>%
  transmute(
    ID = ID,
    Ratio_rank = dplyr::percent_rank(ABC_GH_ratio)
  ) %>%
  mutate(Ratio_rank = pmin(pmax(Ratio_rank, 0), 1))

# optional exports (nice for GitHub reproducibility)
write.csv(dt1 %>% mutate(Phyla = as.character(Phyla)),
          file = file.path(out_layers_dir, "layer_phyla_GLO.csv"),
          row.names = FALSE)
write.csv(dt2,
          file = file.path(out_layers_dir, "layer_ratio_GLO_rank.csv"),
          row.names = FALSE)

# ---------------------------
# 5) prune tree to IDs in dt1
# ---------------------------
ids_to_keep <- intersect(tree$tip.label, dt1$ID)
if (length(ids_to_keep) < 3) {
  stop("Too few overlapping IDs between tree tips and GLO IDs. Likely ID naming mismatch.", call. = FALSE)
}

pruned_tree <- ape::keep.tip(tree, ids_to_keep)
ape::write.tree(pruned_tree, file = out_pruned)

# keep only those on pruned tree
tip_ids <- pruned_tree$tip.label
dt1 <- dt1 %>% filter(ID %in% tip_ids)
dt2 <- dt2 %>% filter(ID %in% tip_ids)

message("Pruned tips: ", length(tip_ids))
message("Saved pruned tree: ", out_pruned)

# ---------------------------
# 6) plot
# ---------------------------
p_tree <- ggtree(pruned_tree, layout = "circular", linewidth = 0.2)

p_tip <- p_tree %<+% dt1 +
  geom_tippoint(
    aes(fill = Phyla),
    shape = 21, size = 0.7, stroke = 0, alpha = 0.9
  ) +
  scale_fill_manual(values = phyla_cols, name = "Phyla")

# ONLY ONE RING: GLO
p_final <- add_ring(
  p_tip,
  dt2,
  value_col = "Ratio_rank",
  legend_title = "GLO cheating index\n(ranked)",
  offset = 0.08,
  pwidth = 0.12
) +
  theme(
    legend.key.height = unit(0.5, "cm"),
    legend.key.width  = unit(0.4, "cm")
  )

ggsave(out_png, p_final, width = 9, height = 9, dpi = 300)
message("Saved: ", out_png)
message("Done ✅")