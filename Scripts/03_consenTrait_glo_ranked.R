#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ape)
  library(phytools)
  library(castor)
})

set.seed(1)

############################################################
## params
############################################################
MIN_FRAC         <- 0.9
COUNT_SINGLETONS <- TRUE
WEIGHTED_CLADES  <- FALSE
N_PERM           <- 1000
N_BOOT           <- 1000

Q_TOP <- 0.90   # top 10% => cheater=1

############################################################
## inputs (edit)
############################################################
tree_in  <- "Input_Data/bacteria_globdb_r226.tree"

img_rds  <- "data/processed/img_input.rds"
loma_rds <- "data/processed/loma_input.rds"
iso_rds  <- "data/processed/isolates_input.rds"
glo_rds  <- "data/processed/glo_input.rds"

out_tree_rooted <- "results/consentrait/consenTrait_rooted.nwk"
out_trait_bin   <- "results/consentrait/consenTrait_trait_top10pct.txt"

############################################################
## helpers
############################################################
stop_if_missing <- function(path, label=NULL){
  if (!file.exists(path)) stop(paste0(label %||% "File", " not found: ", path), call.=FALSE)
}
`%||%` <- function(a,b) if (!is.null(a)) a else b

make_dir <- function(dir) if (!dir.exists(dir)) dir.create(dir, recursive=TRUE, showWarnings=FALSE)

# unify ID + ratio
unify_cols <- function(df){
  df %>%
    mutate(
      taxon_oid = if ("taxon_oid" %in% names(.)) as.character(.data$taxon_oid) else NA_character_,
      id        = if ("id"        %in% names(.)) as.character(.data$id)        else NA_character_,
      Bin.ID    = if ("Bin.ID"    %in% names(.)) as.character(.data$`Bin.ID`)  else NA_character_,
      ID        = dplyr::coalesce(taxon_oid, id, Bin.ID),
      ABC_GH_ratio = suppressWarnings(as.numeric(.data$ABC_GH_ratio))
    ) %>%
    select(ID, ABC_GH_ratio, database) %>%
    filter(!is.na(ID), ID != "", is.finite(ABC_GH_ratio))
}

# tree height (ape-compatible)
get_tree_height_ape <- function(tr){
  d <- ape::node.depth.edgelength(tr)
  max(d[seq_len(length(tr$tip.label))], na.rm=TRUE)
}

# consenTRAIT wrapper: unit-height + bootstrap CI
run_consenTrait <- function(tree_in, tip_states,
                            min_fraction=MIN_FRAC,
                            count_singletons=COUNT_SINGLETONS,
                            weighted=WEIGHTED_CLADES,
                            Npermutations=N_PERM,
                            Nboot=N_BOOT){
  H <- get_tree_height_ape(tree_in)
  stopifnot(is.finite(H), H > 0)
  
  tree_unit <- tree_in
  tree_unit$edge.length <- tree_unit$edge.length / H
  
  has_sr <- "singleton_resolution" %in% names(formals(castor::consentrait_depth))
  args <- list(tree=tree_unit, tip_states=tip_states,
               min_fraction=min_fraction,
               count_singletons=count_singletons,
               weighted=weighted,
               Npermutations=Npermutations)
  if (has_sr) args$singleton_resolution <- 0
  res <- do.call(castor::consentrait_depth, args)
  
  tau_rel <- res$mean_depth
  tau_abs <- tau_rel * H
  
  # stratified bootstrap (pos/neg)
  pos <- which(tip_states==1L); neg <- which(tip_states==0L)
  boot_tau <- replicate(Nboot, {
    yb <- integer(length(tip_states))
    yb[sample(pos, length(pos), replace=TRUE)] <- 1L
    yb[sample(neg, length(neg), replace=TRUE)] <- 0L
    args2 <- list(tree=tree_unit, tip_states=yb,
                  min_fraction=min_fraction,
                  count_singletons=count_singletons,
                  weighted=weighted,
                  Npermutations=0)
    if (has_sr) args2$singleton_resolution <- 0
    do.call(castor::consentrait_depth, args2)$mean_depth
  })
  CI_rel <- stats::quantile(boot_tau, c(0.025, 0.975), na.rm=TRUE)
  
  list(H=H, tau_abs=tau_abs, tau_rel=tau_rel, CI_rel=CI_rel,
       P=res$P, Npositives=res$Npositives, raw=res)
}

############################################################
## 0) checks + dirs
############################################################
stop_if_missing(tree_in, "tree")
stop_if_missing(img_rds, "img_rds")
stop_if_missing(loma_rds,"loma_rds")
stop_if_missing(iso_rds, "iso_rds")
stop_if_missing(glo_rds, "glo_rds")

make_dir(dirname(out_tree_rooted))
make_dir(dirname(out_trait_bin))

############################################################
## 1) load 4 DBs and compute GLOBAL rank across all
############################################################
img  <- readRDS(img_rds)  %>% mutate(database="IMG")       %>% unify_cols()
loma <- readRDS(loma_rds) %>% mutate(database="LOMA")      %>% unify_cols()
iso  <- readRDS(iso_rds)  %>% mutate(database="isolates")  %>% unify_cols()
glo  <- readRDS(glo_rds)  %>% mutate(database="glo")       %>% unify_cols()

all_df <- bind_rows(loma, img, iso, glo) %>%
  mutate(Ratio_rank = dplyr::percent_rank(ABC_GH_ratio))

# threshold on GLOBAL ranks
q_top <- stats::quantile(all_df$Ratio_rank, Q_TOP, na.rm=TRUE)

all_df <- all_df %>%
  mutate(cheater = if_else(Ratio_rank >= q_top, 1L, 0L))

cat("Global threshold q", Q_TOP, "=", q_top, "\n")
cat("Global positives:", sum(all_df$cheater==1L), " / ", nrow(all_df), "\n")

############################################################
## 2) restrict trait to GLO tips on this tree, align, root tree
############################################################
tree <- ape::read.tree(tree_in)

# only use IDs that are in GLO and in tree
glo_trait <- all_df %>%
  filter(database=="glo") %>%
  distinct(ID, .keep_all=TRUE) %>%
  select(ID, cheater)

shared <- intersect(tree$tip.label, glo_trait$ID)
cat("Tree tips:", length(tree$tip.label),
    "| GLO IDs:", nrow(glo_trait),
    "| shared:", length(shared), "\n")

if (length(shared) < 10) {
  stop("Too few shared IDs between tree tips and GLO trait IDs. Check ID naming.", call.=FALSE)
}

tree_use <- drop.tip(tree, setdiff(tree$tip.label, shared))

# midpoint root AFTER pruning
tree_rooted <- phytools::midpoint.root(tree_use)
write.tree(tree_rooted, out_tree_rooted)

# align trait to rooted tree order
trait_vec <- setNames(glo_trait$cheater, glo_trait$ID)
tip_states <- trait_vec[tree_rooted$tip.label]
if (any(is.na(tip_states))) {
  # should not happen because we pruned by shared, but keep safe
  keep <- !is.na(tip_states)
  tree_rooted <- drop.tip(tree_rooted, tree_rooted$tip.label[!keep])
  tip_states <- tip_states[keep]
}

# export trait file (2 cols, no header)
write.table(
  data.frame(id = names(tip_states), trait = as.integer(tip_states)),
  file = out_trait_bin,
  sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
)

cat("Saved rooted tree:", out_tree_rooted, "\n")
cat("Saved trait file :", out_trait_bin, "\n")
cat("GLO tips used:", length(tip_states), "| positives:", sum(tip_states==1L), "\n")

############################################################
## 3) run castor consenTRAIT (main + sensitivity)
############################################################
res_main <- run_consenTrait(tree_rooted, tip_states, count_singletons=TRUE)
cat("\n[Main] Relative tauD:", res_main$tau_rel,
    " 95% CI:", paste(round(res_main$CI_rel,4), collapse=" –"),
    " | Absolute tauD:", res_main$tau_abs,
    " | P:", res_main$P,
    " | Npos:", res_main$Npositives, "\n")

res_sens <- run_consenTrait(tree_rooted, tip_states, count_singletons=FALSE)
cat("[Sensitivity] Relative tauD (no singletons):", res_sens$tau_rel,
    " 95% CI:", paste(round(res_sens$CI_rel,4), collapse=" –"),
    " | P:", res_sens$P,
    " | Npos:", res_sens$Npositives, "\n\n")

cat("Done ✅\n")