#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tools)
  library(broom)
  library(MASS)      # stepAIC
  library(car)       # vif
  library(scales)    # number()
})

# ---------------------------
# Helpers
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

p2star <- function(p){
  dplyr::case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ ""
  )
}

# ---------------------------
# Forest plot function (FIXED)
# ---------------------------
plot_forest_sig <- function(fit, var_labels = NULL, alpha = 0.05,
                            order = c("asc","desc"),
                            title = "Significant predictors",
                            show_stats = TRUE) {
  order <- match.arg(order)
  
  td <- broom::tidy(fit, conf.int = TRUE) %>%
    dplyr::filter(term != "(Intercept)") %>%
    dplyr::mutate(sig = p2star(p.value)) %>%
    dplyr::filter(p.value < alpha)
  
  if (nrow(td) == 0) {
    message("No terms with p < ", alpha)
    return(
      ggplot() + theme_void() +
        annotate("text", x = 0, y = 0,
                 label = paste0("No terms with p < ", alpha), size = 5)
    )
  }
  
  td <- td %>%
    dplyr::arrange(if (order == "asc") estimate else dplyr::desc(estimate))
  
  plot_dat <- td %>%
    dplyr::mutate(
      # IMPORTANT FIX: use ifelse (vectorized) instead of if (scalar)
      var_label = ifelse(
        !is.null(var_labels) & term %in% names(var_labels),
        unname(var_labels[term]),
        term
      ),
      var_label = factor(var_label, levels = rev(unique(var_label))),
      lab = paste0(
        ifelse(abs(estimate) < 1e-2,
               format(estimate, scientific = TRUE, digits = 2),
               scales::number(estimate, accuracy = 0.01)),
        " ", sig
      ),
      col = ifelse(estimate >= 0, "#67a9cf", "#f05b5b")
    )
  
  ms <- summary(fit)
  stat_text <- sprintf(
    "n = %d | R\u00B2 = %.3f | adj.R\u00B2 = %.3f | AIC = %.1f",
    nobs(fit), ms$r.squared, ms$adj.r.squared, AIC(fit)
  )
  
  p <- ggplot(plot_dat, aes(x = estimate, y = var_label)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey55") +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                   height = 0.18, linewidth = 1, color = plot_dat$col) +
    geom_point(size = 3.2, color = plot_dat$col) +
    geom_text(aes(label = lab), nudge_y = 0.32,
              size = 3.6, color = plot_dat$col) +
    scale_x_continuous("Estimate \u00B1 95% CI",
                       expand = expansion(mult = c(0.06, 0.18))) +
    ylab(NULL) +
    ggtitle(title) +
    theme_classic(base_size = 14) +
    theme(
      axis.text = element_text(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, color = "grey30"),
      plot.margin = margin(10, 24, 10, 10)
    )
  
  if (show_stats) p <- p + labs(subtitle = stat_text)
  p
}

# ---------------------------
# Args
# ---------------------------
args <- commandArgs(trailingOnly = TRUE)

in_rds      <- ifelse(length(args) >= 1, args[1], "data/processed/img_input.rds")
climate_csv <- ifelse(length(args) >= 2, args[2], "Input_Data/climate.csv")
soil_csv    <- ifelse(length(args) >= 3, args[3], "Input_Data/soil_pro.csv")
fig_prefix  <- ifelse(length(args) >= 4, args[4], "Fig_3")

fig_dir   <- "results/figures"
table_dir <- "results/tables"
make_dir(fig_dir)
make_dir(table_dir)

message("=== Stepwise LM + forest plots (clean, fixed) ===")
message("Input RDS:   ", in_rds)
message("Climate CSV: ", climate_csv)
message("Soil CSV:    ", soil_csv)
message("Prefix:      ", fig_prefix)

stop_if_missing(in_rds, "Input RDS")
stop_if_missing(climate_csv, "Climate CSV")
stop_if_missing(soil_csv, "Soil CSV")

# ---------------------------
# Load & merge
# ---------------------------
merged_table <- readRDS(in_rds)
if (!("Bin.ID" %in% colnames(merged_table))) stop("Column 'Bin.ID' not found in merged_table.", call. = FALSE)
if (!("ABC_GH_ratio" %in% colnames(merged_table))) stop("Column 'ABC_GH_ratio' not found in merged_table.", call. = FALSE)

climate <- read.csv(climate_csv, stringsAsFactors = FALSE, check.names = FALSE)
soil    <- read.csv(soil_csv, stringsAsFactors = FALSE, check.names = FALSE)

if (!("Bin.ID" %in% colnames(climate))) stop("Column 'Bin.ID' not found in climate.csv.", call. = FALSE)
if (!("Bin.ID" %in% colnames(soil)))    stop("Column 'Bin.ID' not found in soil_pro.csv.", call. = FALSE)

meta <- merge(climate, soil, by = "Bin.ID")
meta <- merge(meta, merged_table, by = "Bin.ID")

# Longitude / Latitude columns (compat)
lon <- if ("Longitude.x" %in% names(meta)) "Longitude.x" else if ("Longitude" %in% names(meta)) "Longitude" else NA
lat <- if ("Latitude.x"  %in% names(meta)) "Latitude.x"  else if ("Latitude"  %in% names(meta)) "Latitude"  else NA
if (is.na(lon) || is.na(lat)) stop("Longitude/Latitude columns not found.", call. = FALSE)

# Aggregate by site (lon/lat), numeric medians
avg_table <- meta %>%
  group_by(.data[[lon]], .data[[lat]]) %>%
  summarise(
    across(where(is.numeric), ~ median(.x, na.rm = TRUE)),
    n_bin = n_distinct(Bin.ID),
    .groups = "drop"
  )

# ---------------------------
# Variable labels 
# ---------------------------
var_labels <- c(
  MAT = "Annual Mean Temperature",
  MAP = "Annual Precipitation",
  TS  = "Temperature Seasonality",
  PS  = "Precipitation Seasonality",
  nitrogen_15cm = "Total Nitrogen",
  phh2o_15cm    = "Soil pH",
  clay_15cm     = "Clay fraction",
  soc_15cm      = "Soil Organic Carbon"
)

# ---------------------------
# Climate model + stepAIC
# ---------------------------
clim_formula <- ABC_GH_ratio ~ MAT + MAP + TS + PS
miss_clim <- setdiff(all.vars(clim_formula), names(avg_table))
if (length(miss_clim) > 0) stop("Missing variables for climate model: ", paste(miss_clim, collapse = ", "), call. = FALSE)

lm_clim <- lm(clim_formula, data = avg_table)
vif_clim <- car::vif(lm_clim)
step_clim <- MASS::stepAIC(lm_clim, direction = "both", trace = TRUE)

# outputs (climate)
clim_sum_path <- file.path(table_dir, paste0(fig_prefix, "_climate_stepAIC_summary.txt"))
sink(clim_sum_path)
cat("== Initial climate model ==\n"); print(summary(lm_clim))
cat("\n== VIF (initial climate model) ==\n"); print(vif_clim)
cat("\n== StepAIC selected climate model ==\n"); print(summary(step_clim))
sink()

write.csv(
  broom::tidy(step_clim, conf.int = TRUE),
  file.path(table_dir, paste0(fig_prefix, "_climate_stepAIC_coefficients.csv")),
  row.names = FALSE
)

p_clim <- plot_forest_sig(step_clim, var_labels, title = "Significant climatic predictors", show_stats = TRUE)
clim_pdf <- file.path(fig_dir, paste0(fig_prefix, "_climate_forest_sig.pdf"))
ggsave(clim_pdf, p_clim, width = 5.5, height = 3)

message("Saved climate: ", clim_pdf)

# ---------------------------
# Soil model + stepAIC 
# ---------------------------
soil_formula <- ABC_GH_ratio ~ nitrogen_15cm + phh2o_15cm + clay_15cm + soc_15cm
miss_soil <- setdiff(all.vars(soil_formula), names(avg_table))
if (length(miss_soil) > 0) stop("Missing variables for soil model: ", paste(miss_soil, collapse = ", "), call. = FALSE)

lm_soil <- lm(soil_formula, data = avg_table)
vif_soil <- car::vif(lm_soil)
step_soil <- MASS::stepAIC(lm_soil, direction = "both", trace = FALSE)

# outputs (soil)
soil_sum_path <- file.path(table_dir, paste0(fig_prefix, "_soil_stepAIC_summary.txt"))
sink(soil_sum_path)
cat("== Initial soil model ==\n"); print(summary(lm_soil))
cat("\n== VIF (initial soil model) ==\n"); print(vif_soil)
cat("\n== StepAIC selected soil model ==\n"); print(summary(step_soil))
sink()

write.csv(
  broom::tidy(step_soil, conf.int = TRUE),
  file.path(table_dir, paste0(fig_prefix, "_soil_stepAIC_coefficients.csv")),
  row.names = FALSE
)

p_soil <- plot_forest_sig(step_soil, var_labels, title = "Significant edaphic predictors", show_stats = TRUE)
soil_pdf <- file.path(fig_dir, paste0(fig_prefix, "_soil_forest_sig.pdf"))
ggsave(soil_pdf, p_soil, width = 4.5, height = 3)

message("Saved soil: ", soil_pdf)
message("Done ✅")