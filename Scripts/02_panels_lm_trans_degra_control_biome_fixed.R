#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(broom)
  library(patchwork)
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
# 单变量：拟合 + 预测区间
# ---------------------------
one_var_fit <- function(df, xvar, yvar, n = 200) {
  d <- df %>% filter(is.finite(.data[[xvar]]), is.finite(.data[[yvar]]))
  if (nrow(d) < 3) return(NULL)
  
  fit <- lm(reformulate(xvar, yvar), data = d)
  
  xseq <- seq(min(d[[xvar]], na.rm = TRUE),
              max(d[[xvar]], na.rm = TRUE),
              length.out = n)
  nd <- setNames(data.frame(xseq), xvar)
  
  pr <- predict(fit, newdata = nd, interval = "confidence")
  
  tibble(
    x = xseq,
    fit = pr[, "fit"],
    lwr = pr[, "lwr"],
    upr = pr[, "upr"]
  )
}

# ---------------------------
# 一个panel里放两个独立x轴（左半边一个，右半边一个）
# 中间没有线，但保留一个小空隙
# ---------------------------
panel_two_half_axes <- function(df, yvar,
                                xvar_left, xvar_right,
                                col_left, col_right,
                                xlab_left, xlab_right,
                                ylab = "trans/degra (raw)",
                                y_expand = 0.08,
                                gap = 0.08) {
  
  res_left  <- one_var_fit(df, xvar_left,  yvar, n = 200)
  res_right <- one_var_fit(df, xvar_right, yvar, n = 200)
  
  if (is.null(res_left))  stop("Not enough data for: ", xvar_left)
  if (is.null(res_right)) stop("Not enough data for: ", xvar_right)
  
  # y 轴按两条线 + CI 的共同范围来定
  y_all <- c(res_left$lwr, res_left$upr, res_right$lwr, res_right$upr)
  y_rng <- range(y_all, na.rm = TRUE)
  
  if (diff(y_rng) == 0) {
    y_pad <- max(abs(y_rng[1]) * y_expand, 0.05)
  } else {
    y_pad <- diff(y_rng) * y_expand
  }
  ylim_use <- c(y_rng[1] - y_pad, y_rng[2] + y_pad)
  
  # 各自原始 x 范围
  xL_rng <- range(res_left$x, na.rm = TRUE)
  xR_rng <- range(res_right$x, na.rm = TRUE)
  
  # 定义左右两块在同一panel中的位置
  # 左块: [0, left_end]
  # 右块: [right_start, 2]
  left_end <- 1 - gap / 2
  right_start <- 1 + gap / 2
  
  # 把左变量映射到 [0, left_end]
  map_left <- function(x) {
    (x - xL_rng[1]) / diff(xL_rng) * left_end
  }
  
  # 把右变量映射到 [right_start, 2]
  map_right <- function(x) {
    right_start + (x - xR_rng[1]) / diff(xR_rng) * (2 - right_start)
  }
  
  left_df <- res_left %>%
    mutate(x_plot = map_left(x))
  
  right_df <- res_right %>%
    mutate(x_plot = map_right(x))
  
  # 自定义底部刻度：左半边和右半边各自一套
  pretty_left_raw  <- pretty(xL_rng, n = 4)
  pretty_left_raw  <- pretty_left_raw[pretty_left_raw >= xL_rng[1] & pretty_left_raw <= xL_rng[2]]
  pretty_left_pos  <- map_left(pretty_left_raw)
  
  pretty_right_raw <- pretty(xR_rng, n = 4)
  pretty_right_raw <- pretty_right_raw[pretty_right_raw >= xR_rng[1] & pretty_right_raw <= xR_rng[2]]
  pretty_right_pos <- map_right(pretty_right_raw)
  
  breaks_all <- c(pretty_left_pos, pretty_right_pos)
  labels_all <- c(pretty_left_raw, pretty_right_raw)
  
  ggplot() +
    # 左半边
    geom_ribbon(
      data = left_df,
      aes(x = x_plot, ymin = lwr, ymax = upr),
      fill = col_left, alpha = 0.18, color = NA
    ) +
    geom_line(
      data = left_df,
      aes(x = x_plot, y = fit),
      color = col_left, linewidth = 1.2
    ) +
    
    # 右半边
    geom_ribbon(
      data = right_df,
      aes(x = x_plot, ymin = lwr, ymax = upr),
      fill = col_right, alpha = 0.18, color = NA
    ) +
    geom_line(
      data = right_df,
      aes(x = x_plot, y = fit),
      color = col_right, linewidth = 1.2
    ) +
    
    coord_cartesian(
      xlim = c(0, 2),
      ylim = ylim_use,
      clip = "off"
    ) +
    scale_x_continuous(
      breaks = breaks_all,
      labels = labels_all,
      expand = c(0.01, 0.01)
    ) +
    labs(
      x = NULL,
      y = ylab
    ) +
    annotate(
      "text",
      x = left_end / 2,
      y = ylim_use[1] - 0.12 * diff(ylim_use),
      label = xlab_left,
      vjust = 1,
      size = 4.5
    ) +
    annotate(
      "text",
      x = (right_start + 2) / 2,
      y = ylim_use[1] - 0.12 * diff(ylim_use),
      label = xlab_right,
      vjust = 1,
      size = 4.5
    ) +
    theme_classic(base_size = 13) +
    theme(
      axis.text = element_text(color = "black"),
      axis.title.y = element_text(color = "black"),
      plot.margin = margin(8, 8, 20, 8)
    )
}

# ---------------------------
# Inputs
# ---------------------------
in_rds <- "data/processed/img_input.rds"
climate_csv <- "Input_Data/img_30s_climate.csv"
soil_csv <- "Input_Data/soil_pro.csv"

stop_if_missing(in_rds, "img_input.rds")
stop_if_missing(climate_csv, "img_30s_climate.csv")
stop_if_missing(soil_csv, "soil_pro.csv")

merged_table <- readRDS(in_rds)
climate <- read.csv(climate_csv, stringsAsFactors = FALSE, check.names = FALSE)
soil <- read.csv(soil_csv, stringsAsFactors = FALSE, check.names = FALSE)

# merge by Bin.ID
meta <- merge(climate, soil, by = "Bin.ID")
meta <- merge(meta, merged_table, by = "Bin.ID")

# longitude/latitude col names compat
lon <- if ("Longitude.x" %in% names(meta)) "Longitude.x" else "Longitude"
lat <- if ("Latitude.x" %in% names(meta))  "Latitude.x"  else "Latitude"

# site-level medians
avg_table <- meta %>%
  group_by(.data[[lon]], .data[[lat]]) %>%
  summarise(
    across(where(is.numeric), ~ median(.x, na.rm = TRUE)),
    n_bin = n_distinct(Bin.ID),
    .groups = "drop"
  )

# ---------------------------
# variables
# ---------------------------
yvar <- "ABC_GH_ratio"

x_clim1 <- "TS"
x_clim2 <- "MAP"

x_soil1 <- "nitrogen_15cm"
x_soil2 <- "clay_15cm"

var_labels <- c(
  TS = "Temperature seasonality",
  MAP = "Annual precipitation (mm)",
  nitrogen_15cm = "Total nitrogen (g/kg)",
  clay_15cm = "Clay content (%)"
)

cols <- c(
  TS = "#2C7FB8",
  MAP = "#41AB5D",
  nitrogen_15cm = "#6A4C93",
  clay_15cm = "#F28E2B"
)

# ---------------------------
# Build 2 panels
# ---------------------------
p_climate <- panel_two_half_axes(
  df = avg_table,
  yvar = yvar,
  xvar_left  = x_clim1,
  xvar_right = x_clim2,
  col_left   = cols[[x_clim1]],
  col_right  = cols[[x_clim2]],
  xlab_left  = var_labels[[x_clim1]],
  xlab_right = var_labels[[x_clim2]],
  ylab = "trans/degra (raw)",
  gap = 0.08
) +
  ggtitle("Climate predictors") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

p_soil <- panel_two_half_axes(
  df = avg_table,
  yvar = yvar,
  xvar_left  = x_soil1,
  xvar_right = x_soil2,
  col_left   = cols[[x_soil1]],
  col_right  = cols[[x_soil2]],
  xlab_left  = var_labels[[x_soil1]],
  xlab_right = var_labels[[x_soil2]],
  ylab = "trans/degra (raw)",
  gap = 0.08
) +
  ggtitle("Soil predictors") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

p_all <- p_climate / p_soil

print(p_all)

# ---------------------------
# Save
# ---------------------------
make_dir("results/figures")
out_pdf <- file.path("results/figures", "Fig3_two_panels_oneframe_two_xaxes_gap.pdf")
out_png <- file.path("results/figures", "Fig3_two_panels_oneframe_two_xaxes_gap.png")

ggsave(out_pdf, p_all, width = 7.5, height = 8.5)
ggsave(out_png, p_all, width = 7.5, height = 8.5, dpi = 300)

message("Saved: ", out_pdf)
message("Saved: ", out_png)
message("Done ✅")