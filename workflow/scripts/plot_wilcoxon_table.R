#!/usr/bin/env Rscript
# Meeting-ready table figure of the Wilcoxon LC-vs-threatened results --
# reads the CSV test_lc_vs_threatened_significance.R already wrote, sorts
# by raw p-value (most promising first), and highlights rows with raw
# p<0.05 -- while still showing the BH-adjusted column plainly next to it,
# since none of those currently survive correction. The highlight marks
# "worth discussing," not "proven."
#
# Usage:
#   Rscript plot_wilcoxon_table.R [significance_csv] [out_path]
#   Rscript plot_wilcoxon_table.R plots/lc_vs_threatened_significance.csv plots/lc_vs_threatened_significance_table.png

args <- commandArgs(trailingOnly = TRUE)
sig_csv <- if (length(args) >= 1 && nzchar(args[1])) args[1] else "plots/lc_vs_threatened_significance.csv"
out_path <- if (length(args) >= 2 && nzchar(args[2])) args[2] else "plots/lc_vs_threatened_significance_table.png"

open_png <- function(path, width, height, res = 300) {
  scale <- res / 72
  png(path, width = width * scale, height = height * scale, res = res)
}

tbl <- read.csv(sig_csv, stringsAsFactors = FALSE)
tbl <- tbl[order(tbl$p_value), ]

fmt_num <- function(x, digits = 3) ifelse(is.na(x), "--", formatC(x, format = "g", digits = digits))
fmt_p <- function(x) ifelse(is.na(x), "--", ifelse(x < 0.001, "<0.001", formatC(x, format = "f", digits = 3)))

rows <- data.frame(
  Comparison = tbl$comparison,
  n_nt = tbl$n_non_threatened,
  n_th = tbl$n_threatened,
  med_nt = fmt_num(tbl$median_non_threatened),
  med_th = fmt_num(tbl$median_threatened),
  p = fmt_p(tbl$p_value),
  p_adj = fmt_p(tbl$p_adj_BH),
  sig_raw = ifelse(!is.na(tbl$p_value) & tbl$p_value < 0.05, TRUE, FALSE),
  sig_bh = ifelse(!is.na(tbl$significant_BH_0.05) & tbl$significant_BH_0.05, TRUE, FALSE),
  stringsAsFactors = FALSE
)

HEADER <- c("Comparison", "n (non-threat.)", "n (threat.)", "Median (non-threat.)",
            "Median (threat.)", "p-value", "p-adj (BH)")
col_x <- c(20, 700, 850, 1010, 1170, 1330, 1450)  # pixel x positions, pre-DPI-scale
col_left_align <- c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)

HIGHLIGHT_BG <- "#FCEFD8"
HIGHLIGHT_BORDER <- "#E8A33D"
BH_BG <- "#E4EEDD"
BH_BORDER <- "#2D5A34"
TEXT_DARK <- "#2A2A2A"
TEXT_MUTED <- "#6B6B6B"

# Fixed pixel row heights (pre-DPI-scale, matching open_png's own scaling)
# -- avoids the row spacing distorting as n_rows changes, unlike a layout
# based on fractions of total plot height.
TITLE_H <- 55
HEADER_H <- 40
ROW_H <- 46
LEGEND_H <- 40
TOP_PAD <- 25
BOTTOM_PAD <- 20
GAP_AFTER_TITLE <- 15
GAP_AFTER_HEADER <- 10
GAP_BEFORE_LEGEND <- 20

n_rows <- nrow(rows)
WIDTH <- 1520
HEIGHT <- TOP_PAD + TITLE_H + GAP_AFTER_TITLE + HEADER_H + GAP_AFTER_HEADER +
  n_rows * ROW_H + GAP_BEFORE_LEGEND + LEGEND_H + BOTTOM_PAD

open_png(out_path, width = WIDTH, height = HEIGHT)
par(mar = c(0, 0, 0, 0), family = "sans")
plot(NA, xlim = c(0, WIDTH), ylim = c(0, HEIGHT), xaxs = "i", yaxs = "i",
     xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")

# y coordinates count DOWN from the top (HEIGHT), matching reading order
y <- HEIGHT - TOP_PAD
title_y <- y - TITLE_H * 0.5
text(col_x[1], title_y, "LC (non-threatened) vs Threatened -- Wilcoxon rank-sum tests",
     adj = c(0, 0.5), cex = 1.35, font = 2, col = TEXT_DARK)
y <- y - TITLE_H - GAP_AFTER_TITLE

header_y <- y - HEADER_H * 0.5
for (i in seq_along(HEADER)) {
  text(col_x[i], header_y, HEADER[i],
       adj = c(if (col_left_align[i]) 0 else 1, 0.5), cex = 0.88, font = 2, col = TEXT_MUTED)
}
y <- y - HEADER_H
segments(20, y, WIDTH - 20, y, col = "grey70", lwd = 1)
y <- y - GAP_AFTER_HEADER

for (i in seq_len(n_rows)) {
  row_top <- y - (i - 1) * ROW_H
  row_bottom <- row_top - ROW_H
  row_center <- (row_top + row_bottom) / 2
  r <- rows[i, ]
  if (r$sig_bh) {
    rect(10, row_bottom + 3, WIDTH - 10, row_top - 3, col = BH_BG, border = BH_BORDER, lwd = 1.2)
  } else if (r$sig_raw) {
    rect(10, row_bottom + 3, WIDTH - 10, row_top - 3, col = HIGHLIGHT_BG, border = HIGHLIGHT_BORDER, lwd = 1.2)
  }
  font_w <- if (r$sig_raw) 2 else 1
  text(col_x[1], row_center, r$Comparison, adj = c(0, 0.5), cex = 0.82, font = font_w, col = TEXT_DARK)
  text(col_x[2], row_center, r$n_nt, adj = c(1, 0.5), cex = 0.82, col = TEXT_DARK)
  text(col_x[3], row_center, r$n_th, adj = c(1, 0.5), cex = 0.82, col = TEXT_DARK)
  text(col_x[4], row_center, r$med_nt, adj = c(1, 0.5), cex = 0.82, col = TEXT_DARK)
  text(col_x[5], row_center, r$med_th, adj = c(1, 0.5), cex = 0.82, col = TEXT_DARK)
  text(col_x[6], row_center, r$p, adj = c(1, 0.5), cex = 0.82, font = font_w, col = TEXT_DARK)
  text(col_x[7], row_center, r$p_adj, adj = c(1, 0.5), cex = 0.82, col = TEXT_DARK)
}
y <- y - n_rows * ROW_H - GAP_BEFORE_LEGEND

legend_y <- y - LEGEND_H * 0.5
rect(col_x[1], legend_y - 8, col_x[1] + 24, legend_y + 8, col = HIGHLIGHT_BG, border = HIGHLIGHT_BORDER)
text(col_x[1] + 32, legend_y, "raw p < 0.05 (not yet significant after multiple-testing correction)",
     adj = c(0, 0.5), cex = 0.72, col = TEXT_MUTED)

dev.off()
cat("Wrote", out_path, "\n")
