#!/usr/bin/env Rscript
# Harmonic mean Ne by full IUCN status (CR/EN/VU/NT/LC), one plot per
# time period (Overall + the 4 fixed calendar bins) per method (PSMC,
# MSMC2) -- lets you see whether historical Ne differs by threat status,
# and whether that pattern is stable across time periods or only shows up
# in one. Same visual grammar as froh_by_iucn_status.png.
#
# Reads the small CSVs that analyze_ne_harmonic_means.R already writes
# (ne_harmonic_means_msmc2.csv / ne_harmonic_means_psmc.csv), so this is
# fast -- no need to re-parse every .psmc/msmc2.final.txt file again.
#
# Usage:
#   Rscript plot_ne_harmonic_by_iucn.R [ne_summary_dir] [iucn_csv] [out_dir]
#   Rscript plot_ne_harmonic_by_iucn.R plots/ne_summary /global/scratch/users/julesperez/snakemake_pipeline/metadata/complete_2026.csv plots/ne_summary

args <- commandArgs(trailingOnly = TRUE)
ne_summary_dir <- if (length(args) >= 1 && nzchar(args[1])) args[1] else "plots/ne_summary"
iucn_csv <- if (length(args) >= 2 && nzchar(args[2])) args[2] else stop("iucn_csv is required (scientificName,redlistCategory columns)")
out_dir  <- if (length(args) >= 3 && nzchar(args[3])) args[3] else ne_summary_dir

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

IUCN_ORDER  <- c("CR", "EN", "VU", "NT", "LC")
IUCN_COLORS <- c(CR = "#B72E3C", EN = "#D9772E", VU = "#C9A227",
                  NT = "#7A9D54", LC = "#2F6B4F")
# 300 DPI instead of R's default 72 -- scale pixel dimensions up
# proportionally to res so physical layout (margins, font size relative to
# plot) stays identical, only pixel density increases.
open_png <- function(path, width, height, res = 300) {
  scale <- res / 72
  png(path, width = width * scale, height = height * scale, res = res)
}

tufte_par <- function(mar = c(2, 5, 4, 2)) {
  par(bty = "n", family = "sans", las = 1, mar = mar,
      tck = -0.015, cex.axis = 0.85, col.axis = "grey30", col.lab = "grey20")
}

normalize_iucn <- function(cat) {
  if (is.na(cat) || !nzchar(cat)) return(NA_character_)
  cl <- tolower(cat)
  if (grepl("critically endangered", cl)) return("CR")
  if (grepl("endangered", cl))            return("EN")
  if (grepl("vulnerable", cl))            return("VU")
  if (grepl("near threatened", cl))       return("NT")
  if (grepl("least concern", cl))         return("LC")
  if (grepl("extinct in the wild", cl))   return("EW")
  if (grepl("extinct", cl))               return("EX")
  if (grepl("data deficient", cl))        return("DD")
  NA_character_
}
iucn_lookup_df <- read.csv(iucn_csv, stringsAsFactors = FALSE)
iucn_for_species <- function(species) {
  sci <- gsub("_", " ", species)
  row <- iucn_lookup_df[tolower(iucn_lookup_df$scientificName) == tolower(sci), ]
  if (nrow(row) == 0) return(NA_character_)
  normalize_iucn(row$redlistCategory[1])
}

plot_by_iucn <- function(values, iucn, out_path, x_label) {
  keep <- !is.na(iucn) & iucn %in% IUCN_ORDER & is.finite(values) & values > 0
  v_all <- values[keep]; iucn_all <- iucn[keep]
  if (length(v_all) == 0) {
    cat("Skipping", out_path, "-- no samples with both a value and an IUCN match\n")
    return(invisible())
  }
  present_order <- IUCN_ORDER[IUCN_ORDER %in% iucn_all]
  n <- length(present_order)
  set.seed(1)
  open_png(out_path, width = 950, height = 550)
  tufte_par()
  xr <- range(v_all)
  plot(NA, xlim = xr, ylim = c(0.3, n + 0.7), log = "x", xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = "")
  axis(3, at = 10^pretty(log10(xr), n = 5), lwd = 0.6, cex.axis = 0.8)
  mtext(x_label, side = 3, line = 2.2, cex = 0.95, col = "grey20")
  for (i in seq_len(n)) {
    cat_code <- present_order[i]
    y <- n - i + 1
    v <- v_all[iucn_all == cat_code]
    col <- IUCN_COLORS[[cat_code]]
    if (length(v) >= 2) {
      qs <- quantile(v, c(0.25, 0.5, 0.75), na.rm = TRUE)
      rect(qs[1], y - 0.22, qs[3], y + 0.22, col = adjustcolor(col, alpha.f = 0.15), border = col, lwd = 1)
      segments(qs[2], y - 0.22, qs[2], y + 0.22, col = col, lwd = 2)
    }
    yj <- y + (stats::runif(length(v)) - 0.5) * 0.32
    points(v, yj, pch = 16, col = adjustcolor(col, alpha.f = 0.75), cex = 1.1)
  }
  axis(2, at = seq_len(n), labels = rev(present_order), lwd = 0, cex.axis = 0.95, font = 2)
  dev.off()
  cat("Wrote", out_path, "(", length(v_all), "samples)\n")
}

run_method <- function(method_name, csv_path) {
  if (!file.exists(csv_path)) {
    cat("Skipping", method_name, "-- no file at", csv_path, "\n")
    return(invisible())
  }
  tbl <- read.csv(csv_path, stringsAsFactors = FALSE)
  tbl$iucn <- sapply(tbl$species, iucn_for_species)
  for (p in unique(tbl$period)) {
    sub <- tbl[tbl$period == p, ]
    period_slug <- gsub("[^A-Za-z0-9]+", "_", tolower(p))
    out_path <- file.path(out_dir, paste0("ne_harmonic_by_iucn_", tolower(method_name), "_", period_slug, ".png"))
    plot_by_iucn(sub$harmonic_mean_ne, sub$iucn, out_path,
                 bquote(.(method_name) ~ N[e] ~ "(harmonic mean)," ~ .(p)))
  }
}

run_method("MSMC2", file.path(ne_summary_dir, "ne_harmonic_means_msmc2.csv"))
run_method("PSMC", file.path(ne_summary_dir, "ne_harmonic_means_psmc.csv"))
