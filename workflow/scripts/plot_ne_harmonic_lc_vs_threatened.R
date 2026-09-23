#!/usr/bin/env Rscript
# Harmonic mean Ne, LC (non-threatened) vs everything else (threatened),
# one plot per time period (Overall + the 4 fixed calendar bins) per method
# (PSMC, MSMC2). Binary grouping instead of the full CR/EN/VU/NT/LC split in
# plot_ne_harmonic_by_iucn.R -- much better balanced sample sizes per group,
# at the cost of losing resolution within the threatened categories. Same
# rangebox visual style and colors as plot_froh_lc_vs_threatened.R.
#
# Reads the small CSVs that analyze_ne_harmonic_means.R already writes, so
# this is fast -- no need to re-parse every .psmc/msmc2.final.txt file again.
#
# Usage:
#   Rscript plot_ne_harmonic_lc_vs_threatened.R [ne_summary_dir] [iucn_csv] [out_dir]
#   Rscript plot_ne_harmonic_lc_vs_threatened.R plots/ne_summary /global/scratch/users/julesperez/snakemake_pipeline/metadata/complete_2026.csv plots/ne_summary

args <- commandArgs(trailingOnly = TRUE)
ne_summary_dir <- if (length(args) >= 1 && nzchar(args[1])) args[1] else "plots/ne_summary"
iucn_csv <- if (length(args) >= 2 && nzchar(args[2])) args[2] else stop("iucn_csv is required (scientificName,redlistCategory columns)")
out_dir  <- if (length(args) >= 3 && nzchar(args[3])) args[3] else ne_summary_dir

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

open_png <- function(path, width, height, res = 300) {
  scale <- res / 72
  png(path, width = width * scale, height = height * scale, res = res)
}
tufte_par <- function(mar = c(3, 4, 3, 2)) {
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

CUTOFF_ORANGE <- "#E8A33D"
CUTOFF_GREEN  <- "#2D5A34"
rangebox_row <- function(y_center, values, color, height = 0.32) {
  values <- values[is.finite(values)]
  if (length(values) == 0) return(invisible())
  rng <- range(values)
  segments(rng[1], y_center, rng[2], y_center, col = color, lwd = 1)
  if (length(values) >= 2) {
    qs <- quantile(values, c(0.25, 0.5, 0.75))
    rect(qs[1], y_center - height, qs[3], y_center + height,
         col = adjustcolor(color, alpha.f = 0.18), border = color, lwd = 1.2)
    segments(qs[2], y_center - height, qs[2], y_center + height, col = color, lwd = 2.2)
  }
  yj <- y_center + (stats::runif(length(values)) - 0.5) * height * 1.7
  points(values, yj, pch = 21, bg = adjustcolor(color, alpha.f = 0.55),
         col = adjustcolor("black", alpha.f = 0.5), cex = 1.15, lwd = 0.6)
}

plot_lc_vs_threatened <- function(values, group, out_path, x_label, log_scale = TRUE) {
  keep <- !is.na(group) & is.finite(values) & values > 0
  v <- values[keep]; g <- group[keep]
  if (length(v) == 0 || length(unique(g)) < 2) {
    cat("Skipping", out_path, "-- not enough samples in both groups\n")
    return(invisible())
  }
  xr <- range(v, na.rm = TRUE)
  set.seed(1)
  open_png(out_path, width = 950, height = 560)
  # mar left = 13 (not 11) -- "Threatened (NT/VU/EN/CR/EW)" was clipping
  # against the plot edge at the 300 DPI scale.
  tufte_par(mar = c(5, 13, 2, 2))
  plot(NA, xlim = xr, ylim = c(0.4, 2.6), log = if (log_scale) "x" else "",
       xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = "")
  rangebox_row(2, v[g == "Threatened (NT/VU/EN/CR/EW)"], CUTOFF_ORANGE)
  rangebox_row(1, v[g == "Non-threatened (LC)"], CUTOFF_GREEN)
  axis(2, at = c(2, 1), labels = c("Threatened (NT/VU/EN/CR/EW)", "Non-threatened (LC)"),
       lwd = 0, cex.axis = 1, las = 1)
  if (log_scale) {
    axis(1, at = 10^pretty(log10(xr), n = 5), lwd = 0.6, cex.axis = 0.85)
  } else {
    axis(1, at = pretty(xr, n = 5), lwd = 0.6, cex.axis = 0.85)
  }
  mtext(x_label, side = 1, line = 3, font = 2, cex = 1, col = "grey20")
  dev.off()
  cat("Wrote", out_path, "\n")
}

run_method <- function(method_name, csv_path) {
  if (!file.exists(csv_path)) {
    cat("Skipping", method_name, "-- no file at", csv_path, "\n")
    return(invisible())
  }
  tbl <- read.csv(csv_path, stringsAsFactors = FALSE)
  tbl$iucn <- sapply(tbl$species, iucn_for_species)
  tbl$threat_group <- ifelse(tbl$iucn == "LC", "Non-threatened (LC)",
                              ifelse(is.na(tbl$iucn), NA, "Threatened (NT/VU/EN/CR/EW)"))
  for (p in unique(tbl$period)) {
    sub <- tbl[tbl$period == p, ]
    period_slug <- gsub("[^A-Za-z0-9]+", "_", tolower(p))
    base_name <- paste0("ne_harmonic_lc_vs_threatened_", tolower(method_name), "_", period_slug)
    label <- paste0(method_name, " Ne (harmonic mean), ", p)
    plot_lc_vs_threatened(sub$harmonic_mean_ne, sub$threat_group,
                           file.path(out_dir, paste0(base_name, ".png")), label, log_scale = TRUE)
    plot_lc_vs_threatened(sub$harmonic_mean_ne, sub$threat_group,
                           file.path(out_dir, paste0(base_name, "_linear.png")), label, log_scale = FALSE)
  }
}

run_method("MSMC2", file.path(ne_summary_dir, "ne_harmonic_means_msmc2.csv"))
run_method("PSMC", file.path(ne_summary_dir, "ne_harmonic_means_psmc.csv"))
