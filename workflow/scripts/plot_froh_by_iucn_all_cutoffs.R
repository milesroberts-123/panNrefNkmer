#!/usr/bin/env Rscript
# FROH by full IUCN status (CR/EN/VU/NT/LC, not just LC-vs-everything-else),
# one plot per ROH-segment cutoff definition -- lets you see whether the
# threat-status pattern holds up consistently across cutoffs, or is an
# artifact of any one particular definition. Same visual grammar as
# froh_by_iucn_status.png: top-mounted axis, IQR box + median tick,
# jittered raw points, categories colored and labeled directly (no legend).
#
# Reads the small froh_comparison.csv that plot_msmc_roh.R already writes,
# so this is fast -- no need to re-run the expensive per-sample scan.
#
# Usage:
#   Rscript plot_froh_by_iucn_all_cutoffs.R [froh_csv] [iucn_csv] [out_dir]
#   Rscript plot_froh_by_iucn_all_cutoffs.R plots/roh/froh_comparison.csv /global/scratch/users/julesperez/snakemake_pipeline/metadata/complete_2026.csv plots/roh

args <- commandArgs(trailingOnly = TRUE)
froh_csv <- if (length(args) >= 1 && nzchar(args[1])) args[1] else "plots/roh/froh_comparison.csv"
iucn_csv <- if (length(args) >= 2 && nzchar(args[2])) args[2] else stop("iucn_csv is required (scientificName,redlistCategory columns)")
out_dir  <- if (length(args) >= 3 && nzchar(args[3])) args[3] else "plots/roh"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
froh_df <- read.csv(froh_csv, stringsAsFactors = FALSE)

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
froh_df$iucn <- sapply(froh_df$species, iucn_for_species)

plot_by_iucn <- function(values, out_path, x_label) {
  keep <- !is.na(froh_df$iucn) & froh_df$iucn %in% IUCN_ORDER & is.finite(values)
  v_all <- values[keep]; iucn_all <- froh_df$iucn[keep]
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
  xr <- c(max(0, xr[1] - 0.02 * diff(xr)), xr[2] + 0.05 * diff(xr))
  plot(NA, xlim = xr, ylim = c(0.3, n + 0.7), xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = "")
  axis(3, at = pretty(xr), lwd = 0.6, cex.axis = 0.85)
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
  cat("Wrote", out_path, "(", length(v_all), "of", nrow(froh_df), "samples had an IUCN match)\n")
}

CUTOFFS <- list(
  list(col = "pct_genome_in_roh",              file = "froh_by_iucn_status_nofilter.png",   label = "% of genome in ROH (FROH), no length cutoff"),
  list(col = "pct_genome_in_roh_1pctchrom",     file = "froh_by_iucn_status_1pctchrom.png",  label = "% of genome in ROH (FROH), ROH >= 1% of chromosome"),
  list(col = "pct_genome_in_roh_1mb",           file = "froh_by_iucn_status_1mb.png",        label = "% of genome in ROH (FROH), ROH >= 1 Mb"),
  list(col = "pct_genome_in_roh_2_5mb",         file = "froh_by_iucn_status_2_5mb.png",      label = "% of genome in ROH (FROH), ROH >= 2.5 Mb"),
  list(col = "pct_genome_in_roh_1pctgenome",    file = "froh_by_iucn_status_1pctgenome.png", label = "% of genome in ROH (FROH), ROH >= 1% of genome")
)
for (cfg in CUTOFFS) {
  if (!(cfg$col %in% names(froh_df))) {
    cat("Skipping", cfg$file, "-- column", cfg$col, "not in", froh_csv, "\n")
    next
  }
  plot_by_iucn(froh_df[[cfg$col]], file.path(out_dir, cfg$file), cfg$label)
}
