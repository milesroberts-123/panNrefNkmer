#!/usr/bin/env Rscript
# FROH spread, LC vs everything else, in the "Bray-Curtis k-mer
# dissimilarity" reference style: one row per group, a full min-max range
# line, an IQR box with median tick, and jittered raw points on top.
# Produced once per cutoff definition (>=1 Mb absolute, >=1% of genome) so
# the two can be compared side by side.
#
# Reads the small froh_comparison.csv that plot_msmc_roh.R already writes,
# so iterating on this plot's style never requires re-running the expensive
# per-sample scan (painting, contig-size-gap detection, .fai reads across
# every species) that produces that CSV.
#
# Deliberately a green/orange palette (not the red/green IUCN palette used
# for froh_by_iucn_status.png) so the two plot families read as distinct.
#
# Usage:
#   Rscript plot_froh_lc_vs_threatened.R [froh_csv] [iucn_csv] [out_dir]
#   Rscript plot_froh_lc_vs_threatened.R plots/roh/froh_comparison.csv /global/scratch/users/julesperez/snakemake_pipeline/metadata/complete_2026.csv plots/roh
#
# iucn_csv: "scientificName","redlistCategory" columns (e.g. complete_2026.csv).

args <- commandArgs(trailingOnly = TRUE)
froh_csv <- if (length(args) >= 1 && nzchar(args[1])) args[1] else "plots/roh/froh_comparison.csv"
iucn_csv <- if (length(args) >= 2 && nzchar(args[2])) args[2] else stop("iucn_csv is required (scientificName,redlistCategory columns)")
out_dir  <- if (length(args) >= 3 && nzchar(args[3])) args[3] else "plots/roh"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
froh_df <- read.csv(froh_csv, stringsAsFactors = FALSE)

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

froh_df$iucn <- sapply(froh_df$species, iucn_for_species)
threat_group <- ifelse(froh_df$iucn == "LC", "Non-threatened (LC)",
                        ifelse(is.na(froh_df$iucn), NA, "Threatened (NT/VU/EN/CR/EW)"))

CUTOFF_ORANGE <- "#D97B29"
CUTOFF_GREEN  <- "#4E8B3B"
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

plot_lc_vs_threatened <- function(values, out_path, x_label) {
  keep <- !is.na(threat_group) & is.finite(values)
  v <- values[keep]; g <- threat_group[keep]
  if (length(v) == 0 || length(unique(g)) < 2) {
    cat("Skipping", out_path, "-- not enough samples in both groups\n")
    return(invisible())
  }
  xr_rb <- range(v, na.rm = TRUE)
  xr_rb <- c(0, xr_rb[2] + 0.05 * diff(xr_rb))
  set.seed(1)
  png(out_path, width = 950, height = 560)
  tufte_par(mar = c(5, 11, 2, 2))
  plot(NA, xlim = xr_rb, ylim = c(0.4, 2.6), xaxt = "n", yaxt = "n",
       xlab = "", ylab = "", main = "")
  rangebox_row(2, v[g == "Threatened (NT/VU/EN/CR/EW)"], CUTOFF_ORANGE)
  rangebox_row(1, v[g == "Non-threatened (LC)"], CUTOFF_GREEN)
  axis(2, at = c(2, 1), labels = c("Threatened (NT/VU/EN/CR/EW)", "Non-threatened (LC)"),
       lwd = 0, cex.axis = 1, las = 1)
  axis(1, at = pretty(xr_rb), lwd = 0.6, cex.axis = 0.85)
  mtext(x_label, side = 1, line = 3, font = 2, cex = 1, col = "grey20")
  dev.off()
  cat("Wrote", out_path, "\n")
}

if ("pct_genome_in_roh_1mb" %in% names(froh_df)) {
  plot_lc_vs_threatened(froh_df$pct_genome_in_roh_1mb,
                         file.path(out_dir, "froh_lc_vs_threatened_1mb.png"),
                         "% of genome in ROH (FROH), ROH >= 1 Mb")
}
if ("pct_genome_in_roh_1pctgenome" %in% names(froh_df)) {
  plot_lc_vs_threatened(froh_df$pct_genome_in_roh_1pctgenome,
                         file.path(out_dir, "froh_lc_vs_threatened_1pctgenome.png"),
                         "% of genome in ROH (FROH), ROH >= 1% of genome")
}
