#!/usr/bin/env Rscript
# Tufte slopegraph comparing two FROH segment-inclusion cutoffs (>=1% of
# chromosome vs >=1 Mb absolute). Reads the small froh_comparison.csv that
# plot_msmc_roh.R already writes, so iterating on this plot's style never
# requires re-running the expensive per-sample scan (painting, contig-size-
# gap detection, .fai reads across every species) that produces that CSV.
#
# Usage:
#   Rscript plot_froh_cutoff_slopegraph.R [froh_csv] [out_dir]
#   Rscript plot_froh_cutoff_slopegraph.R plots/roh/froh_comparison.csv plots/roh

args <- commandArgs(trailingOnly = TRUE)
froh_csv <- if (length(args) >= 1 && nzchar(args[1])) args[1] else "plots/roh/froh_comparison.csv"
out_dir  <- if (length(args) >= 2 && nzchar(args[2])) args[2] else "plots/roh"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
froh_df <- read.csv(froh_csv, stringsAsFactors = FALSE)

MUTED_PALETTE <- c("#2F6B4F", "#7A9D54", "#C9A227", "#D9772E", "#B72E3C",
                    "#3D6B94", "#6B4C7A", "#8C8C8C")
tufte_par <- function(mar = c(3, 4, 3, 2)) {
  par(bty = "n", family = "sans", las = 1, mar = mar,
      tck = -0.015, cex.axis = 0.85, col.axis = "grey30", col.lab = "grey20")
}
tufte_title <- function(main) {
  mtext(main, side = 3, line = 1, adj = 0, cex = 1.05, font = 2, col = "grey15")
}

stopifnot(all(c("pct_genome_in_roh_1pctchrom", "pct_genome_in_roh_1mb") %in% names(froh_df)))

# One dot per species per cutoff, connected by a line, species named
# directly at both ends instead of a legend -- the two vertical dot
# columns plus their slopes ARE the comparison.
cmp <- froh_df[order(-froh_df$pct_genome_in_roh_1mb), ]
x1 <- 1; x2 <- 2
yr_cmp <- range(c(cmp$pct_genome_in_roh_1pctchrom, cmp$pct_genome_in_roh_1mb), na.rm = TRUE)
label_ok <- nrow(cmp) <= 30
png(file.path(out_dir, "froh_cutoff_comparison.png"),
    width = if (label_ok) 1100 else 700, height = max(500, nrow(cmp) * 16))
tufte_par(mar = c(3, 3, 4, 3))
plot(NA, xlim = c(if (label_ok) 0.3 else 0.8, if (label_ok) 2.7 else 2.2), ylim = yr_cmp,
     xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = "")
axis(2, lwd = 0.6)
axis(1, at = c(x1, x2), labels = c("ROH >= 1% of chromosome", "ROH >= 1 Mb"),
     lwd = 0, cex.axis = 0.9, padj = -1)
tufte_title(paste0("FROH -- relative (1% of chromosome) vs absolute (1 Mb) segment cutoff, ",
                    nrow(cmp), " sample(s)"))
for (i in seq_len(nrow(cmp))) {
  y1 <- cmp$pct_genome_in_roh_1pctchrom[i]
  y2 <- cmp$pct_genome_in_roh_1mb[i]
  segments(x1, y1, x2, y2, col = adjustcolor(MUTED_PALETTE[1], alpha.f = 0.45), lwd = 1)
  points(c(x1, x2), c(y1, y2), pch = 16, col = MUTED_PALETTE[1], cex = 0.8)
}
if (label_ok) {
  text(x1 - 0.05, cmp$pct_genome_in_roh_1pctchrom, cmp$species, adj = 1, cex = 0.55, col = "grey30")
  text(x2 + 0.05, cmp$pct_genome_in_roh_1mb, cmp$species, adj = 0, cex = 0.55, col = "grey30")
}
dev.off()
cat("Wrote", file.path(out_dir, "froh_cutoff_comparison.png"), "\n")
