#!/usr/bin/env Rscript
# Harmonic mean Ne over fixed calendar-year time periods, for every finished
# PSMC and MSMC2 sample with a known generation time. Harmonic mean is the
# right summary for a coalescent Ne trajectory -- it weights toward periods
# of LOW Ne, which is exactly where most coalescence (and most information
# about a bottleneck) actually happens; an arithmetic mean would be
# dominated by whatever the highest Ne excursion happened to be.
#
# Reuses the exact same read_msmc()/read_psmc() parsing+scaling logic as
# plot_msmc_roh.R (same formulas, same boundary handling) -- duplicated here
# rather than sourced, matching the same self-contained-script pattern as
# plot_froh_cutoff_slopegraph.R / plot_froh_lc_vs_threatened.R. This script
# only reads small finished .psmc/msmc2.final.txt files directly (not the
# heavy per-sample .fai/painting scan that plot_msmc_roh.R does), so it's
# fast enough to run standalone.
#
# Usage:
#   Rscript analyze_ne_harmonic_means.R [results_dir] [mu] [gentime_csv] [samples_tsv] [chromosome_level_tsv] [sample_ids_file] [out_dir]
#
# Time periods (fixed calendar-year bins, same edges for every species):
#   0-10kya, 10k-100kya, 100k-1Mya, >1Mya -- plus "Overall" (full observed
#   range). A curve's own last time boundary is used as the effective upper
#   edge for both "Overall" and the >1Mya bin (rather than true infinity),
#   since PSMC/MSMC2's final interval is open-ended and integrating it out
#   to infinity would make that single interval's own Ne dominate the
#   harmonic mean completely, discarding everything else in the curve.

args <- commandArgs(trailingOnly = TRUE)
results_dir  <- if (length(args) >= 1 && nzchar(args[1])) args[1] else "results"
mu           <- if (length(args) >= 2 && nzchar(args[2])) as.numeric(args[2]) else 1.25e-8
gentime_csv  <- if (length(args) >= 3 && nzchar(args[3])) args[3] else stop("gentime_csv is required (phylo_name,gen_time columns)")
samples_tsv  <- if (length(args) >= 4 && nzchar(args[4])) args[4] else "../config/samples_medium.tsv"
arg_or_null <- function(i) if (length(args) >= i && nzchar(args[i])) args[i] else NULL
chromosome_level_tsv <- arg_or_null(5)
sample_ids <- { f <- arg_or_null(6); if (!is.null(f)) readLines(f) else NULL }
out_dir <- if (length(args) >= 7 && nzchar(args[7])) args[7] else "plots/ne_summary"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

MUTED_PALETTE <- c("#2F6B4F", "#7A9D54", "#C9A227", "#D9772E", "#B72E3C",
                    "#3D6B94", "#6B4C7A", "#8C8C8C")
# 300 DPI instead of R's default 72 -- scale pixel dimensions up
# proportionally to res so physical layout (margins, font size relative to
# plot) stays identical, only pixel density increases.
open_png <- function(path, width, height, res = 300) {
  scale <- res / 72
  png(path, width = width * scale, height = height * scale, res = res)
}

tufte_par <- function(mar = c(3, 4, 3, 2)) {
  par(bty = "n", family = "sans", las = 1, mar = mar,
      tck = -0.015, cex.axis = 0.85, col.axis = "grey30", col.lab = "grey20")
}

chrom_level_species <- if (!is.null(chromosome_level_tsv)) {
  read.table(chromosome_level_tsv, header = TRUE, sep = "\t", stringsAsFactors = FALSE)$Species
} else {
  NULL
}
is_chrom_level_species <- function(species) {
  if (is.null(chrom_level_species)) TRUE else species %in% chrom_level_species
}

gentimes <- read.csv(gentime_csv, stringsAsFactors = FALSE)
sample_sheet <- read.table(samples_tsv, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

gen_for_run <- function(run_id) {
  species <- sample_sheet$Species[sample_sheet$Run == run_id]
  if (length(species) == 0) return(NA)
  gt <- gentimes$gen_time[gentimes$phylo_name == species[1]]
  if (length(gt) == 0) return(NA)
  gt[1]
}
species_for_run <- function(run_id) {
  species <- sample_sheet$Species[sample_sheet$Run == run_id]
  if (length(species) == 0) run_id else species[1]
}

# --- Same conversions as plot_msmc_roh.R (verbatim) -------------------------
read_msmc <- function(path, mu, gen) {
  d <- read.table(path, header = TRUE, sep = "\t")
  x <- d[[2]] * gen / mu
  y <- (1 / d[[4]]) / (2 * mu)
  data.frame(x = x, y = y)
}
read_psmc <- function(path, mu, gen, bin_size = 100) {
  lines <- readLines(path)
  rd_idx <- grep("^RD", lines)
  if (length(rd_idx) == 0) return(NULL)
  block <- lines[rd_idx[length(rd_idx)]:length(lines)]
  tr_line <- block[grepl("^TR", block)][1]
  if (is.na(tr_line)) return(NULL)
  theta0 <- as.numeric(strsplit(tr_line, "\t")[[1]][2])
  rs_lines <- block[grepl("^RS", block)]
  if (length(rs_lines) == 0) return(NULL)
  parts <- strsplit(rs_lines, "\t")
  t_k <- sapply(parts, function(p) as.numeric(p[3]))
  lambda_k <- sapply(parts, function(p) as.numeric(p[4]))
  N0 <- theta0 / (4 * mu * bin_size)
  data.frame(x = 2 * N0 * t_k * gen, y = N0 * lambda_k)
}

# --- Harmonic mean of a step-function Ne curve over [lo, hi) ---------------
# x = interval start boundaries (ascending), y = Ne held constant from x[i]
# to x[i+1] (the last interval runs to +Inf in the raw data, but is capped
# at `hi` by the caller -- see note above).
harmonic_mean_ne <- function(x, y, lo, hi) {
  if (length(x) == 0 || !is.finite(hi) || hi <= lo) return(NA_real_)
  x_end <- c(x[-1], Inf)
  ov_lo <- pmax(x, lo)
  ov_hi <- pmin(x_end, hi)
  dur <- ov_hi - ov_lo
  dur[dur < 0] <- 0
  valid <- dur > 0 & is.finite(y) & y > 0
  if (!any(valid)) return(NA_real_)
  total_dur <- sum(dur[valid])
  total_dur / sum(dur[valid] / y[valid])
}

PERIOD_EDGES <- c(0, 1e4, 1e5, 1e6)  # lower edge of each of the 4 fixed bins
PERIOD_LABELS <- c("0-10kya", "10k-100kya", "100k-1Mya", ">1Mya")

periods_for_curve <- function(curve) {
  x_max <- max(curve$x)
  edges <- c(PERIOD_EDGES, x_max)
  rows <- lapply(seq_along(PERIOD_LABELS), function(i) {
    data.frame(period = PERIOD_LABELS[i],
               harmonic_mean_ne = harmonic_mean_ne(curve$x, curve$y, edges[i], edges[i + 1]))
  })
  rows[[length(rows) + 1]] <- data.frame(period = "Overall",
                                          harmonic_mean_ne = harmonic_mean_ne(curve$x, curve$y, 0, x_max))
  do.call(rbind, rows)
}

# --- Build the table for one method (msmc2 or psmc) -------------------------
build_method_table <- function(files, id_of, reader) {
  results <- list()
  for (f in files) {
    srr <- id_of(f)
    if (!is.null(sample_ids) && !(srr %in% sample_ids)) next
    species <- species_for_run(srr)
    if (!is_chrom_level_species(species)) next
    gen <- gen_for_run(srr)
    if (is.na(gen)) next
    curve <- reader(f, mu = mu, gen = gen)
    if (is.null(curve) || nrow(curve) == 0) next
    per <- periods_for_curve(curve)
    per$species <- species
    per$srr <- srr
    results[[srr]] <- per
  }
  if (length(results) == 0) return(NULL)
  do.call(rbind, results)
}

msmc_files <- Sys.glob(file.path(results_dir, "msmc2", "*", "msmc2.final.txt"))
msmc_table <- build_method_table(msmc_files, function(f) basename(dirname(f)), read_msmc)
if (!is.null(msmc_table)) {
  msmc_table$method <- "MSMC2"
  write.csv(msmc_table, file.path(out_dir, "ne_harmonic_means_msmc2.csv"), row.names = FALSE)
  cat("Wrote", file.path(out_dir, "ne_harmonic_means_msmc2.csv"), "(", length(unique(msmc_table$srr)), "samples)\n")
} else {
  cat("No MSMC2 samples with a known generation time -- skipping.\n")
}

psmc_files <- Sys.glob(file.path(results_dir, "psmc", "*.psmc"))
psmc_table <- build_method_table(psmc_files, function(f) sub("\\.psmc$", "", basename(f)), read_psmc)
if (!is.null(psmc_table)) {
  psmc_table$method <- "PSMC"
  write.csv(psmc_table, file.path(out_dir, "ne_harmonic_means_psmc.csv"), row.names = FALSE)
  cat("Wrote", file.path(out_dir, "ne_harmonic_means_psmc.csv"), "(", length(unique(psmc_table$srr)), "samples)\n")
} else {
  cat("No PSMC samples with a known generation time -- skipping.\n")
}

# --- Colored block/range-box plot: one row per period, harmonic mean Ne
# spread across species -- same visual grammar as the FROH range-box plots
# (full min-max line, IQR box, median tick, jittered points), one color per
# row so the periods are visually distinguishable at a glance.
plot_periods <- function(tbl, out_path, title) {
  if (is.null(tbl)) {
    cat("Skipping", out_path, "-- no data\n")
    return(invisible())
  }
  row_order <- c("Overall", rev(PERIOD_LABELS))
  cols <- setNames(muted <- colorRampPalette(MUTED_PALETTE)(length(row_order)), row_order)
  v_all <- tbl$harmonic_mean_ne[is.finite(tbl$harmonic_mean_ne) & tbl$harmonic_mean_ne > 0]
  if (length(v_all) == 0) {
    cat("Skipping", out_path, "-- no finite values\n")
    return(invisible())
  }
  xr <- range(v_all)
  n <- length(row_order)
  open_png(out_path, width = 950, height = 560)
  tufte_par(mar = c(4, 11, 3, 2))
  plot(NA, xlim = xr, ylim = c(0.4, n + 0.6), log = "x", xaxt = "n", yaxt = "n",
       xlab = "", ylab = "", main = "")
  axis(1, lwd = 0.6, cex.axis = 0.85)
  mtext(expression(N[e]~"(harmonic mean)"), side = 1, line = 2.6, font = 2, cex = 1, col = "grey20")
  mtext(title, side = 3, line = 1, adj = 0, cex = 1.05, font = 2, col = "grey15")
  for (i in seq_along(row_order)) {
    p <- row_order[i]
    y <- n - i + 1
    v <- tbl$harmonic_mean_ne[tbl$period == p & is.finite(tbl$harmonic_mean_ne) & tbl$harmonic_mean_ne > 0]
    if (length(v) == 0) next
    col <- cols[[p]]
    rng <- range(v)
    segments(rng[1], y, rng[2], y, col = col, lwd = 1)
    if (length(v) >= 2) {
      qs <- quantile(v, c(0.25, 0.5, 0.75))
      rect(qs[1], y - 0.28, qs[3], y + 0.28, col = adjustcolor(col, alpha.f = 0.18), border = col, lwd = 1.2)
      segments(qs[2], y - 0.28, qs[2], y + 0.28, col = col, lwd = 2.2)
    }
    yj <- y + (stats::runif(length(v)) - 0.5) * 0.5
    points(v, yj, pch = 21, bg = adjustcolor(col, alpha.f = 0.55),
           col = adjustcolor("black", alpha.f = 0.5), cex = 1.1, lwd = 0.5)
  }
  axis(2, at = seq_len(n), labels = rev(row_order), lwd = 0, cex.axis = 0.95, las = 1)
  dev.off()
  cat("Wrote", out_path, "\n")
}

plot_periods(msmc_table, file.path(out_dir, "ne_harmonic_by_period_msmc2.png"),
             "MSMC2 -- harmonic mean Ne by time period")
plot_periods(psmc_table, file.path(out_dir, "ne_harmonic_by_period_psmc.png"),
             "PSMC -- harmonic mean Ne by time period")
