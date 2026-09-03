#!/usr/bin/env Rscript
# Same-day progress plots for whatever ROH/MSMC2 samples have finished so far.
# Base R only (no ggplot2/tidyverse) since no R conda env exists in this repo yet.
#
# Usage:
#   Rscript plot_msmc_roh.R [results_dir] [mu] [gentime_csv] [sample_ids_file] [samples_tsv]
#   Rscript plot_msmc_roh.R results 1.25e-8 gentimes.csv top10_ids.txt ../config/samples_medium.tsv
#
# gentime_csv: two columns "species,gentime" (header required) -- generation
# time varies per species, so each sample's years-ago axis is scaled using
# its own species' value, looked up via samples_tsv's Run->Species mapping.
#
# sample_ids_file (optional): one Run ID per line -- restricts plotting to
# just those samples instead of every finished one found under results_dir.
# Generate one for e.g. the first 10 samples in the sheet with:
#   tail -n +2 ../config/samples_medium.tsv | cut -f1 | head -10 > top10_ids.txt

args <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1) args[1] else "results"
mu           <- if (length(args) >= 2) as.numeric(args[2]) else 1.25e-8
gentime_csv  <- if (length(args) >= 3) args[3] else stop("gentime_csv is required (species,gentime columns)")
sample_ids   <- if (length(args) >= 4) readLines(args[4]) else NULL
samples_tsv  <- if (length(args) >= 5) args[5] else "../config/samples_medium.tsv"

gentimes <- read.csv(gentime_csv, stringsAsFactors = FALSE)
sample_sheet <- read.table(samples_tsv, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Run ID -> generation time, via Run -> Species -> gentime
gen_for_run <- function(run_id) {
  species <- sample_sheet$Species[sample_sheet$Run == run_id]
  if (length(species) == 0) {
    warning(paste0("No Species found for Run ", run_id, " in ", samples_tsv))
    return(NA)
  }
  gt <- gentimes$gentime[gentimes$species == species[1]]
  if (length(gt) == 0) {
    warning(paste0("No gentime found for species '", species[1], "' (sample ", run_id, ") in ", gentime_csv))
    return(NA)
  }
  gt[1]
}

dir.create("plots/msmc2", recursive = TRUE, showWarnings = FALSE)
dir.create("plots/roh", recursive = TRUE, showWarnings = FALSE)

# --- MSMC2 ---
# Same conversion as msmc-tools' plot_utils.py popSizeStepPlot: x = left_time
# boundary scaled to years via generation time / mu, y = Ne from the
# coalescence rate (first lambda column). Read by position, not by column
# name -- msmc2's single-sample output names that column "lambda", not
# "lambda_00" (the "_00" naming only shows up for multi-population runs).
read_msmc <- function(path, mu, gen) {
  d <- read.table(path, header = TRUE, sep = "\t")
  x <- d[[2]] * gen / mu       # column 2 = left_time_boundary
  y <- (1 / d[[4]]) / (2 * mu) # column 4 = first lambda column
  data.frame(x = x, y = y)
}

msmc_files <- Sys.glob(file.path(results_dir, "msmc2", "*", "msmc2.final.txt"))
if (!is.null(sample_ids)) {
  msmc_files <- msmc_files[basename(dirname(msmc_files)) %in% sample_ids]
}
if (length(msmc_files) == 0) {
  cat("No finished msmc2.final.txt files found.\n")
} else {
  srr_ids <- basename(dirname(msmc_files))
  gens <- setNames(sapply(srr_ids, gen_for_run), srr_ids)
  keep <- !is.na(gens)
  if (any(!keep)) {
    cat("Skipping (no gentime found):", paste(srr_ids[!keep], collapse = ", "), "\n")
  }
  msmc_files <- msmc_files[keep]
  srr_ids <- srr_ids[keep]
  gens <- gens[keep]

  curves <- Map(function(f, g) read_msmc(f, mu = mu, gen = g), msmc_files, gens)
  names(curves) <- srr_ids

  # individual plots
  for (srr in srr_ids) {
    d <- curves[[srr]]
    png(file.path("plots/msmc2", paste0(srr, ".png")), width = 900, height = 700)
    plot(d$x, d$y, type = "s", log = "xy",
         xlab = "Years ago", ylab = "Effective population size (Ne)",
         main = paste0("MSMC2 -- ", srr, " (gen=", gens[[srr]], ")"))
    dev.off()
  }

  # overlay plot, all finished samples together
  cols <- rainbow(length(curves))
  xr <- range(unlist(lapply(curves, function(d) d$x[d$x > 0])))
  yr <- range(unlist(lapply(curves, function(d) d$y)))
  png("plots/msmc2/all_samples_overlay.png", width = 1100, height = 800)
  plot(NA, xlim = xr, ylim = yr, log = "xy",
       xlab = "Years ago", ylab = "Effective population size (Ne)",
       main = paste0("MSMC2 -- ", length(curves), " sample(s) completed"))
  for (i in seq_along(curves)) {
    lines(curves[[i]]$x, curves[[i]]$y, type = "s", col = cols[i])
  }
  legend("topright", legend = names(curves), col = cols, lty = 1, cex = 0.6, ncol = 2)
  dev.off()

  cat("Wrote", length(msmc_files), "individual MSMC2 plots + 1 overlay to plots/msmc2/\n")
}

# --- ROH ---
# bcftools roh output: "RG" data lines have columns
# RG, Sample, Chromosome, Start, End, Length(bp), nMarkers, Quality
roh_files <- Sys.glob(file.path(results_dir, "roh", "*_ROH.txt"))
if (!is.null(sample_ids)) {
  roh_files <- roh_files[sub("_ROH.txt$", "", basename(roh_files)) %in% sample_ids]
}
if (length(roh_files) == 0) {
  cat("No finished *_ROH.txt files found.\n")
} else {
  all_roh <- do.call(rbind, lapply(roh_files, function(f) {
    srr <- sub("_ROH.txt$", "", basename(f))
    lines <- readLines(f)
    rg <- lines[startsWith(lines, "RG")]
    if (length(rg) == 0) return(NULL)
    parts <- strsplit(rg, "\t")
    data.frame(
      srr = srr,
      chrom = sapply(parts, `[`, 3),
      length_bp = as.numeric(sapply(parts, `[`, 6))
    )
  }))

  if (is.null(all_roh) || nrow(all_roh) == 0) {
    cat("ROH files found but contained no RG segments.\n")
  } else {
    # top 30 scaffolds by total ROH length summed across all finished samples
    by_chrom <- aggregate(length_bp ~ chrom, data = all_roh, sum)
    by_chrom <- by_chrom[order(-by_chrom$length_bp), ]
    top30 <- head(by_chrom, 30)

    png("plots/roh/top30_scaffolds_total_roh.png", width = 1200, height = 700)
    par(mar = c(8, 5, 4, 2))
    barplot(top30$length_bp / 1e6, names.arg = top30$chrom, las = 2, cex.names = 0.7,
            ylab = "Total ROH length (Mb)",
            main = paste0("Top ", nrow(top30), " scaffolds by total ROH -- ",
                           length(unique(all_roh$srr)), " sample(s)"))
    dev.off()

    # per-sample total ROH, for a quick sample-level overview
    by_sample <- aggregate(length_bp ~ srr, data = all_roh, sum)
    by_sample <- by_sample[order(-by_sample$length_bp), ]
    png("plots/roh/total_roh_per_sample.png",
        width = max(900, nrow(by_sample) * 40), height = 700)
    par(mar = c(10, 5, 4, 2))
    barplot(by_sample$length_bp / 1e6, names.arg = by_sample$srr, las = 2, cex.names = 0.7,
            ylab = "Total ROH length (Mb)",
            main = paste0("Total ROH per sample -- ", nrow(by_sample), " sample(s) completed"))
    dev.off()

    cat("Wrote plots/roh/top30_scaffolds_total_roh.png and plots/roh/total_roh_per_sample.png\n")
  }
}
