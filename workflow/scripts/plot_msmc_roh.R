#!/usr/bin/env Rscript
# Same-day progress plots for whatever ROH/MSMC2 samples have finished so far.
# Base R only (no ggplot2/tidyverse) since no R conda env exists in this repo yet.
#
# Usage:
#   Rscript plot_msmc_roh.R [results_dir] [mu] [gentime_csv] [sample_ids_file] [samples_tsv]
#   Rscript plot_msmc_roh.R results 1.25e-8 gentimes.csv top10_ids.txt ../config/samples_medium.tsv
#
# gentime_csv: expects "phylo_name","gen_time" columns (e.g.
# gen_time_estimates.csv) -- generation time varies per species, so each
# sample's years-ago axis is scaled using its own species' value, looked up
# via samples_tsv's Run->Species mapping matched against phylo_name.
#
# sample_ids_file (optional): one Run ID per line -- restricts plotting to
# just those samples instead of every finished one found under results_dir.
# Generate one for e.g. the first 10 samples in the sheet with:
#   tail -n +2 ../config/samples_medium.tsv | cut -f1 | head -10 > top10_ids.txt
#
# ref_genome_path (optional): matches config.yaml's reference_genome_path --
# needed for chromosome-painting/FROH, which read each species' .fasta.fai
# for real scaffold lengths (bcftools roh output alone has no genome-length
# info, only ROH segment coordinates).

args <- commandArgs(trailingOnly = TRUE)
results_dir  <- if (length(args) >= 1) args[1] else "results"
mu           <- if (length(args) >= 2) as.numeric(args[2]) else 1.25e-8
gentime_csv  <- if (length(args) >= 3) args[3] else stop("gentime_csv is required (species,gentime columns)")
sample_ids   <- if (length(args) >= 4) readLines(args[4]) else NULL
samples_tsv  <- if (length(args) >= 5) args[5] else "../config/samples_medium.tsv"
ref_genome_path <- if (length(args) >= 6) args[6] else "/global/scratch/projects/fc_moilab/julesperez/post_rot/new_refgenomes/"

gentimes <- read.csv(gentime_csv, stringsAsFactors = FALSE)
sample_sheet <- read.table(samples_tsv, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Run ID -> generation time, via Run -> Species -> gentime
gen_for_run <- function(run_id) {
  species <- sample_sheet$Species[sample_sheet$Run == run_id]
  if (length(species) == 0) {
    warning(paste0("No Species found for Run ", run_id, " in ", samples_tsv))
    return(NA)
  }
  gt <- gentimes$gen_time[gentimes$phylo_name == species[1]]
  if (length(gt) == 0) {
    warning(paste0("No gentime found for species '", species[1], "' (sample ", run_id, ") in ", gentime_csv))
    return(NA)
  }
  gt[1]
}

species_for_run <- function(run_id) {
  species <- sample_sheet$Species[sample_sheet$Run == run_id]
  if (length(species) == 0) run_id else species[1]
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
  species_labels <- sapply(srr_ids, species_for_run)

  curves <- Map(function(f, g) read_msmc(f, mu = mu, gen = g), msmc_files, gens)
  names(curves) <- srr_ids

  # individual plots -- filenames/titles/legend all use species name, not
  # the Run ID, per request. Falls back to srr_id only if a species somehow
  # maps to itself (see species_for_run).
  for (srr in srr_ids) {
    d <- curves[[srr]]
    label <- species_labels[[srr]]
    png(file.path("plots/msmc2", paste0(label, ".png")), width = 900, height = 700)
    plot(d$x, d$y, type = "s", log = "xy",
         xlab = "Years ago", ylab = "Effective population size (Ne)",
         main = paste0("MSMC2 -- ", label, " (gen=", gens[[srr]], ")"))
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
  legend("topright", legend = species_labels, col = cols, lty = 1, cex = 0.6, ncol = 2)
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
      species = species_for_run(srr),
      chrom = sapply(parts, `[`, 3),
      start = as.numeric(sapply(parts, `[`, 4)),
      end = as.numeric(sapply(parts, `[`, 5)),
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

    # per-sample total ROH, for a quick sample-level overview -- labeled by
    # species name, not Run ID, same as the MSMC2 plots above
    by_sample <- aggregate(length_bp ~ species, data = all_roh, sum)
    by_sample <- by_sample[order(-by_sample$length_bp), ]
    png("plots/roh/total_roh_per_sample.png",
        width = max(900, nrow(by_sample) * 40), height = 700)
    par(mar = c(10, 5, 4, 2))
    barplot(by_sample$length_bp / 1e6, names.arg = by_sample$species, las = 2, cex.names = 0.7,
            ylab = "Total ROH length (Mb)",
            main = paste0("Total ROH per sample -- ", nrow(by_sample), " sample(s) completed"))
    dev.off()

    cat("Wrote plots/roh/top30_scaffolds_total_roh.png and plots/roh/total_roh_per_sample.png\n")

    # --- Chromosome painting + genome-wide FROH ---
    # bcftools roh's own output has no scaffold-length info, only ROH
    # segment coordinates -- pull real lengths from each species' .fasta.fai
    # (the same file jules_msmc2_contigs/jules_samtools_faidx already
    # generate) so both plots are scaled to real genome size, not just
    # segment counts.
    dir.create("plots/roh/painting", recursive = TRUE, showWarnings = FALSE)
    fai_for_species <- function(species) {
      file.path(ref_genome_path, species, paste0(species, ".fasta.fai"))
    }

    froh_rows <- list()
    for (srr in unique(all_roh$srr)) {
      species <- unique(all_roh$species[all_roh$srr == srr])[1]
      fai_path <- fai_for_species(species)
      if (!file.exists(fai_path)) {
        warning(paste0("No .fai found for species '", species, "' (sample ", srr, ") at ", fai_path))
        next
      }
      fai <- read.table(fai_path, sep = "\t", stringsAsFactors = FALSE)
      colnames(fai)[1:2] <- c("chrom", "len")
      # >500000bp matches jules_msmc2_contigs' own real-chromosome cutoff --
      # painting every scaffold (not just a fixed top-N) is what was asked
      # for, but a raw .fai can carry thousands of tiny unplaced fragments,
      # which is both unreadable and can blow past Cairo's ~32767px max PNG
      # dimension (confirmed: this failed outright on a fragmented assembly
      # before this filter was added).
      top_scaffolds <- fai[fai$len > 500000, ]
      top_scaffolds <- top_scaffolds[order(-top_scaffolds$len), ]
      if (nrow(top_scaffolds) == 0) {
        top_scaffolds <- fai[order(-fai$len), ][1, ]
      }

      sample_roh <- all_roh[all_roh$srr == srr, ]

      # chromosome painting: one horizontal track per scaffold (grey = full
      # length), red blocks = ROH segments drawn to scale within it.
      # Height clamped to Cairo's device limit as a hard safety net even
      # after the >500000bp filter above.
      png(file.path("plots/roh/painting", paste0(species, "_painting.png")),
          width = 1200, height = min(30000, max(400, nrow(top_scaffolds) * 25)))
      par(mar = c(5, 10, 4, 2))
      plot(NA, xlim = c(0, max(top_scaffolds$len)), ylim = c(0, nrow(top_scaffolds) + 1),
           yaxt = "n", xlab = "Position (bp)", ylab = "",
           main = paste0("ROH painting -- ", species, " (", nrow(top_scaffolds), " scaffolds >500kb)"))
      axis(2, at = seq_len(nrow(top_scaffolds)), labels = rev(top_scaffolds$chrom),
           las = 2, cex.axis = 0.6)
      for (i in seq_len(nrow(top_scaffolds))) {
        y <- nrow(top_scaffolds) - i + 1
        chrom <- top_scaffolds$chrom[i]
        rect(0, y - 0.3, top_scaffolds$len[i], y + 0.3, col = "grey85", border = NA)
        segs <- sample_roh[sample_roh$chrom == chrom, ]
        if (nrow(segs) > 0) {
          rect(segs$start, y - 0.3, segs$end, y + 0.3, col = "red", border = NA)
        }
      }
      dev.off()

      froh_rows[[srr]] <- data.frame(
        species = species,
        genome_len_bp = sum(fai$len),
        roh_len_bp = sum(sample_roh$length_bp)
      )
    }

    if (length(froh_rows) > 0) {
      froh_df <- do.call(rbind, froh_rows)
      froh_df$pct_genome_in_roh <- 100 * froh_df$roh_len_bp / froh_df$genome_len_bp
      froh_df <- froh_df[order(-froh_df$pct_genome_in_roh), ]

      png("plots/roh/froh_percent_genome.png",
          width = max(900, nrow(froh_df) * 60), height = 700)
      par(mar = c(10, 5, 4, 2))
      barplot(froh_df$pct_genome_in_roh, names.arg = froh_df$species, las = 2, cex.names = 0.8,
              ylab = "% of genome in ROH (FROH)",
              main = paste0("Genome-wide FROH -- ", nrow(froh_df), " sample(s)"))
      dev.off()

      cat("Wrote", length(froh_rows), "chromosome-painting plots to plots/roh/painting/",
          "and plots/roh/froh_percent_genome.png\n")
    }
  }
}
