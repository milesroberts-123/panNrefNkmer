#!/usr/bin/env Rscript
# Same-day progress plots for whatever ROH/MSMC2/PSMC samples have finished
# so far. Base R only (no ggplot2/tidyverse) since no R conda env exists in
# this repo yet. Styled loosely after Tufte's "Visual Display of Quantitative
# Information": no chart box, minimal sparse ticks, direct labeling instead
# of legends where practical, muted (not rainbow) color, ink kept
# proportional to information.
#
# Usage:
#   Rscript plot_msmc_roh.R [results_dir] [mu] [gentime_csv] [sample_ids_file] [samples_tsv] [ref_genome_path] [chromosome_level_tsv] [iucn_csv]
#
# gentime_csv: expects "phylo_name","gen_time" columns (e.g.
# gen_time_estimates.csv) -- generation time varies per species, so each
# sample's years-ago axis is scaled using its own species' value, looked up
# via samples_tsv's Run->Species mapping matched against phylo_name.
#
# sample_ids_file (optional): one Run ID per line -- restricts plotting to
# just those samples instead of every finished one found under results_dir.
#
# ref_genome_path (optional): matches config.yaml's reference_genome_path --
# needed for chromosome-painting/FROH, which read each species' .fasta.fai
# for real scaffold lengths.
#
# chromosome_level_tsv (optional): "Species","SampleID" columns -- the
# ground-truth list of species with an actual chromosome-level assembly.
# Applied to every plot so a scaffold-level species never silently ends up
# presented as if it were chromosome-level.
#
# iucn_csv (optional): "scientificName","redlistCategory" columns (e.g.
# complete_2026.csv). redlistCategory is free text ("Endangered",
# "Lower Risk/least concern", etc.) and gets normalized to the standard
# CR/EN/VU/NT/LC/EW/EX/DD codes. Needed for the FROH-by-IUCN-status plot.

args <- commandArgs(trailingOnly = TRUE)
results_dir  <- if (length(args) >= 1) args[1] else "results"
mu           <- if (length(args) >= 2) as.numeric(args[2]) else 1.25e-8
gentime_csv  <- if (length(args) >= 3) args[3] else stop("gentime_csv is required (species,gentime columns)")
arg_or_null <- function(i) if (length(args) >= i && nzchar(args[i])) args[i] else NULL
sample_ids   <- { f <- arg_or_null(4); if (!is.null(f)) readLines(f) else NULL }
samples_tsv  <- if (length(args) >= 5 && nzchar(args[5])) args[5] else "../config/samples_medium.tsv"
ref_genome_path <- if (length(args) >= 6 && nzchar(args[6])) args[6] else "/global/scratch/projects/fc_moilab/julesperez/post_rot/new_refgenomes/"
chromosome_level_tsv <- arg_or_null(7)
iucn_csv     <- arg_or_null(8)

# --- Tufte-ish shared style ---------------------------------------------
IUCN_ORDER  <- c("CR", "EN", "VU", "NT", "LC")
IUCN_COLORS <- c(CR = "#B72E3C", EN = "#D9772E", VU = "#C9A227",
                  NT = "#7A9D54", LC = "#2F6B4F")
MUTED_PALETTE <- c("#2F6B4F", "#7A9D54", "#C9A227", "#D9772E", "#B72E3C",
                    "#3D6B94", "#6B4C7A", "#8C8C8C")
muted_colors <- function(n) {
  if (n <= length(MUTED_PALETTE)) return(MUTED_PALETTE[seq_len(n)])
  colorRampPalette(MUTED_PALETTE)(n)
}
tufte_par <- function(mar = c(3, 4, 3, 2)) {
  par(bty = "n", family = "sans", las = 1, mar = mar,
      tck = -0.015, cex.axis = 0.85, col.axis = "grey30", col.lab = "grey20")
}
tufte_title <- function(main) {
  mtext(main, side = 3, line = 1, adj = 0, cex = 1.05, font = 2, col = "grey15")
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

# --- IUCN status lookup ---------------------------------------------------
# redlistCategory is free text -- normalize "Lower Risk/least concern" etc.
# (older IUCN category names) down to the standard 5-letter codes.
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

iucn_lookup_df <- if (!is.null(iucn_csv)) read.csv(iucn_csv, stringsAsFactors = FALSE) else NULL
iucn_for_species <- function(species) {
  if (is.null(iucn_lookup_df)) return(NA_character_)
  sci <- gsub("_", " ", species)
  row <- iucn_lookup_df[tolower(iucn_lookup_df$scientificName) == tolower(sci), ]
  if (nrow(row) == 0) return(NA_character_)
  normalize_iucn(row$redlistCategory[1])
}

dir.create("plots/msmc2", recursive = TRUE, showWarnings = FALSE)
dir.create("plots/roh", recursive = TRUE, showWarnings = FALSE)
dir.create("plots/comparison", recursive = TRUE, showWarnings = FALSE)

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

# --- PSMC ---
# Standard psmc_plot.pl conversion: N0 = theta0 / (4 * mu * bin_size), using
# the LAST restart ("RD") block's TR (theta0) and RS (t_k, lambda_k) lines.
# bin_size matches jules_psmc_run_psmc's default -s100 (unset in the rule ->
# psmc's own default of 100bp bins). Time/Ne are rescaled to real units the
# same way psmc_plot.pl does: t (years) = 2*N0*t_k*gen, Ne = N0*lambda_k.
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

msmc_files <- Sys.glob(file.path(results_dir, "msmc2", "*", "msmc2.final.txt"))
if (!is.null(sample_ids)) {
  msmc_files <- msmc_files[basename(dirname(msmc_files)) %in% sample_ids]
}
not_chrom_level <- msmc_files[!sapply(basename(dirname(msmc_files)), function(s) is_chrom_level_species(species_for_run(s)))]
if (length(not_chrom_level) > 0) {
  cat("Excluding (not chromosome-level assembly):", paste(basename(dirname(not_chrom_level)), collapse = ", "), "\n")
}
msmc_files <- msmc_files[sapply(basename(dirname(msmc_files)), function(s) is_chrom_level_species(species_for_run(s)))]

curves <- list()
srr_ids <- character(0)
species_labels <- character(0)
gens <- numeric(0)

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

  # individual plots -- filenames/titles use species name, not the Run ID.
  for (srr in srr_ids) {
    d <- curves[[srr]]
    label <- species_labels[[srr]]
    png(file.path("plots/msmc2", paste0(label, ".png")), width = 900, height = 650)
    tufte_par()
    plot(d$x, d$y, type = "s", log = "xy", lwd = 1.6, col = MUTED_PALETTE[1],
         xlab = "Years ago", ylab = expression(N[e]), main = "", axes = FALSE)
    axis(1, lwd = 0.6); axis(2, lwd = 0.6)
    tufte_title(paste0(label, "  (MSMC2, generation time = ", round(gens[[srr]], 1), "y)"))
    dev.off()
  }

  # overlay plot, all finished samples together. A per-species legend only
  # stays legible up to ~15 curves -- past that, individually identifying a
  # species by line color is a lost cause regardless of palette, so instead
  # show the ensemble as a thin translucent "cloud" (the actual informative
  # signal at this n is the shape of the bundle, not any one species' line)
  # and point the reader to the per-species small multiples already in
  # plots/msmc2/ for individual identification.
  LEGEND_MAX <- 15
  xr <- range(unlist(lapply(curves, function(d) d$x[d$x > 0])))
  yr <- range(unlist(lapply(curves, function(d) d$y)))
  png("plots/msmc2/all_samples_overlay.png", width = 1100, height = 800)
  tufte_par()
  plot(NA, xlim = xr, ylim = yr, log = "xy",
       xlab = "Years ago", ylab = expression(N[e]), main = "", axes = FALSE)
  axis(1, lwd = 0.6); axis(2, lwd = 0.6)
  if (length(curves) <= LEGEND_MAX) {
    cols <- muted_colors(length(curves))
    tufte_title(paste0("MSMC2 -- ", length(curves), " sample(s) completed"))
    for (i in seq_along(curves)) {
      lines(curves[[i]]$x, curves[[i]]$y, type = "s", col = cols[i], lwd = 1.4)
    }
    legend("topright", legend = species_labels, col = cols, lty = 1, lwd = 1.4,
           cex = 0.6, ncol = 2, bty = "n")
  } else {
    tufte_title(paste0("MSMC2 -- ", length(curves),
                        " samples (individual plots in plots/msmc2/)"))
    for (i in seq_along(curves)) {
      lines(curves[[i]]$x, curves[[i]]$y, type = "s",
            col = adjustcolor(MUTED_PALETTE[1], alpha.f = 0.25), lwd = 1)
    }
  }
  dev.off()

  cat("Wrote", length(msmc_files), "individual MSMC2 plots + 1 overlay to plots/msmc2/\n")
}

# --- PSMC vs MSMC2 comparison ---
# Only for samples where both a finished .psmc and msmc2.final.txt exist --
# auto-discovered, no need to hand-pick a sample.
psmc_files <- Sys.glob(file.path(results_dir, "psmc", "*.psmc"))
psmc_ids <- sub("\\.psmc$", "", basename(psmc_files))
if (!is.null(sample_ids)) psmc_ids <- psmc_ids[psmc_ids %in% sample_ids]
psmc_ids <- psmc_ids[sapply(psmc_ids, function(s) is_chrom_level_species(species_for_run(s)))]

common_ids <- intersect(srr_ids, psmc_ids)
if (length(common_ids) == 0) {
  cat("No samples with both finished PSMC and MSMC2 output yet -- skipping comparison plots.\n")
} else {
  for (srr in common_ids) {
    label <- species_labels[[srr]]
    psmc_path <- file.path(results_dir, "psmc", paste0(srr, ".psmc"))
    pd <- read_psmc(psmc_path, mu = mu, gen = gens[[srr]])
    md <- curves[[srr]]
    if (is.null(pd) || is.null(md)) next

    x_vals <- c(pd$x[pd$x > 0], md$x[md$x > 0])
    y_vals <- c(pd$y[pd$y > 0], md$y[md$y > 0])
    x_vals <- x_vals[is.finite(x_vals)]
    y_vals <- y_vals[is.finite(y_vals)]
    if (length(x_vals) == 0 || length(y_vals) == 0) {
      cat("Skipping comparison plot for", label, "-- non-finite/empty PSMC or MSMC2 values (degenerate fit, e.g. theta0=0)\n")
      next
    }
    xr <- range(x_vals)
    # PSMC's earliest/latest time bins are known boundary artifacts that can
    # blow up Ne by 1-2 orders of magnitude -- left uncorrected, that alone
    # sets the axis scale and squashes the actual informative part of the
    # curve into a sliver at the bottom. Clip to the 1st-99th percentile of
    # observed Ne instead of the raw min/max; the underlying data points are
    # still drawn (and can run past the frame), only the axis range is robust.
    yr <- quantile(y_vals, c(0.01, 0.99), na.rm = TRUE)
    if (!is.finite(diff(yr)) || diff(yr) <= 0) yr <- range(y_vals)
    png(file.path("plots/comparison", paste0(label, "_psmc_vs_msmc2.png")), width = 950, height = 700)
    tufte_par()
    plot(NA, xlim = xr, ylim = yr, log = "xy",
         xlab = "Years ago", ylab = expression(N[e]), main = "", axes = FALSE)
    axis(1, lwd = 0.6); axis(2, lwd = 0.6)
    tufte_title(paste0(label, " -- PSMC vs MSMC2"))
    lines(pd$x, pd$y, type = "s", col = MUTED_PALETTE[5], lwd = 1.6, lty = 2)
    lines(md$x, md$y, type = "s", col = MUTED_PALETTE[1], lwd = 1.6, lty = 1)
    legend("topright", legend = c("PSMC", "MSMC2"), col = c(MUTED_PALETTE[5], MUTED_PALETTE[1]),
           lty = c(2, 1), lwd = 1.6, bty = "n", cex = 0.85)
    dev.off()
  }
  cat("Wrote", length(common_ids), "PSMC-vs-MSMC2 comparison plots to plots/comparison/\n")
}

# --- ROH ---
# bcftools roh output: "RG" data lines have columns
# RG, Sample, Chromosome, Start, End, Length(bp), nMarkers, Quality
roh_files <- Sys.glob(file.path(results_dir, "roh", "*_ROH.txt"))
if (!is.null(sample_ids)) {
  roh_files <- roh_files[sub("_ROH.txt$", "", basename(roh_files)) %in% sample_ids]
}
not_chrom_level_roh <- roh_files[!sapply(sub("_ROH.txt$", "", basename(roh_files)), function(s) is_chrom_level_species(species_for_run(s)))]
if (length(not_chrom_level_roh) > 0) {
  cat("Excluding (not chromosome-level assembly):", paste(sub("_ROH.txt$", "", basename(not_chrom_level_roh)), collapse = ", "), "\n")
}
roh_files <- roh_files[sapply(sub("_ROH.txt$", "", basename(roh_files)), function(s) is_chrom_level_species(species_for_run(s)))]
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
    tufte_par(mar = c(8, 5, 3, 2))
    barplot(top30$length_bp / 1e6, names.arg = top30$chrom, las = 2, cex.names = 0.7,
            col = MUTED_PALETTE[1], border = NA,
            ylab = "Total ROH length (Mb)", main = "")
    tufte_title(paste0("Top ", nrow(top30), " scaffolds by total ROH -- ",
                        length(unique(all_roh$srr)), " sample(s)"))
    dev.off()

    # per-sample total ROH -- labeled by species name, not Run ID
    by_sample <- aggregate(length_bp ~ species, data = all_roh, sum)
    by_sample <- by_sample[order(-by_sample$length_bp), ]
    png("plots/roh/total_roh_per_sample.png",
        width = max(900, nrow(by_sample) * 40), height = 700)
    tufte_par(mar = c(10, 5, 3, 2))
    barplot(by_sample$length_bp / 1e6, names.arg = by_sample$species, las = 2, cex.names = 0.7,
            col = MUTED_PALETTE[1], border = NA,
            ylab = "Total ROH length (Mb)", main = "")
    tufte_title(paste0("Total ROH per sample -- ", nrow(by_sample), " sample(s) completed"))
    dev.off()

    cat("Wrote plots/roh/top30_scaffolds_total_roh.png and plots/roh/total_roh_per_sample.png\n")

    # --- Chromosome painting + genome-wide FROH ---
    dir.create("plots/roh/painting", recursive = TRUE, showWarnings = FALSE)
    fai_for_species <- function(species) {
      file.path(ref_genome_path, species, paste0(species, ".fasta.fai"))
    }

    # NOTE: a "CM" GenBank-accession filter was tried here and reverted --
    # it incorrectly excluded species (e.g. Arabidopsis) whose reference
    # uses a different chromosome-naming convention. No single naming
    # pattern reliably separates real chromosomes from scaffolds across
    # every species' source here, so this paints whatever ROH was actually
    # called on.

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

      sample_roh <- all_roh[all_roh$srr == srr, ]

      # paint exactly the chromosome set present in this sample's own
      # ROH output, looked up in its .fai for real lengths
      called_chroms <- unique(sample_roh$chrom)
      top_scaffolds <- fai[fai$chrom %in% called_chroms, ]
      top_scaffolds <- top_scaffolds[order(-top_scaffolds$len), ]
      if (nrow(top_scaffolds) == 0) {
        warning(paste0("No .fai lengths matched the chromosomes ROH was called on for ", srr, " -- skipping painting"))
        next
      }

      png(file.path("plots/roh/painting", paste0(species, "_painting.png")),
          width = 1200, height = min(30000, max(400, nrow(top_scaffolds) * 25)))
      tufte_par(mar = c(4, 10, 3, 2))
      plot(NA, xlim = c(0, max(top_scaffolds$len)), ylim = c(0, nrow(top_scaffolds) + 1),
           yaxt = "n", xaxt = "n", xlab = "Position (bp)", ylab = "", main = "")
      axis(1, lwd = 0.6)
      axis(2, at = seq_len(nrow(top_scaffolds)), labels = rev(top_scaffolds$chrom),
           las = 2, cex.axis = 0.6, lwd = 0)
      tufte_title(paste0(species, " -- ROH painting (", nrow(top_scaffolds), " chromosomes)"))
      for (i in seq_len(nrow(top_scaffolds))) {
        y <- nrow(top_scaffolds) - i + 1
        chrom <- top_scaffolds$chrom[i]
        rect(0, y - 0.3, top_scaffolds$len[i], y + 0.3, col = "grey88", border = NA)
        segs <- sample_roh[sample_roh$chrom == chrom, ]
        if (nrow(segs) > 0) {
          rect(segs$start, y - 0.3, segs$end, y + 0.3, col = IUCN_COLORS[["CR"]], border = NA)
        }
      }
      dev.off()

      # Two alternative ROH-segment inclusion criteria for FROH, compared
      # side by side below: a RELATIVE cutoff (segment >= 1% of the length
      # of the chromosome/scaffold it's on -- scales with assembly, so a
      # tiny unplaced scaffold's near-zero-length ROH calls mostly fail it
      # too) and an ABSOLUTE cutoff (segment >= 1 Mb, the standard
      # "long ROH" threshold used to flag recent/close inbreeding in the
      # conservation genomics literature, independent of scaffold size).
      chrom_len_lookup <- setNames(fai$len, fai$chrom)
      sample_roh$chrom_len <- chrom_len_lookup[sample_roh$chrom]
      qualifying_pct1chrom <- sample_roh[!is.na(sample_roh$chrom_len) &
                                            sample_roh$length_bp >= 0.01 * sample_roh$chrom_len, ]
      qualifying_1mb <- sample_roh[sample_roh$length_bp >= 1e6, ]

      froh_rows[[srr]] <- data.frame(
        species = species,
        genome_len_bp = sum(fai$len),
        roh_len_bp = sum(sample_roh$length_bp),
        roh_len_bp_1pctchrom = sum(qualifying_pct1chrom$length_bp),
        roh_len_bp_1mb = sum(qualifying_1mb$length_bp)
      )
    }

    if (length(froh_rows) > 0) {
      froh_df <- do.call(rbind, froh_rows)
      froh_df$pct_genome_in_roh <- 100 * froh_df$roh_len_bp / froh_df$genome_len_bp
      froh_df$pct_genome_in_roh_1pctchrom <- 100 * froh_df$roh_len_bp_1pctchrom / froh_df$genome_len_bp
      froh_df$pct_genome_in_roh_1mb <- 100 * froh_df$roh_len_bp_1mb / froh_df$genome_len_bp
      froh_df <- froh_df[order(-froh_df$pct_genome_in_roh), ]
      write.csv(froh_df, "plots/roh/froh_comparison.csv", row.names = FALSE)

      png("plots/roh/froh_percent_genome.png",
          width = max(900, nrow(froh_df) * 60), height = 700)
      tufte_par(mar = c(10, 5, 3, 2))
      barplot(froh_df$pct_genome_in_roh, names.arg = froh_df$species, las = 2, cex.names = 0.8,
              col = MUTED_PALETTE[1], border = NA,
              ylab = "% of genome in ROH (FROH)", main = "")
      tufte_title(paste0("Genome-wide FROH -- ", nrow(froh_df), " sample(s)"))
      dev.off()

      cat("Wrote", length(froh_rows), "chromosome-painting plots to plots/roh/painting/",
          "and plots/roh/froh_percent_genome.png\n")

      # --- FROH by IUCN status ---
      # Styled after the "Raw genetic diversity by threat status" panel:
      # top-mounted numeric axis, IQR box + median tick + jittered raw
      # points per category, colored directly by IUCN status (no separate
      # legend needed since the category IS the row label).
      if (!is.null(iucn_lookup_df)) {
        froh_df$iucn <- sapply(froh_df$species, iucn_for_species)
        missing_iucn <- unique(froh_df$species[is.na(froh_df$iucn)])
        if (length(missing_iucn) > 0) {
          cat("No IUCN status found for:", paste(missing_iucn, collapse = ", "), "\n")
        }
        froh_iucn <- froh_df[!is.na(froh_df$iucn) & froh_df$iucn %in% IUCN_ORDER, ]

        if (nrow(froh_iucn) > 0) {
          present_order <- IUCN_ORDER[IUCN_ORDER %in% froh_iucn$iucn]
          n <- length(present_order)
          set.seed(1) # stable jitter across reruns
          png("plots/roh/froh_by_iucn_status.png", width = 950, height = 550)
          tufte_par(mar = c(2, 5, 4, 2))
          xr <- range(froh_iucn$pct_genome_in_roh, na.rm = TRUE)
          xr <- c(max(0, xr[1] - 0.02 * diff(xr)), xr[2] + 0.05 * diff(xr))
          plot(NA, xlim = xr, ylim = c(0.3, n + 0.7), xaxt = "n", yaxt = "n",
               xlab = "", ylab = "", main = "")
          axis(3, at = pretty(xr), lwd = 0.6, cex.axis = 0.85)
          mtext("% of genome in ROH (FROH)", side = 3, line = 2.2, cex = 0.95, col = "grey20")
          for (i in seq_len(n)) {
            cat_code <- present_order[i]
            y <- n - i + 1
            v <- froh_iucn$pct_genome_in_roh[froh_iucn$iucn == cat_code]
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
          cat("Wrote plots/roh/froh_by_iucn_status.png (", nrow(froh_iucn), "of", nrow(froh_df), "samples had an IUCN match)\n")
        } else {
          cat("No samples matched an IUCN status -- skipping froh_by_iucn_status.png\n")
        }
      }

      # --- Compare the two FROH cutoff implementations ---
      # Tufte-style slopegraph: one dot per species per cutoff, connected by
      # a line, species named directly at both ends instead of a legend --
      # the two vertical dot columns plus their slopes ARE the comparison.
      cmp <- froh_df[order(-froh_df$pct_genome_in_roh_1mb), ]
      x1 <- 1; x2 <- 2
      yr_cmp <- range(c(cmp$pct_genome_in_roh_1pctchrom, cmp$pct_genome_in_roh_1mb), na.rm = TRUE)
      label_ok <- nrow(cmp) <= 30
      png("plots/roh/froh_cutoff_comparison.png",
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
      cat("Wrote plots/roh/froh_cutoff_comparison.png and plots/roh/froh_comparison.csv\n")
    }
  }
}
