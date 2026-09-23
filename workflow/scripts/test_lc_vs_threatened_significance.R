#!/usr/bin/env Rscript
# Wilcoxon rank-sum test (Mann-Whitney U) for LC (non-threatened) vs
# everything else (threatened), across every FROH cutoff and every Ne
# harmonic-mean time period/method. Non-parametric on purpose -- FROH is a
# bounded percentage with a lot of near-zero values, and Ne spans many
# orders of magnitude, so neither is well-behaved enough to trust a t-test's
# normality assumption. Same test family as the "*" significance marker on
# the Wilder et al. reference panel (threatened vs non-threatened Ne).
#
# Running ~15 tests at once (5 FROH cutoffs + 5 periods x 2 methods) means a
# naive p<0.05 cutoff will produce false positives by chance alone -- this
# reports both the raw p-value and a Benjamini-Hochberg FDR-adjusted one
# (p.adjust(method="BH")), and you should read the adjusted column, not the
# raw one, when deciding what's actually worth trusting.
#
# Reads the small CSVs that plot_msmc_roh.R and analyze_ne_harmonic_means.R
# already write, so this is fast.
#
# Usage:
#   Rscript test_lc_vs_threatened_significance.R [froh_csv] [ne_summary_dir] [iucn_csv] [out_csv]
#   Rscript test_lc_vs_threatened_significance.R plots/roh/froh_comparison.csv plots/ne_summary /global/scratch/users/julesperez/snakemake_pipeline/metadata/complete_2026.csv plots/lc_vs_threatened_significance.csv

args <- commandArgs(trailingOnly = TRUE)
froh_csv <- if (length(args) >= 1 && nzchar(args[1])) args[1] else "plots/roh/froh_comparison.csv"
ne_summary_dir <- if (length(args) >= 2 && nzchar(args[2])) args[2] else "plots/ne_summary"
iucn_csv <- if (length(args) >= 3 && nzchar(args[3])) args[3] else stop("iucn_csv is required (scientificName,redlistCategory columns)")
out_csv <- if (length(args) >= 4 && nzchar(args[4])) args[4] else "plots/lc_vs_threatened_significance.csv"

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
threat_group_of <- function(iucn) {
  ifelse(iucn == "LC", "Non-threatened", ifelse(is.na(iucn), NA, "Threatened"))
}

run_test <- function(values, group, label) {
  keep <- !is.na(group) & is.finite(values)
  v <- values[keep]; g <- group[keep]
  n_nt <- sum(g == "Non-threatened")
  n_th <- sum(g == "Threatened")
  if (n_nt < 2 || n_th < 2) {
    return(data.frame(comparison = label, n_non_threatened = n_nt, n_threatened = n_th,
                       median_non_threatened = NA, median_threatened = NA,
                       p_value = NA, note = "too few samples in one group"))
  }
  med_nt <- median(v[g == "Non-threatened"])
  med_th <- median(v[g == "Threatened"])
  test <- suppressWarnings(wilcox.test(v[g == "Non-threatened"], v[g == "Threatened"]))
  data.frame(comparison = label, n_non_threatened = n_nt, n_threatened = n_th,
             median_non_threatened = med_nt, median_threatened = med_th,
             p_value = test$p.value, note = "")
}

results <- list()

# --- FROH cutoffs ---
if (file.exists(froh_csv)) {
  froh_df <- read.csv(froh_csv, stringsAsFactors = FALSE)
  froh_df$iucn <- sapply(froh_df$species, iucn_for_species)
  froh_df$threat_group <- threat_group_of(froh_df$iucn)

  froh_cols <- list(
    pct_genome_in_roh = "FROH, no length cutoff",
    pct_genome_in_roh_1pctchrom = "FROH, ROH >= 1% of chromosome",
    pct_genome_in_roh_1mb = "FROH, ROH >= 1 Mb",
    pct_genome_in_roh_2_5mb = "FROH, ROH >= 2.5 Mb",
    pct_genome_in_roh_1pctgenome = "FROH, ROH >= 1% of genome"
  )
  for (col in names(froh_cols)) {
    if (col %in% names(froh_df)) {
      results[[froh_cols[[col]]]] <- run_test(froh_df[[col]], froh_df$threat_group, froh_cols[[col]])
    }
  }
} else {
  cat("No froh_comparison.csv found at", froh_csv, "-- skipping FROH tests\n")
}

# --- Ne harmonic means, by method x period ---
for (method_name in c("MSMC2", "PSMC")) {
  csv_path <- file.path(ne_summary_dir, paste0("ne_harmonic_means_", tolower(method_name), ".csv"))
  if (!file.exists(csv_path)) {
    cat("No", csv_path, "-- skipping", method_name, "Ne tests\n")
    next
  }
  tbl <- read.csv(csv_path, stringsAsFactors = FALSE)
  tbl$iucn <- sapply(tbl$species, iucn_for_species)
  tbl$threat_group <- threat_group_of(tbl$iucn)
  for (p in unique(tbl$period)) {
    sub <- tbl[tbl$period == p, ]
    label <- paste0(method_name, " Ne (harmonic mean), ", p)
    results[[label]] <- run_test(sub$harmonic_mean_ne, sub$threat_group, label)
  }
}

if (length(results) == 0) {
  stop("No comparisons could be run -- check that froh_csv/ne_summary_dir point at real output")
}

out <- do.call(rbind, results)
rownames(out) <- NULL
out$p_adj_BH <- NA_real_
has_p <- !is.na(out$p_value)
out$p_adj_BH[has_p] <- p.adjust(out$p_value[has_p], method = "BH")
out$significant_BH_0.05 <- ifelse(has_p, out$p_adj_BH < 0.05, NA)

write.csv(out, out_csv, row.names = FALSE)
cat("\n")
print(out[, c("comparison", "n_non_threatened", "n_threatened", "p_value", "p_adj_BH", "significant_BH_0.05")], row.names = FALSE)
cat("\nWrote", out_csv, "\n")
