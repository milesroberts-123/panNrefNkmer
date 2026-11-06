#!/usr/bin/env bash
#SBATCH --job-name=collect_for_paper
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=collect_for_paper_%j.out

# Gather everything the paper analyses need from the Savio runs into one
# dated folder, then copy it to Google Drive.
#
#   sbatch scripts/collect_for_paper.sh              # preferred: the ROH files can be 10-60 GB each
#   bash   scripts/collect_for_paper.sh --no-upload  # skip the Drive copy
#
# Run from workflow/. Safe to rerun: every run writes a new dated folder.
# Source mtimes are preserved (cp -p) and listed in MANIFEST.tsv so freshness can be checked.

set -u

WF=/global/scratch/projects/fc_moilab/julesperez/snakemake_pipeline/panNrefNkmer/workflow
PROJ=/global/scratch/projects/fc_moilab/julesperez/snakemake_pipeline
REFROOT=/global/scratch/projects/fc_moilab/julesperez/post_rot/new_refgenomes
POSTROT_ROOTS="/global/scratch/projects/fc_moilab/julesperez/post_rot /global/scratch/users/julesperez/post_rot"
STAMP=$(date +%Y-%m-%d_%H%M)
OUT=/global/scratch/projects/fc_moilab/julesperez/paper_pull/$STAMP
DRIVE="gdrive_moilab:PROJECTS/redlist_and_genomics/paper_pull/$STAMP"
UPLOAD=1
[[ "${1:-}" == "--no-upload" ]] && UPLOAD=0

cd "$WF" || { echo "cannot cd to $WF" >&2; exit 1; }

SHEET=""
for c in ../config/samples_medium_chromosome.tsv ../samples_medium_chromosome.tsv; do
    [[ -s "$c" ]] && { SHEET="$c"; break; }
done
[[ -n "$SHEET" ]] || { echo "sample sheet not found" >&2; exit 1; }

# Read Run / BioSample / Species by header name, not by position.
col() { head -1 "$SHEET" | tr -d '\r' | tr '\t' '\n' | grep -n -x "$1" | head -1 | cut -d: -f1; }
C_RUN=$(col Run); C_BIO=$(col BioSample); C_SP=$(col Species)
if [[ -z "$C_RUN" || -z "$C_SP" ]]; then
    echo "sample sheet $SHEET needs 'Run' and 'Species' columns; header is:" >&2
    head -1 "$SHEET" >&2
    exit 1
fi
sheet_rows() { awk -F'\t' -v r="$C_RUN" -v b="${C_BIO:-0}" -v s="$C_SP" 'NR>1 && $r != "" {gsub(/\r/,""); print $r "\t" (b>0 ? $b : "NA") "\t" $s}' "$SHEET"; }

mkdir -p "$OUT"/{roh_rg,fai,ne/psmc,ne/msmc2,coverage,picard,pi,het,summary_tables,metadata,audit}
echo "Writing to $OUT"

# ---- provenance: when, where, and whether the pipeline is still changing ----
{
    echo "collected:        $(date -Is)"
    echo "host:             $(hostname)"
    echo "sample sheet:     $SHEET ($(($(wc -l < "$SHEET") - 1)) samples)"
    latest_log=$(ls -t .snakemake/log/*.log 2>/dev/null | head -1)
    echo "newest driver log: ${latest_log:-none} ($(stat -c '%y' "$latest_log" 2>/dev/null | cut -d. -f1))"
    echo "slurm jobs queued or running: $(squeue -u "$USER" -h 2>/dev/null | wc -l)"
    echo "pipeline repo commit: $(git -C .. rev-parse --short HEAD 2>/dev/null)"
} > "$OUT/PROVENANCE.txt"
n_jobs=$(squeue -u "$USER" -h 2>/dev/null | wc -l)
if [[ "$n_jobs" -gt 0 ]]; then
    echo "WARNING: $n_jobs SLURM jobs still queued/running -- some samples will be missing; rerun this script when the pipeline finishes." | tee -a "$OUT/PROVENANCE.txt"
fi

# ---- ROH: RG (segment) lines only, plus a per-sample provenance table ----
# ST lines (per-site) are dropped: they make these files huge and are not needed.
# Reuse ROH extractions from ANY earlier collector folder: the extraction reads each 10-60 GB file in
# full. A copy is reused only if it is newer than its source file, and each folder's newest file is
# skipped because that job may have been cut off while writing it.
PREVS=(); PREV_NEWEST=()
while IFS= read -r d; do
    d=${d%/}
    [[ "$d" == "$OUT" || ! -d "$d/roh_rg" ]] && continue
    PREVS+=("$d")
    PREV_NEWEST+=("$(ls -t "$d/roh_rg" 2>/dev/null | head -1)")
done < <(ls -dt /global/scratch/projects/fc_moilab/julesperez/paper_pull/*/ 2>/dev/null)
echo "Reusing finished ROH extractions from ${#PREVS[@]} earlier folder(s) where possible"
printf 'Run\tBioSample\tSpecies\troh_mtime\tn_RG_segments\tG\tAF_dflt\tab_mask_in_filtered_vcf\tmean_depth\n' > "$OUT/roh_provenance.tsv"
: > "$OUT/missing_roh.txt"
sheet_rows | while IFS=$'\t' read -r run bios sp; do
    f="results/roh/${run}_ROH.txt"
    if [[ ! -s "$f" ]]; then echo "$run" >> "$OUT/missing_roh.txt"; continue; fi
    dest="$OUT/roh_rg/${run}_ROH_RG.txt"
    reused=0
    for i in ${!PREVS[@]}; do
        prev="${PREVS[$i]}/roh_rg/${run}_ROH_RG.txt"
        if [[ -s "$prev" && "$prev" -nt "$f" && "$(basename "$prev")" != "${PREV_NEWEST[$i]}" ]]; then
            cp -p "$prev" "$dest"; reused=1; break
        fi
    done
    if [[ "$reused" -eq 0 ]]; then
        grep -E $'^(#|RG\t)' "$f" > "$dest.tmp" && mv "$dest.tmp" "$dest"
    fi
    cl=$(grep -m1 -i 'command line' "$dest")
    g=$(echo "$cl" | sed -nE 's/.*-G ?([0-9.]+).*/\1/p'); g=${g:-NA}
    af=$(echo "$cl" | sed -nE 's/.*--AF-dflt ?([0-9.]+).*/\1/p'); af=${af:-NA}
    nrg=$(grep -c $'^RG\t' "$dest")
    vcf="results/bcfvcfs/${run}_filtered.vcf.gz"
    mask=NA
    if [[ -s "$vcf" ]]; then
        hdr=$(gzip -dc "$vcf" 2>/dev/null | awk '/^#CHROM/{exit} /^##bcftools_filterCommand/{print}')
        if echo "$hdr" | grep -q 'FMT/AD'; then mask=yes; else mask=no; fi
    fi
    depth=$(awk '{print $3}' "results/coverages/${run}.50k.coverage.txt" 2>/dev/null); depth=${depth:-NA}
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$run" "$bios" "$sp" "$(stat -c '%y' "$f" | cut -d. -f1)" "$nrg" "$g" "$af" "$mask" "$depth" >> "$OUT/roh_provenance.tsv"
done

# ---- reference .fai for chromosome lengths ----
: > "$OUT/missing_fai.txt"
sheet_rows | cut -f3 | sort -u | while read -r sp; do
    [[ -z "$sp" ]] && continue
    if [[ -s "$REFROOT/$sp/$sp.fasta.fai" ]]; then cp -p "$REFROOT/$sp/$sp.fasta.fai" "$OUT/fai/"; else echo "$sp" >> "$OUT/missing_fai.txt"; fi
done

# ---- mean depth per sample, and the sample sheet itself ----
printf 'Run\tmean_depth\n' > "$OUT/coverage/coverage_per_sample.tsv"
for f in results/coverages/*.50k.coverage.txt; do
    [[ -s "$f" ]] || continue
    printf '%s\t%s\n' "$(basename "$f" .50k.coverage.txt)" "$(awk '{print $3}' "$f")" >> "$OUT/coverage/coverage_per_sample.tsv"
done
cp -p "$SHEET" "$OUT/metadata/"

# ---- historical Ne: raw PSMC and MSMC2 curves (small; mu and generation time are applied later) ----
find results/psmc -maxdepth 1 -name '*.psmc' -exec cp -p {} "$OUT/ne/psmc/" \; 2>/dev/null
for d in results/msmc2/*/; do
    id=$(basename "$d")
    [[ -s "$d/msmc2.final.txt" ]] && cp -p "$d/msmc2.final.txt" "$OUT/ne/msmc2/${id}.msmc2.final.txt"
done
find results/picard_metrics -maxdepth 1 -name '*.txt' -exec cp -p {} "$OUT/picard/" \; 2>/dev/null

# ---- IUCN, generation times, audit tables ----
for f in "$PROJ"/metadata/complete_2026.csv "$PROJ"/metadata/gen_times.csv; do [[ -s "$f" ]] && cp -p "$f" "$OUT/metadata/"; done
for f in "$PROJ"/audit/*.tsv "$PROJ"/audit/*.txt; do [[ -s "$f" ]] && cp -p "$f" "$OUT/audit/"; done

# ---- pi and heterozygosity: SPECIES_PROCESSING_GOAT is the newest source and takes priority;
# ---- the older gatk_results outputs only fill species (pi) or samples (heterozygosity) GOAT lacks ----
sp_from_path() {
    local p="$1" c
    local IFS=/
    for c in $p; do
        case "$c" in
            angsd_heterozygosity) continue ;;
            angsd_pi_level_final_*) c=${c#angsd_pi_level_final_}; echo "${c%%_mincov*}"; return ;;
            *_aligned_to_*) echo "${c%%_aligned_to_*}"; return ;;
            *_heterozygosity) echo "${c%_heterozygosity}"; return ;;
        esac
    done
    echo UNKNOWN
}

# $1 = file name, $2 = gatk_results subfolder, $3 = find depth inside it.
# Output: priority (0 = GOAT, 1 = GATK), epoch, date, path -- GOAT first, newest first within each.
find_sources() {
    local R
    for R in $POSTROT_ROOTS; do
        [[ -d "$R/SPECIES_PROCESSING_GOAT" ]] && find "$R/SPECIES_PROCESSING_GOAT" -maxdepth 6 -name "$1" -type f \
            -printf '0\t%T@\t%TY-%Tm-%Td %TH:%TM\t%p\n' 2>/dev/null
        [[ -d "$R/gatk_results/$2" ]] && find "$R/gatk_results/$2" -maxdepth "$3" -name "$1" -type f \
            -printf '1\t%T@\t%TY-%Tm-%Td %TH:%TM\t%p\n' 2>/dev/null
    done | sort -t$'\t' -k1,1n -k2,2nr
}

# pi: one file per species, from GOAT if GOAT has it, otherwise from gatk_results
declare -A seen_pi
printf 'Species\tsource\tmtime\tchosen_file\n' > "$OUT/pi/PI_SOURCES.tsv"
while IFS=$'\t' read -r prio ep mt p; do
    sp=$(sp_from_path "$p")
    [[ -n "${seen_pi[$sp]:-}" ]] && continue
    seen_pi[$sp]=1
    src=GOAT; [[ "$prio" == 1 ]] && src=GATK
    cp -p "$p" "$OUT/pi/${sp}_pi_summary.txt"
    printf '%s\t%s\t%s\t%s\n' "$sp" "$src" "$mt" "$p" >> "$OUT/pi/PI_SOURCES.tsv"
done < <(find_sources pi_summary.txt low_coverage_angsd 2)
for f in "$OUT"/pi/*_pi_summary.txt; do
    [[ -s "$f" ]] || continue
    tail -n +2 "$f" | awk -F, '!k[$1]++'
done | awk -F, '$5 != "" && $5+0 == $5 && tolower($5) !~ /nan|inf/ {s=$2; n[s]++; sum[s]+=$5; sq[s]+=$5*$5}
    END {print "Species,Replicates,Pi_Average,Pi_StdDev";
         for (s in n) {m=sum[s]/n[s]; v=(n[s]>1)?(sq[s]-n[s]*m*m)/(n[s]-1):0; if (v<0) v=0;
                       printf "%s,%d,%.9g,%.9g\n", s, n[s], m, sqrt(v)}}' > "$OUT/pi/pi_species_means_fresh.csv"

# heterozygosity: per-sample summaries. Every GOAT row is kept; a gatk_results row is added only
# when GOAT has no row for that species + sample.
mkdir -p "$OUT/het/goat" "$OUT/het/gatk"
declare -A seen_het
printf 'Species\tsource\tmtime\tchosen_file\n' > "$OUT/het/HET_SOURCES.tsv"
while IFS=$'\t' read -r prio ep mt p; do
    sp=$(sp_from_path "$p")
    src=GOAT; d=goat; [[ "$prio" == 1 ]] && { src=GATK; d=gatk; }
    [[ -n "${seen_het[$sp,$src]:-}" ]] && continue
    seen_het[$sp,$src]=1
    cp -p "$p" "$OUT/het/$d/${sp}_heterozygosity_summary.txt"
    printf '%s\t%s\t%s\t%s\n' "$sp" "$src" "$mt" "$p" >> "$OUT/het/HET_SOURCES.tsv"
done < <(find_sources heterozygosity_summary.txt angsd_heterozygosity 3)
goat_files=(); gatk_files=()
for f in "$OUT"/het/goat/*_heterozygosity_summary.txt; do [[ -e "$f" ]] && goat_files+=("$f"); done
for f in "$OUT"/het/gatk/*_heterozygosity_summary.txt; do [[ -e "$f" ]] && gatk_files+=("$f"); done
if (( ${#goat_files[@]} + ${#gatk_files[@]} > 0 )); then
    awk -F, -v OFS=, -v hdr="Species,Sample,AvgCoverage,MinCovThreshold,Heterozygosity" -v goatdir="$OUT/het/goat/" '
        BEGIN { print hdr ",source" }
        FNR==1 { ok = ($0 == hdr); if (!ok) print "SKIPPED (unexpected header): " FILENAME > "/dev/stderr"
                 goat = (index(FILENAME, goatdir) == 1); next }
        !ok { next }
        { key = $1 SUBSEP $2
          if (goat) have[key] = 1
          else if (key in have) next
          k4 = key SUBSEP $4 SUBSEP goat
          if (k4 in dup) next
          dup[k4] = 1
          print $0, (goat ? "GOAT" : "GATK") }' \
        ${goat_files[@]+"${goat_files[@]}"} ${gatk_files[@]+"${gatk_files[@]}"} > "$OUT/het/het_per_sample_merged.csv"
fi

# raw per-sample heterozygosity SFS (one line each), bundled into one table; GOAT first, GATK fills what GOAT lacks
printf 'Species\tSample\tMinCov\tsource\tsfs\n' > "$OUT/het/het_sfs_all.tsv"
declare -A have_sfs
while IFS=$'\t' read -r prio ep mt p; do
    sp=$(sp_from_path "$p")
    b=$(basename "$p" .sfs)
    mc=${b##*_mincov}; sample=${b%_mincov*}
    key="$sp|$sample|$mc"
    [[ -n "${have_sfs[$key]:-}" ]] && continue
    have_sfs[$key]=1
    src=GOAT; [[ "$prio" == 1 ]] && src=GATK
    vals=$(head -n 1 "$p" | tr -s ' \t' ',' | sed 's/^,//; s/,$//')
    printf '%s\t%s\t%s\t%s\t%s\n' "$sp" "$sample" "$mc" "$src" "$vals" >> "$OUT/het/het_sfs_all.tsv"
done < <(find_sources '*_mincov*.sfs' angsd_heterozygosity 4)

# raw per-replicate pi SFS (folded), from the same source that supplied each species' pi_summary
printf 'Species\tReplicate\tsource\tfolded_sfs\n' > "$OUT/pi/pi_sfs_all.tsv"
tail -n +2 "$OUT/pi/PI_SOURCES.tsv" | while IFS=$'\t' read -r sp src mt chosen; do
    d=$(dirname "$chosen")
    for f in "$d"/n_*_rep*/*_n*_rep*.sfs; do
        [[ -s "$f" ]] || continue
        rep=$(basename "$(dirname "$f")")
        vals=$(head -n 1 "$f" | tr -s ' \t' ',' | sed 's/^,//; s/,$//')
        printf '%s\t%s\t%s\t%s\n' "$sp" "$rep" "$src" "$vals" >> "$OUT/pi/pi_sfs_all.tsv"
    done
done

# ---- consolidated pi / heterozygosity tables if they exist on the cluster (newest copy of each) ----
{
    for R in $POSTROT_ROOTS "$HOME"; do
        [[ -d "$R" ]] && find "$R" -maxdepth 3 -type f \( -name 'pi_averages*' -o -name 'het_per_sample*' -o -name 'consolidated_heterozygosity*' \) \
            -printf '%T@\t%TY-%Tm-%Td %TH:%TM\t%s\t%p\n' 2>/dev/null
    done
} | sort -k1,1nr > "$OUT/summary_tables/CANDIDATES.tsv"
for key in pi_averages het_per_sample consolidated_heterozygosity; do
    p=$(awk -F'\t' -v k="$key" 'index($4, "/" k) {print $4; exit}' "$OUT/summary_tables/CANDIDATES.tsv")
    [[ -n "$p" ]] && cp -p "$p" "$OUT/summary_tables/"
done

# ---- manifest ----
(cd "$OUT" && find . -type f ! -name MANIFEST.tsv -printf '%P\t%s\t%TY-%Tm-%Td %TH:%TM\n' | sort) > "$OUT/MANIFEST.tsv"

{
    echo
    echo "ROH segment files:   $(ls "$OUT/roh_rg" | wc -l)   (samples with no ROH output: $(wc -l < "$OUT/missing_roh.txt"))"
    echo "  without the 0.30/0.70 allele-balance mask in their filtered VCF: $(awk -F'\t' 'NR>1 && $8=="no"' "$OUT/roh_provenance.tsv" | wc -l)   (mask not checkable: $(awk -F'\t' 'NR>1 && $8=="NA"' "$OUT/roh_provenance.tsv" | wc -l))"
    echo "reference .fai:      $(ls "$OUT/fai" | wc -l)   (missing: $(wc -l < "$OUT/missing_fai.txt"))"
    echo "PSMC curves:         $(ls "$OUT/ne/psmc" | wc -l)"
    echo "MSMC2 curves:        $(ls "$OUT/ne/msmc2" | wc -l)"
    echo "pi species:          $(($(wc -l < "$OUT/pi/pi_species_means_fresh.csv") - 1))   (from GOAT: $(awk -F'\t' 'NR>1 && $2=="GOAT"' "$OUT/pi/PI_SOURCES.tsv" | wc -l), filled from gatk_results: $(awk -F'\t' 'NR>1 && $2=="GATK"' "$OUT/pi/PI_SOURCES.tsv" | wc -l))"
    echo "het SFS rows:        $(($(wc -l < "$OUT/het/het_sfs_all.tsv") - 1))   pi SFS rows: $(($(wc -l < "$OUT/pi/pi_sfs_all.tsv") - 1))"
    echo "heterozygosity rows: $([[ -s "$OUT/het/het_per_sample_merged.csv" ]] && echo $(($(wc -l < "$OUT/het/het_per_sample_merged.csv") - 1)) || echo 0)   (GOAT: $(awk -F, 'NR>1 && $NF=="GOAT"' "$OUT/het/het_per_sample_merged.csv" 2>/dev/null | wc -l), gatk_results fill: $(awk -F, 'NR>1 && $NF=="GATK"' "$OUT/het/het_per_sample_merged.csv" 2>/dev/null | wc -l))"
} | tee -a "$OUT/PROVENANCE.txt"
du -sh "$OUT"

# ---- copy to Google Drive ----
if [[ "$UPLOAD" -eq 1 ]]; then
    command -v module >/dev/null 2>&1 && module load rclone 2>/dev/null
    if command -v rclone >/dev/null 2>&1; then
        rclone copy "$OUT" "$DRIVE" --progress && echo "Uploaded to $DRIVE"
    else
        echo "rclone not available; upload manually:  rclone copy $OUT $DRIVE --progress"
    fi
fi
touch "$OUT/COMPLETE"
echo "Done: $OUT"
