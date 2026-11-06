#!/usr/bin/env bash
#SBATCH --job-name=collect_for_paper
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=co_moilab
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

mkdir -p "$OUT"/{roh_rg,fai,ne/psmc,ne/msmc2,coverage,picard,pi,summary_tables,metadata,audit}
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
printf 'Run\tBioSample\tSpecies\troh_mtime\tn_RG_segments\tG\tAF_dflt\tab_mask_in_filtered_vcf\tmean_depth\n' > "$OUT/roh_provenance.tsv"
: > "$OUT/missing_roh.txt"
tail -n +2 "$SHEET" | while IFS=$'\t' read -r run bios sp rest; do
    f="results/roh/${run}_ROH.txt"
    if [[ ! -s "$f" ]]; then echo "$run" >> "$OUT/missing_roh.txt"; continue; fi
    dest="$OUT/roh_rg/${run}_ROH_RG.txt"
    grep -E $'^(#|RG\t)' "$f" > "$dest"
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
tail -n +2 "$SHEET" | cut -f3 | sort -u | while read -r sp; do
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

# ---- pi: newest pi_summary.txt per species, then species means ----
for R in $POSTROT_ROOTS; do
    d="$R/gatk_results/low_coverage_angsd"
    [[ -d "$d" ]] && find "$d" -maxdepth 2 -name pi_summary.txt -printf '%T@\t%TY-%Tm-%Td %TH:%TM\t%p\n'
done | sort -k1,1nr > "$OUT/pi/.candidates.tsv"
declare -A seen
printf 'Species\tchosen_file\tmtime\tolder_copies_ignored\n' > "$OUT/pi/PI_SOURCES.tsv"
while IFS=$'\t' read -r t mt p; do
    dn=$(basename "$(dirname "$p")")
    sp=${dn#angsd_pi_level_final_}; sp=${sp%%_mincov*}
    if [[ -n "${seen[$sp]:-}" ]]; then seen[$sp]=$((seen[$sp] + 1)); continue; fi
    seen[$sp]=0
    cp -p "$p" "$OUT/pi/${sp}_pi_summary.txt"
    printf '%s\t%s\t%s\n' "$sp" "$p" "$mt" >> "$OUT/pi/.chosen.tsv"
done < "$OUT/pi/.candidates.tsv"
if [[ -s "$OUT/pi/.chosen.tsv" ]]; then
    while IFS=$'\t' read -r sp p mt; do
        printf '%s\t%s\t%s\t%s\n' "$sp" "$p" "$mt" "${seen[$sp]}" >> "$OUT/pi/PI_SOURCES.tsv"
    done < "$OUT/pi/.chosen.tsv"
fi
for f in "$OUT"/pi/*_pi_summary.txt; do
    [[ -s "$f" ]] || continue
    tail -n +2 "$f" | awk -F, '!k[$1]++'
done | awk -F, '$5 != "" && $5+0 == $5 && tolower($5) !~ /nan|inf/ {s=$2; n[s]++; sum[s]+=$5; sq[s]+=$5*$5}
    END {print "Species,Replicates,Pi_Average,Pi_StdDev";
         for (s in n) {m=sum[s]/n[s]; v=(n[s]>1)?(sq[s]-n[s]*m*m)/(n[s]-1):0; if (v<0) v=0;
                       printf "%s,%d,%.9g,%.9g\n", s, n[s], m, sqrt(v)}}' > "$OUT/pi/pi_species_means_fresh.csv"
rm -f "$OUT/pi/.candidates.tsv" "$OUT/pi/.chosen.tsv"

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
    echo "pi species:          $(($(wc -l < "$OUT/pi/pi_species_means_fresh.csv") - 1))"
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
echo "Done: $OUT"
