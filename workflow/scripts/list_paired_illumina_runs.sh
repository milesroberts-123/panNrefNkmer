#!/usr/bin/env bash
# For every Run/BioSample ID in a sample sheet, ask NCBI which constituent SRA
# runs are paired-end Illumina WGS (Platform == ILLUMINA, LibraryLayout ==
# PAIRED, LibraryStrategy == WGS). jules_fastq_dump applies only the first two
# (and only to BioSample IDs); this applies all three to every ID, including
# direct SRR/ERR/DRR entries, which the pipeline does not platform-check, and
# so also shows whether any non-WGS (e.g. RNA-seq) runs slipped in.
#
# Usage: list_paired_illumina_runs.sh <samples.tsv> <out.tsv>
#   samples.tsv: tab-separated, first column = Run/BioSample ID, header row.
#   out.tsv: one row per run: ID, Run, Platform, LibraryLayout, LibraryStrategy,
#            LibrarySource, PASS/FAIL
#            (ID with no SRA records -> "no_records").
# Needs internet; throttled to stay under NCBI's 3 requests/second limit.
# Uses the same naive comma split as the pipeline, so a quoted comma inside a
# runinfo field could in rare cases misread a row -- spot-check FAIL rows.

set -euo pipefail
sheet="$1"
out="$2"
printf "ID\tRun\tPlatform\tLibraryLayout\tLibraryStrategy\tLibrarySource\tfilter\n" > "$out"

tail -n +2 "$sheet" | cut -f1 | sort -u | while read -r id; do
    uids=$(curl -s --max-time 120 \
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=${id}&retmax=500" \
        | grep -oE "<Id>[0-9]+</Id>" | sed 's/<[^>]*>//g' | paste -sd, - || true)
    sleep 0.4
    if [[ -z "$uids" ]]; then
        printf "%s\tNA\tNA\tNA\tNA\tNA\tno_records\n" "$id" >> "$out"
        continue
    fi
    curl -s --max-time 120 \
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id=${uids}&rettype=runinfo&retmode=text" \
        | awk -F',' -v id="$id" '
            NR==1 {for(i=1;i<=NF;i++){if($i=="Run")r=i; if($i=="Platform")p=i; if($i=="LibraryLayout")l=i; if($i=="LibraryStrategy")s=i; if($i=="LibrarySource")o=i}; next}
            $r ~ /^[SED]RR/ && (id !~ /^[SED]RR/ || $r==id) {
                ok = ($p=="ILLUMINA" && $l=="PAIRED" && $s=="WGS") ? "PASS" : "FAIL"
                printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n", id, $r, $p, $l, $s, $o, ok
            }' >> "$out"
    sleep 0.4
done
echo "wrote $out"
