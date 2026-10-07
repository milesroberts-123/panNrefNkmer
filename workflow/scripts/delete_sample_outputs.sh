#!/usr/bin/env bash
# Remove every pipeline output for the Run/BioSample IDs in a list, so the next
# Snakemake invocation rebuilds those samples from scratch (e.g. after a
# download filter change made their old reads invalid).
#
# Usage (run from the workflow/ directory, with no Snakemake driver running):
#   bash scripts/delete_sample_outputs.sh ids.txt            # dry run: only lists
#   bash scripts/delete_sample_outputs.sh ids.txt --delete   # actually removes
#
# ids.txt: one ID per line (blank lines ignored). Anything that is not purely
# letters/digits is skipped, so a stray glob character can never widen a match.
# Matching is delimiter-anchored ("ID.*", "ID_*", or exactly "ID") so one ID can
# never match another ID it is a prefix of. Refuses to touch anything outside
# ./results/.

set -euo pipefail
ids="${1:-}"
mode="${2:-}"
[[ -n "$ids" && -f "$ids" ]] || { echo "usage: $0 <ids.txt> [--delete]" >&2; exit 1; }
[[ -d results ]] || { echo "no ./results here -- run from the workflow/ directory" >&2; exit 1; }

shopt -s nullglob
list=$(mktemp)
trap 'rm -f "$list"' EXIT

while read -r id; do
    [[ -z "$id" ]] && continue
    if [[ ! "$id" =~ ^[A-Za-z0-9]+$ ]]; then
        echo "skipping odd ID: $id" >&2
        continue
    fi
    for p in results/*/"$id".* results/*/"$id"_* results/msmc2/"$id" \
             results/raw_reads/.staging_"$id" results/tmp_picard_"$id"; do
        if [[ -e "$p" ]]; then echo "$p"; fi
    done
done < "$ids" | sort -u > "$list"

n=$(wc -l < "$list")
echo "$n paths for $(grep -c . "$ids") IDs (first 40):"
head -40 "$list"

if grep -qv '^results/' "$list"; then
    echo "REFUSING: a path outside results/ is in the list" >&2
    exit 1
fi

if [[ "$mode" == "--delete" ]]; then
    tr '\n' '\0' < "$list" | xargs -0 rm -rf
    echo "deleted $n paths"
else
    echo "dry run only -- rerun with --delete to remove these"
fi
