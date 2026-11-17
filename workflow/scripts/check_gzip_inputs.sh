#!/usr/bin/env bash
#
# Cheap structural validation for gzipped FASTQ inputs.
#
# Catches the failure mode where an upstream rule declares a file as an output
# but the tool wrote nothing useful (0 bytes, or a non-gzip file). Such a file
# satisfies Snakemake's existence check and only breaks several rules later,
# inside KMC, as "Error: Unexpected end of gzip file".
#
# Per file this reads 2 bytes, so it is safe to run on hundreds of GB of
# inputs. Pass --deep to additionally run `pigz -t`, which verifies the full
# gzip stream and CRCs but has to decompress the whole file.
#
# Usage:
#   check_gzip_inputs.sh FILE [FILE ...]
#   check_gzip_inputs.sh --deep FILE [FILE ...]
#
# Exit status: 0 if every file is well formed, 1 otherwise. Invalid files are
# reported on stderr, one line each.

set -uo pipefail

deep=0
if [ "${1:-}" = "--deep" ]; then
    deep=1
    shift
fi

if [ "$#" -eq 0 ]; then
    echo "check_gzip_inputs.sh: no input files given" >&2
    exit 1
fi

fail=0
for f in "$@"; do
    if [ ! -e "$f" ]; then
        printf 'INVALID %s: file does not exist\n' "$f" >&2
        fail=1
        continue
    fi
    if [ ! -s "$f" ]; then
        printf 'INVALID %s: zero bytes (no gzip header written)\n' "$f" >&2
        fail=1
        continue
    fi

    # gzip streams start with 1f 8b (RFC 1952)
    magic=$(od -An -tx1 -N2 "$f" | tr -d ' \n')
    if [ "$magic" != "1f8b" ]; then
        printf 'INVALID %s: not a gzip file (first bytes: %s)\n' "$f" "$magic" >&2
        fail=1
        continue
    fi

    if [ "$deep" -eq 1 ]; then
        # pigz where available, plain gzip otherwise; a missing tool must not
        # be reported as a corrupt file
        if command -v pigz >/dev/null 2>&1; then
            if ! pigz -t -p 4 "$f" 2>/dev/null; then
                printf 'INVALID %s: gzip integrity check failed (truncated or corrupt)\n' "$f" >&2
                fail=1
            fi
        elif command -v gzip >/dev/null 2>&1; then
            if ! gzip -t "$f" 2>/dev/null; then
                printf 'INVALID %s: gzip integrity check failed (truncated or corrupt)\n' "$f" >&2
                fail=1
            fi
        else
            printf 'check_gzip_inputs.sh: --deep needs pigz or gzip on PATH\n' >&2
            exit 2
        fi
    fi
done

exit "$fail"
