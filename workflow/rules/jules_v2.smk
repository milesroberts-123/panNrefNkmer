# NOT temp() -- these live in the shared reference genome directory, not in
# results/. They're durable artifacts other scripts depend on, and the bwa
# index takes hours to rebuild.
rule jules_bwa_index:
    input:
        config["reference_genome_path"] + "{ref}/{ref}.fasta",
    output:
        amb=config["reference_genome_path"] + "{ref}/{ref}.fasta.amb",
        ann=config["reference_genome_path"] + "{ref}/{ref}.fasta.ann",
        bwt=config["reference_genome_path"] + "{ref}/{ref}.fasta.bwt",
        pac=config["reference_genome_path"] + "{ref}/{ref}.fasta.pac",
        sa=config["reference_genome_path"] + "{ref}/{ref}.fasta.sa",
    conda:
        "../envs/bwa.yaml"
    shell:
        """
        bwa index {input}
        """

rule jules_samtools_faidx:
    input:
        config["reference_genome_path"] + "{ref}/{ref}.fasta"
    output:
        config["reference_genome_path"] + "{ref}/{ref}.fasta.fai"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        samtools faidx {input}
        """

rule jules_fastq_dump:
    output:
        r1=temp("results/raw_reads/{ID}_1.fastq.gz"),
        r2=temp("results/raw_reads/{ID}_2.fastq.gz")
    # Caps how many of these run concurrently against the global ncbi_slots
    # pool (set in profiles/default/config.yaml), independent of the overall
    # `jobs: 10000` limit -- without this, a batch-wide rerun (e.g. triggered
    # by editing this rule's own code, which Snakemake's default `code`
    # rerun-trigger treats as invalidating every sample's existing download)
    # submits every sample's SRA download at once, flooding NCBI and failing
    # nearly the whole batch instead of a targeted few.
    resources:
        ncbi_slots=1
    conda:
        "../envs/sra.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p results/raw_reads

        if [[ "{wildcards.ID}" =~ ^SAM ]]; then
            # BioSample accession -- resolve to its constituent SRA runs,
            # download+dump each, concatenate into one R1/R2 pair.
            #
            # Uses NCBI's eutils REST API directly rather than edirect's
            # esearch/efetch: edirect's efetch shells out to an `xtract`
            # helper that isn't installed here, and instead of failing it
            # feeds its own error text forward as the query, producing a
            # garbage request that NCBI 400s. curl has no such hidden
            # dependency, and --max-time bounds it so a stalled request
            # can't silently burn the whole job walltime.
            echo "{wildcards.ID} is a BioSample -- resolving constituent SRA runs..."

            # Fetch the run metadata with retries. A rate-limited or failed
            # eutils reply (NCBI allows ~3 requests/s; a restart launches many
            # of these jobs at once) is NOT valid runinfo CSV, and used to be
            # indistinguishable from "this BioSample has no WGS runs" -- it
            # wrongly rejected SAMN43038876, which has two NovaSeq WGS runs.
            # A real runinfo reply starts with the header "Run,"; anything
            # else is retried with a jittered, growing delay. `|| true` keeps
            # a failed grep/curl from killing the script under pipefail.
            runinfo=""
            for attempt in 1 2 3 4 5; do
                sleep $(( (RANDOM % 10) + 1 ))
                uids=$(curl -s --max-time 120 \\
                    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term={wildcards.ID}&retmax=500" \\
                    | grep -oE "<Id>[0-9]+</Id>" | sed 's/<[^>]*>//g' | paste -sd, - || true)
                if [[ -n "$uids" ]]; then
                    sleep 1
                    runinfo=$(curl -s --max-time 120 \\
                        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id=${{uids}}&rettype=runinfo&retmode=text" || true)
                    if [[ "$runinfo" == Run,* ]]; then
                        break
                    fi
                fi
                runinfo=""
                echo "NCBI metadata for {wildcards.ID} unavailable or invalid (attempt $attempt of 5) -- retrying" >&2
                sleep $(( attempt * 10 ))
            done
            if [[ -z "$runinfo" ]]; then
                echo "Error: could not retrieve NCBI run metadata for BioSample {wildcards.ID} after 5 attempts (NCBI rate limiting or outage?) -- this says nothing about the data itself; rerun to retry" >&2
                exit 1
            fi

            # Some BioSamples mix platforms and/or layouts (e.g. an Illumina
            # short-read assembly plus an Oxford Nanopore long-read run
            # submitted under the same BioSample; or a paired-end population
            # run alongside an unrelated single-end pilot/control lane). This
            # pipeline is built for Illumina short-read PAIRED data -- bwa
            # mapping, mpileup, ROH, PSMC, MSMC2 all assume it -- so runs are
            # filtered on both the runinfo CSV's Platform and LibraryLayout
            # columns rather than relying on read-layout heuristics
            # downstream. Confirmed two distinct real cases: SAMN18024572
            # (Acanthochlamys bracteata) was entirely OXFORD_NANOPORE/
            # PromethION (caught by the Platform check); a broad sweep of
            # other BioSamples' staging dirs turned up single, unsuffixed
            # fastq.gz files (no _1/_2) even though Platform=="ILLUMINA" --
            # i.e. genuinely single-end Illumina runs mixed into otherwise
            # paired-end BioSamples (caught only by the LibraryLayout check).
            # LibraryStrategy == WGS is required too: the input criteria are
            # paired-end Illumina short-read WGS, and a BioSample can also
            # carry RNA-seq/amplicon/ATAC runs that are Illumina and paired
            # but would make ROH/PSMC/MSMC2 meaningless if concatenated in.
            runs=$(echo "$runinfo" \\
                | awk -F',' 'NR==1{{for(i=1;i<=NF;i++){{if($i=="Run")rcol=i; if($i=="Platform")pcol=i; if($i=="LibraryLayout")lcol=i; if($i=="LibraryStrategy")scol=i}} next}} $pcol=="ILLUMINA" && $lcol=="PAIRED" && $scol=="WGS"{{print $rcol}}' \\
                | grep -E '^[SED]RR' || true)

            # Say explicitly, in the job log, which runs were dropped and
            # which of the three criteria (platform/layout/strategy) each one
            # failed, so a rebuilt sample's composition can be verified.
            excluded=$(echo "$runinfo" \\
                | awk -F',' 'NR==1{{for(i=1;i<=NF;i++){{if($i=="Run")rcol=i; if($i=="Platform")pcol=i; if($i=="LibraryLayout")lcol=i; if($i=="LibraryStrategy")scol=i}} next}} $rcol ~ /^[SED]RR/ && !($pcol=="ILLUMINA" && $lcol=="PAIRED" && $scol=="WGS"){{print $rcol ":" $pcol "/" $lcol "/" $scol}}' || true)
            if [[ -n "$excluded" ]]; then
                echo "Excluded from {wildcards.ID} (not ILLUMINA/PAIRED/WGS) -- run:platform/layout/strategy: $(echo "$excluded" | paste -sd' ' -)"
            fi

            if [[ -z "$runs" ]]; then
                echo "Error: no paired-end Illumina WGS SRA runs found for BioSample {wildcards.ID} (other platforms/layouts/library strategies may exist but are filtered out)" >&2
                exit 1
            fi
            echo "Found paired-end Illumina WGS runs for {wildcards.ID}: ${{runs}}"

            # Stable, resumable staging dir on shared scratch (not mktemp'd
            # /tmp, which wouldn't survive a retry on a different node).
            # Only removed after full success; a failed attempt leaves it
            # for the next retry to resume from.
            staging_dir="results/raw_reads/.staging_{wildcards.ID}"
            mkdir -p "$staging_dir"

            # Bounded parallel download (cap 4, matches cpus_per_task).
            # Each run is idempotent: skips if already downloaded.
            download_one_run() {{
                local run="$1"
                local dir="$2"
                local r1="${{dir}}/${{run}}_1.fastq.gz"
                local r2="${{dir}}/${{run}}_2.fastq.gz"

                # A run already found unusable on a previous attempt (it
                # downloaded fine but never produced a valid R1/R2 pair) stays
                # skipped on retry -- don't re-download hundreds of GB to
                # reach the same conclusion.
                if [[ -e "${{dir}}/${{run}}.skip" ]]; then
                    echo "--- run ${{run}} previously marked unusable (no valid R1/R2 pair), skipping ---"
                    return 0
                fi

                if [[ -s "$r1" && -s "$r2" ]] && gzip -t "$r1" 2>/dev/null && gzip -t "$r2" 2>/dev/null; then
                    echo "--- run ${{run}} already downloaded (from a prior attempt), skipping ---"
                    return 0
                elif [[ -s "$r1" || -s "$r2" ]]; then
                    # leftover from a prior attempt that got cut short (e.g. a SLURM
                    # walltime kill mid-fastq-dump) -- non-empty but not a valid
                    # gzip stream, so `-s` alone would have wrongly treated it as
                    # complete. Clear it so this run re-downloads from scratch.
                    echo "--- run ${{run}}: found a truncated/corrupt file from a prior attempt, re-downloading ---"
                    rm -f "$r1" "$r2"
                fi

                # If a prior attempt at this specific run was killed mid-download
                # (walltime, node failure, etc.), prefetch leaves its own
                # ${{run}}/${{run}}.sra.lock file behind under $dir and refuses
                # to redownload -- every retry fails immediately with "lock
                # exists ... download canceled" even though nothing is still
                # running. Clear any leftover state for this run before
                # starting, independent of whether r1/r2 above needed clearing.
                if [[ -d "${{dir}}/${{run}}" ]]; then
                    echo "--- run ${{run}}: found a leftover ${{dir}}/${{run}} directory from a prior attempt -- removing before redownloading (likely a stale prefetch lock) ---"
                    rm -rf "${{dir}}/${{run}}"
                fi

                echo "--- downloading run ${{run}} (part of BioSample {wildcards.ID}) ---"
                # A failing prefetch/fastq-dump is a DOWNLOAD problem (network,
                # disk, killed job) -- fail the whole BioSample so a retry
                # picks it up, rather than silently shrinking the sample.
                if ! {{ prefetch --max-size 5000G -O "$dir" "$run" \\
                        && fastq-dump --gzip --clip --outdir "$dir" --split-3 --skip-technical "${{dir}}/${{run}}/${{run}}.sra"; }}; then
                    echo "Error: prefetch/fastq-dump failed for run ${{run}} (download problem, not skipped)" >&2
                    return 1
                fi

                # Tools succeeded but there is no valid pair: the run's
                # metadata says PAIRED but the data are not (e.g. unpaired
                # spots), so fastq-dump --split-3 wrote a single unsuffixed
                # file. That is a permanent property of the run, so skip just
                # this run (loudly) instead of failing the whole BioSample.
                if [[ ! -s "$r1" || ! -s "$r2" ]]; then
                    echo "WARNING: run ${{run}} downloaded but did not produce both non-empty R1/R2 files -- marking it unusable and skipping it; the rest of {wildcards.ID} continues" >&2
                    rm -f "$r1" "$r2" "${{dir}}/${{run}}.fastq.gz"
                    rm -rf "${{dir}}/${{run}}"
                    touch "${{dir}}/${{run}}.skip"
                    return 0
                fi
                if ! gzip -t "$r1" 2>/dev/null || ! gzip -t "$r2" 2>/dev/null; then
                    echo "Error: run ${{run}} produced a truncated/corrupt gzip file (likely killed mid-download, e.g. by a SLURM walltime limit) -- deleting partial output so a retry starts clean" >&2
                    rm -f "$r1" "$r2"
                    return 1
                fi
            }}
            export -f download_one_run

            echo "$runs" | xargs -P 4 -I{{}} bash -c 'download_one_run "$@"' _ {{}} "$staging_dir" \\
                || {{ echo "Error: one or more constituent SRA run downloads failed for BioSample {wildcards.ID} -- staging dir left in place at ${{staging_dir}} for the next retry to resume from" >&2; exit 1; }}

            # Build the merge list from ONLY the runs selected above -- never
            # glob the staging dir. It can hold stale per-run files from
            # earlier attempts made under a different filter (including runs
            # the current filter excludes, e.g. RNA-seq/Hi-C), which a glob
            # would silently merge into the sample or trip the corruption
            # check on. Runs marked .skip (downloaded but no valid R1/R2
            # pair) are dropped here, loudly.
            good_runs=()
            skipped_runs=()
            for run in $runs; do
                if [[ -e "${{staging_dir}}/${{run}}.skip" ]]; then
                    skipped_runs+=("$run")
                else
                    good_runs+=("$run")
                fi
            done
            if [[ ${{#skipped_runs[@]}} -gt 0 ]]; then
                echo "WARNING: {wildcards.ID}: skipped ${{#skipped_runs[@]}} unusable run(s) (no valid R1/R2 pair): ${{skipped_runs[*]}}"
            fi
            if [[ ${{#good_runs[@]}} -eq 0 ]]; then
                echo "Error: no usable paired-end WGS runs left for BioSample {wildcards.ID} after skipping unusable ones. Staging dir left in place for inspection." >&2
                exit 1
            fi

            # Re-validate every selected per-run file immediately before
            # merging, independent of download_one_run's own check. `cat` has
            # no way to tell a truncated gzip member from a complete one, so
            # this is the last point where a corrupt run can be named
            # individually, rather than surfacing later as an unexplained
            # coverage/read-count shortfall on the merged {wildcards.ID} file.
            r1_list=()
            r2_list=()
            for run in "${{good_runs[@]}}"; do
                for f in "${{staging_dir}}/${{run}}_1.fastq.gz" "${{staging_dir}}/${{run}}_2.fastq.gz"; do
                    if [[ ! -s "$f" ]] || ! gzip -t "$f" 2>/dev/null; then
                        echo "Error: ${{f}} is missing or a truncated/corrupt gzip file -- refusing to merge into {wildcards.ID}. Staging dir left in place for inspection." >&2
                        exit 1
                    fi
                done
                r1_list+=("${{staging_dir}}/${{run}}_1.fastq.gz")
                r2_list+=("${{staging_dir}}/${{run}}_2.fastq.gz")
            done

            cat "${{r1_list[@]}}" > results/raw_reads/{wildcards.ID}_1.fastq.gz
            cat "${{r2_list[@]}}" > results/raw_reads/{wildcards.ID}_2.fastq.gz
            echo "Merged ${{#good_runs[@]}} run(s) into {wildcards.ID}: ${{good_runs[*]}}"

            # Verify the combined output before deleting the only copies of
            # the source data. gzip -t on a `cat`-concatenated multi-member
            # gzip stream still validates fine (gzip supports concatenated
            # members), so this also catches a truncated per-run file that
            # somehow made it past download_one_run's own check.
            if [[ ! -s results/raw_reads/{wildcards.ID}_1.fastq.gz || ! -s results/raw_reads/{wildcards.ID}_2.fastq.gz ]]; then
                echo "Error: concatenation produced an empty/missing output for {wildcards.ID} -- staging dir left in place at ${{staging_dir}}" >&2
                exit 1
            fi
            if ! gzip -t results/raw_reads/{wildcards.ID}_1.fastq.gz 2>/dev/null || ! gzip -t results/raw_reads/{wildcards.ID}_2.fastq.gz 2>/dev/null; then
                echo "Error: concatenated output for {wildcards.ID} is a truncated/corrupt gzip stream -- staging dir left in place at ${{staging_dir}} for inspection" >&2
                rm -f results/raw_reads/{wildcards.ID}_1.fastq.gz results/raw_reads/{wildcards.ID}_2.fastq.gz
                exit 1
            fi
            rm -rf "$staging_dir"
        else
            # Plain Run accession (SRR/ERR/DRR) -- single download.
            #
            # Same input criteria as the BioSample branch (paired-end
            # Illumina short-read WGS), checked BEFORE downloading so a
            # non-conforming run is rejected up front with a clear message
            # instead of downloading hundreds of GB first. Previously this
            # branch relied on the R1/R2 check below, which catches
            # single-end/long-read data but not paired-end non-Illumina
            # (e.g. DNBSEQ) or non-WGS (e.g. RNA-seq) runs.
            # Fails closed: if NCBI metadata can't be retrieved after 3
            # tries the job exits rather than download unverified data (a
            # rerun just retries). The 3 tries are jittered so a batch of
            # jobs starting together doesn't hammer eutils in lockstep.
            meta=""
            for attempt in 1 2 3; do
                uid=$(curl -s --max-time 120 \\
                    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term={wildcards.ID}&retmax=500" \\
                    | grep -oE "<Id>[0-9]+</Id>" | sed 's/<[^>]*>//g' | paste -sd, - || true)
                if [[ -n "$uid" ]]; then
                    meta=$(curl -s --max-time 120 \\
                        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id=${{uid}}&rettype=runinfo&retmode=text" \\
                        | awk -F',' -v want="{wildcards.ID}" 'NR==1{{for(i=1;i<=NF;i++){{if($i=="Run")rcol=i; if($i=="Platform")pcol=i; if($i=="LibraryLayout")lcol=i; if($i=="LibraryStrategy")scol=i}} next}} $rcol==want {{print $pcol "/" $lcol "/" $scol; exit}}' || true)
                fi
                [[ -n "$meta" ]] && break
                sleep $(( (RANDOM % 5) + 3 ))
            done
            if [[ -z "$meta" ]]; then
                echo "Error: could not retrieve NCBI run metadata for {wildcards.ID} after 3 attempts -- refusing to download unverified data (rerun to retry)" >&2
                exit 1
            fi
            if [[ "$meta" != "ILLUMINA/PAIRED/WGS" ]]; then
                echo "Error: {wildcards.ID} is $meta (platform/layout/strategy), not ILLUMINA/PAIRED/WGS -- excluded by the pipeline's input criteria" >&2
                exit 1
            fi
            echo "{wildcards.ID} verified as ILLUMINA/PAIRED/WGS"

            # If a prior attempt got killed mid-download (walltime, node
            # failure, etc.), prefetch leaves its own
            # {wildcards.ID}/{wildcards.ID}.sra.lock file behind and refuses
            # to redownload -- every retry fails
            # immediately with "lock exists ... download canceled" even
            # though nothing is actually still running. Clear any leftover
            # state from a prior attempt before starting, so a killed job
            # doesn't permanently wedge this sample.
            if [[ -d "{wildcards.ID}" ]]; then
                echo "Found a leftover {wildcards.ID}/ directory from a prior attempt -- removing before redownloading (likely a stale prefetch lock from a job that was killed mid-download)" >&2
                rm -rf "{wildcards.ID}"
            fi

            # --max-size raised from sra-tools' 20G default; some samples
            # here exceed that and were being silently skipped by prefetch.
            prefetch --max-size 5000G {wildcards.ID}
            fastq-dump --gzip --clip --outdir ./results/raw_reads --split-3 --skip-technical ./{wildcards.ID}

            # Same failure mode as the BioSample branch: a killed job (e.g. a
            # SLURM walltime limit on a run with a large spot count) can leave
            # a non-empty but truncated gzip file that would otherwise be
            # silently accepted by every downstream rule.
            if [[ ! -s results/raw_reads/{wildcards.ID}_1.fastq.gz || ! -s results/raw_reads/{wildcards.ID}_2.fastq.gz ]]; then
                echo "Error: {wildcards.ID} finished without producing both non-empty R1/R2 files" >&2
                exit 1
            fi
            if ! gzip -t results/raw_reads/{wildcards.ID}_1.fastq.gz 2>/dev/null || ! gzip -t results/raw_reads/{wildcards.ID}_2.fastq.gz 2>/dev/null; then
                echo "Error: {wildcards.ID} produced a truncated/corrupt gzip file -- deleting partial output so a retry starts clean" >&2
                rm -f results/raw_reads/{wildcards.ID}_1.fastq.gz results/raw_reads/{wildcards.ID}_2.fastq.gz
                exit 1
            fi
        fi
        """

# Standalone gate between download and trimming. Separated from
# jules_fastq_dump itself so it can be run as its own target across every
# already-downloaded sample (`snakemake jules_verify_downloads_only
# --keep-going`) without re-triggering a full download -- useful for auditing
# raw_reads that were produced before the gzip-integrity checks in
# jules_fastq_dump existed. On failure it deletes the corrupt r1/r2 pair, so
# they no longer exist as jules_fastq_dump's declared outputs -- the next
# snakemake invocation targeting anything downstream sees them missing and
# reruns jules_fastq_dump to redownload just that sample, with no manual
# bookkeeping of which samples failed.
rule jules_verify_raw_reads:
    input:
        r1="results/raw_reads/{srr}_1.fastq.gz",
        r2="results/raw_reads/{srr}_2.fastq.gz"
    output:
        touch("results/raw_reads/{srr}.verified")
    shell:
        """
        if ! gzip -t {input.r1} 2>/dev/null || ! gzip -t {input.r2} 2>/dev/null; then
            echo "Error: {wildcards.srr} raw reads are a truncated/corrupt gzip stream -- deleting so jules_fastq_dump reruns on the next invocation" >&2
            rm -f {input.r1} {input.r2}
            exit 1
        fi
        """

rule jules_trimmomatic:
    input:
        r1="results/raw_reads/{srr}_1.fastq.gz",
        r2="results/raw_reads/{srr}_2.fastq.gz",
        verified="results/raw_reads/{srr}.verified",
        adapters=config["adapter_path"]
    output:
        r1=temp("results/trim_reads/{srr}_r1.fastq.gz"),
        r2=temp("results/trim_reads/{srr}_r2.fastq.gz"),
        u1=temp("results/trim_reads/{srr}_u1.fastq.gz"),
        u2=temp("results/trim_reads/{srr}_u2.fastq.gz")
    conda:
        "../envs/trimmomatic.yaml"
    shell:
        """
        trimmomatic PE -threads {threads} {input.r1} {input.r2} {output.r1} {output.u1} {output.r2} {output.u2} ILLUMINACLIP:{input.adapters}:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36
        """

rule jules_bwa_mem:
    input:
        r1="results/trim_reads/{srr}_r1.fastq.gz",
        r2="results/trim_reads/{srr}_r2.fastq.gz",
        ref=expand("{path_start}{ref}/{ref}.fasta", path_start = config["reference_genome_path"], ref = lookup(query = "Run == '{srr}'", within = reads, cols = "Species")),
        # read implicitly by bwa mem -- declared so jules_bwa_index isn't orphaned
        idx=expand("{path_start}{ref}/{ref}.fasta.{ext}", path_start = config["reference_genome_path"], ref = lookup(query = "Run == '{srr}'", within = reads, cols = "Species"), ext = ["amb", "ann", "bwt", "pac", "sa"])
    output:
        temp("results/align_reads/{srr}.bam")
    conda:
        "../envs/bwa.yaml"
    shell:
        """
        bwa mem -t {threads} -M -R '@RG\\tID:{wildcards.srr}\\tSM:{wildcards.srr}\\tPL:ILLUMINA\\tLB:lib1' {input.ref} {input.r1} {input.r2} | samtools sort -@ {threads} -o {output}
        """

rule jules_picard_mark_dup:
    input:
        "results/align_reads/{srr}.bam",
    output:
        bam="results/mark_reads/{srr}.bam",
        metrics="results/picard_metrics/{srr}.txt",
        tmp_dir=temp(directory("results/tmp_picard_{srr}"))
    conda:
        "../envs/picard.yaml"
    shell:
        """
        picard -Xmx28g MarkDuplicates I={input} O={output.bam} M={output.metrics} VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true TMP_DIR={output.tmp_dir}
        """

rule jules_samtools_index:
    input:
        "results/mark_reads/{srr}.bam"
    output:
        "results/mark_reads/{srr}.bam.bai"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        samtools index -@ {threads} {input}
        """

rule jules_sample_coverage:
    input:
        "results/mark_reads/{srr}.bam"
    output:
        "results/coverages/{srr}.50k.coverage.txt"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        samtools depth -a {input} | awk '{{sum+=$3}} END {{print "Average =", sum/NR}}' > {output}
        """

rule jules_bcftools_mpileup:
    input:
        ref=expand("{path_start}{ref}/{ref}.fasta", path_start = config["reference_genome_path"], ref = lookup(query = "Run == '{srr}'", within = reads, cols = "Species")),
        fai=expand("{path_start}{ref}/{ref}.fasta.fai", path_start = config["reference_genome_path"], ref = lookup(query = "Run == '{srr}'", within = reads, cols = "Species")),
        bam="results/mark_reads/{srr}.bam"
    output:
        vcf="results/bcfvcfs/{srr}_raw.vcf.gz",
        tbi="results/bcfvcfs/{srr}_raw.vcf.gz.tbi"
    conda:
        "../envs/bcftools.yaml"
    params:
        q=config["q"],
        Q=config["Q"]
    shell:
        """
        bcftools mpileup -f {input.ref} -a "FORMAT/AD,FORMAT/DP,FORMAT/SP" -q {params.q} -Q {params.Q} {input.bam} | bcftools call -mv -Oz -o {output.vcf}

        tabix {output.vcf}
        """


MIN_ACCEPTABLE_AVG_DEPTH = 10

def get_avg_cov_value(srr):
    """Average depth for a sample, parsed from its coverage file."""
    return float(open(f"results/coverages/{srr}.50k.coverage.txt").read().split()[2])

rule jules_bcftools_filter:
    input:
        vcf="results/bcfvcfs/{srr}_raw.vcf.gz",
        cov="results/coverages/{srr}.50k.coverage.txt"
    output:
        "results/bcfvcfs/{srr}_filtered.vcf.gz"
    conda:
        "../envs/bcftools.yaml"
    params:
        avg_int=lambda wc: int(get_avg_cov_value(wc.srr)),
        mincov=lambda wc: int(get_avg_cov_value(wc.srr)) // 3,
        maxcov=lambda wc: int(get_avg_cov_value(wc.srr)) * 2,
        ab_low=0.30,
        ab_high=0.70,
        min_depth=MIN_ACCEPTABLE_AVG_DEPTH
    shell:
        # Depth floor moved here from a DAG-construction-time Python
        # exception (formerly raised inside a params lambda) to a plain job
        # failure: a Python exception while Snakemake is still building the
        # DAG is fatal to the whole run regardless of --keep-going, since
        # that flag only isolates failures of jobs that actually execute.
        # At 200+ species, one too-shallow sample shouldn't take down every
        # other species' run -- exiting here instead lets --keep-going skip
        # just this sample's ROH output (and downstream jules_bcftools_roh
        # for it) while everything else proceeds.
        """
        if [ {params.avg_int} -lt {params.min_depth} ]; then
            echo "Sample {wildcards.srr} has average depth {params.avg_int}x, below the {params.min_depth}x floor -- too shallow to trust for ROH/PSMC (undercalled heterozygosity at low depth inflates F_ROH)." >&2
            exit 1
        fi

        bcftools filter -i 'QUAL>=30 && FORMAT/DP>={params.mincov} && FORMAT/DP<={params.maxcov} && INFO/DP>={params.mincov} && INFO/MQ>=30 && FORMAT/SP<60' {input.vcf} \\
            | bcftools view -v snps -m2 -M2 \\
            | bcftools filter -S . -e '(GT="het") && ((FMT/AD[0:1])/(FMT/AD[0:0]+FMT/AD[0:1]) < {params.ab_low} || (FMT/AD[0:1])/(FMT/AD[0:0]+FMT/AD[0:1]) > {params.ab_high})' \\
            -Oz -o {output}
        """

rule jules_bcftools_roh:
    input:
        "results/bcfvcfs/{srr}_filtered.vcf.gz"
    output:
        "results/roh/{srr}_ROH.txt"
    conda:
        "../envs/bcftools.yaml"
    params:
        G=config["G"],
        AFdflt=config["AFdflt"]
    shell:
        """
        bcftools roh -G{params.G} --AF-dflt {params.AFdflt} -o {output} {input}
        """

rule jules_psmc_50k_bed:
    input:
        config["reference_genome_path"] + "{ref}/{ref}.fasta.fai"
    output:
        "results/psmc_bed/{ref}.50k.bed"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        cat {input} | awk '$2>50000 {{print $1, "0", $2}}' > {output}
        """

rule jules_psmc_subset_bam:
    input:
        bam="results/mark_reads/{srr}.bam",
        bed=expand("results/psmc_bed/{species}.50k.bed", species=lookup(query="Run == '{srr}'", within=reads, cols="Species"))
    output:
        "results/psmc/{srr}.50k.bam"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools view -@ {threads} -bh -L {input.bed} \
        -o {output} {input.bam}
        """


rule jules_psmc_gen_consensus:
    input:
        ref=expand("{path_start}{ref}/{ref}.fasta", path_start = config["reference_genome_path"], ref = lookup(query = "Run == '{srr}'", within = reads, cols = "Species")),
        fai=expand("{path_start}{ref}/{ref}.fasta.fai", path_start = config["reference_genome_path"], ref = lookup(query = "Run == '{srr}'", within = reads, cols = "Species")),
        bam="results/psmc/{srr}.50k.bam",
        cov="results/coverages/{srr}.50k.coverage.txt"
    output:
        "results/psmc/{srr}.con.fq.gz"
    conda:
        "../envs/psmc_legacy.yaml"
    params:
        mincov=lambda wc: int(get_avg_cov_value(wc.srr)) // 3,
        maxcov=lambda wc: int(get_avg_cov_value(wc.srr)) * 2,
        min_depth=MIN_ACCEPTABLE_AVG_DEPTH
    shell:
        # Same depth floor as jules_bcftools_filter (ROH), so ROH, PSMC and
        # MSMC2 all share one threshold: a too-shallow genome undercalls
        # heterozygous sites, which biases PSMC's Ne estimates just as it
        # inflates FROH. Read at run time (not in params) so it cannot
        # crash DAG construction; an unreadable coverage file fails closed.
        """
        avg_depth=$(awk '{{print int($3)}}' {input.cov})
        if [[ -z "$avg_depth" ]]; then
            echo "Error: could not read average depth for {wildcards.srr} from {input.cov}" >&2
            exit 1
        fi
        if [[ "$avg_depth" -lt {params.min_depth} ]]; then
            echo "Sample {wildcards.srr} has average depth ${{avg_depth}}x, below the {params.min_depth}x floor -- too shallow to trust for PSMC (undercalled heterozygosity biases Ne)." >&2
            exit 1
        fi

        samtools mpileup -C50 -uf {input.ref} {input.bam} | \
            bcftools call -c - | \
            vcfutils.pl vcf2fq -d {params.mincov} -D {params.maxcov} | \
            gzip > {output}
        """

rule jules_psmc_gen_input:
    input:
        "results/psmc/{srr}.con.fq.gz"
    output:
        "results/psmc/{srr}.psmcfa"
    conda:
        "../envs/psmc_legacy.yaml"
    params:
        psmc_path=config["psmc_path"]
    shell:
        """
        {params.psmc_path}/utils/fq2psmcfa -q20 {input} \
        > {output}
        """

rule jules_psmc_run_psmc:
    input:
         "results/psmc/{srr}.psmcfa"
    output:
         "results/psmc/{srr}.psmc"
    conda:
         "../envs/psmc_legacy.yaml"
    params:
        psmc_path=config["psmc_path"]
    shell:
        """
        {params.psmc_path}/psmc -N25 -t15 -r5 -p "1+1+1+1+25*2+4+6" \
        -o {output} {input}
        """

rule jules_backup_bam_globus:
    input:
        bam="results/mark_reads/{srr}.bam",
        bai="results/mark_reads/{srr}.bam.bai",
        # Unused by name -- ensures this rule waits until every other
        # consumer of mark_reads/{srr}.bam has already run.
        cov="results/coverages/{srr}.50k.coverage.txt",
        vcf="results/bcfvcfs/{srr}_raw.vcf.gz",
        psmc_subset="results/psmc/{srr}.50k.bam"
    output:
        touch("results/mark_reads/{srr}.bam.backed_up")
    conda:
        "../envs/globus.yaml"
    params:
        src_endpoint=config["globus_src_endpoint"],
        dst_endpoint=config["globus_dst_endpoint"],
        dst_base=config["globus_dst_base"],
        species=lookup(query="Run == '{srr}'", within=reads, cols="Species")
    shell:
        """
        set -euo pipefail

        dst_dir="{params.dst_base}/{params.species}_aligned_to_{params.species}/mark_reads"

        # Verify each transfer's status explicitly before deleting -- task
        # wait returning just means the task finished, not that it succeeded.
        bam_task_id=$(globus transfer \\
            "{params.src_endpoint}:{input.bam}" \\
            "{params.dst_endpoint}:${{dst_dir}}/{wildcards.srr}.bam" \\
            --label "backup_{wildcards.srr}_bam" \\
            --sync-level checksum \\
            --jmespath 'task_id' --format=UNIX)
        echo "Submitted BAM transfer task ${{bam_task_id}}, waiting..."
        globus task wait "${{bam_task_id}}" --polling-interval 60 --timeout 21600
        bam_status=$(globus task show "${{bam_task_id}}" --jmespath 'status' --format=UNIX)
        if [ "$bam_status" != "SUCCEEDED" ]; then
            echo "ERROR: Globus BAM transfer ${{bam_task_id}} status=${{bam_status}} -- refusing to delete local BAM." >&2
            exit 1
        fi

        bai_task_id=$(globus transfer \\
            "{params.src_endpoint}:{input.bai}" \\
            "{params.dst_endpoint}:${{dst_dir}}/{wildcards.srr}.bam.bai" \\
            --label "backup_{wildcards.srr}_bai" \\
            --sync-level checksum \\
            --jmespath 'task_id' --format=UNIX)
        echo "Submitted BAI transfer task ${{bai_task_id}}, waiting..."
        globus task wait "${{bai_task_id}}" --polling-interval 60 --timeout 3600
        bai_status=$(globus task show "${{bai_task_id}}" --jmespath 'status' --format=UNIX)
        if [ "$bai_status" != "SUCCEEDED" ]; then
            echo "ERROR: Globus BAI transfer ${{bai_task_id}} status=${{bai_status}} -- refusing to delete local BAM." >&2
            exit 1
        fi

        echo "Both transfers verified SUCCEEDED -- deleting local BAM+index."
        rm -f {input.bam} {input.bai}
        """

