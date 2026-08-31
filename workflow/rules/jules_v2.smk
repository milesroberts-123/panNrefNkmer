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
            uids=$(curl -s --max-time 120 \\
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term={wildcards.ID}&retmax=500" \\
                | grep -oE "<Id>[0-9]+</Id>" | sed 's/<[^>]*>//g' | paste -sd, -)

            if [[ -z "$uids" ]]; then
                echo "Error: no SRA records found for BioSample {wildcards.ID}" >&2
                exit 1
            fi

            runs=$(curl -s --max-time 120 \\
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id=${{uids}}&rettype=runinfo&retmode=text" \\
                | cut -d',' -f1 | grep -E '^[SED]RR' || true)

            if [[ -z "$runs" ]]; then
                echo "Error: no SRA runs found for BioSample {wildcards.ID}" >&2
                exit 1
            fi
            echo "Found runs for {wildcards.ID}: ${{runs}}"

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

                if [[ -s "$r1" && -s "$r2" ]]; then
                    echo "--- run ${{run}} already downloaded (from a prior attempt), skipping ---"
                    return 0
                fi

                echo "--- downloading run ${{run}} (part of BioSample {wildcards.ID}) ---"
                prefetch --max-size 5000G -O "$dir" "$run" \\
                    && fastq-dump --gzip --clip --outdir "$dir" --split-3 --skip-technical "${{dir}}/${{run}}/${{run}}.sra"

                if [[ ! -s "$r1" || ! -s "$r2" ]]; then
                    echo "Error: run ${{run}} finished without producing both non-empty R1/R2 files" >&2
                    return 1
                fi
            }}
            export -f download_one_run

            echo "$runs" | xargs -P 4 -I{{}} bash -c 'download_one_run "$@"' _ {{}} "$staging_dir" \\
                || {{ echo "Error: one or more constituent SRA run downloads failed for BioSample {wildcards.ID} -- staging dir left in place at ${{staging_dir}} for the next retry to resume from" >&2; exit 1; }}

            # Every run must have produced both R1 and R2 -- catches a run
            # turning out to be single-end, which would otherwise silently
            # desync R1/R2 pairing downstream.
            n_r1=$(find "$staging_dir" -maxdepth 1 -name "*_1.fastq.gz" | wc -l)
            n_r2=$(find "$staging_dir" -maxdepth 1 -name "*_2.fastq.gz" | wc -l)
            if [[ "$n_r1" -ne "$n_r2" ]]; then
                echo "Error: mismatched R1 (${{n_r1}}) vs R2 (${{n_r2}}) file counts in ${{staging_dir}} for BioSample {wildcards.ID} -- likely a single-end run mixed into paired-end data. Staging dir left in place for inspection." >&2
                exit 1
            fi

            cat "${{staging_dir}}"/*_1.fastq.gz > results/raw_reads/{wildcards.ID}_1.fastq.gz
            cat "${{staging_dir}}"/*_2.fastq.gz > results/raw_reads/{wildcards.ID}_2.fastq.gz

            # Verify the combined output before deleting the only copies of
            # the source data.
            if [[ ! -s results/raw_reads/{wildcards.ID}_1.fastq.gz || ! -s results/raw_reads/{wildcards.ID}_2.fastq.gz ]]; then
                echo "Error: concatenation produced an empty/missing output for {wildcards.ID} -- staging dir left in place at ${{staging_dir}}" >&2
                exit 1
            fi
            rm -rf "$staging_dir"
        else
            # Plain Run accession (SRR/ERR/DRR) -- single download.
            # --max-size raised from sra-tools' 20G default; some samples
            # here exceed that and were being silently skipped by prefetch.
            prefetch --max-size 5000G {wildcards.ID}
            fastq-dump --gzip --clip --outdir ./results/raw_reads --split-3 --skip-technical ./{wildcards.ID}
        fi
        """

rule jules_trimmomatic:
    input:
        r1="results/raw_reads/{srr}_1.fastq.gz",
        r2="results/raw_reads/{srr}_2.fastq.gz",
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
        maxcov=lambda wc: int(get_avg_cov_value(wc.srr)) * 2
    shell:
        """
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

