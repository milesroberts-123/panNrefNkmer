# MSMC2 for the jules (single-sample, linear-reference) track.
# DRAFT -- not yet included in the Snakefile.
#
# Structural echo of msmc2.smk (which compares a vg pangenome "variant"
# genome against a "backbone" genome via a graph-derived VCF), but every
# rule here is keyed on {srr} against the same linear reference the rest of
# jules_v2.smk already builds -- no vg/pggb/degenotate inputs required.
#
# Chain follows msmc-tools' documented convention (see its README):
# per-sample, per-chromosome calling mask + VCF from bamCaller.py, fed
# directly from `samtools mpileup | bcftools call -c` (NOT the
# SNP-filtered/AB-masked VCF the ROH/PSMC branch builds -- bamCaller.py
# needs the unfiltered all-sites consensus stream to do its own
# covered-vs-variant accounting) -> per-contig generate_multihetsep ->
# msmc2 run.
#
# -q/-Q reuse config["q"]/config["Q"] (the same thresholds
# jules_bcftools_mpileup uses) so quality filtering stays consistent across
# branches if those config values are ever retuned. -C50 matches both the
# msmc-tools README's own documented example and this repo's existing
# jules_psmc_gen_consensus rule, which already uses the same BAQ-downgrade
# coefficient for a consensus-calling pass over the same BWA-aligned BAMs.
#
# resources/msmc-tools is NOT fetched by a rule here -- it was cloned once
# manually (git clone --depth 1 https://github.com/stschiff/msmc-tools.git
# resources/msmc-tools) and is treated like psmc_path elsewhere in this
# repo: a pre-existing external resource referenced directly, not something
# Snakemake regenerates.

checkpoint jules_msmc2_contigs:
    input:
        fai=expand(
            "{path_start}{ref}/{ref}.fasta.fai",
            path_start=config["reference_genome_path"],
            ref=lookup(query="Run == '{srr}'", within=reads, cols="Species"),
        )
    output:
        "results/msmc2/{srr}_contigs.txt"
    shell:
        # Raised from the >50000bp floor jules_psmc_50k_bed uses --
        # confirmed on a real reference (SRR28361932's species, 795.7Mb
        # genome) that >50000bp let through 381 contigs, most of them small
        # unplaced scaffolds, and produced ~20,000 DAG steps for the msmc2
        # track alone. >500000bp gives the same 14 contigs as >1000000bp and
        # >5000000bp on that same genome (a real gap in the size
        # distribution, not an arbitrary cutoff) while still covering 93.8%
        # of total genome bp (746.4Mb of 795.7Mb) -- PSMC's own per-sample
        # cost is far lower per contig, so it can afford to keep more small
        # scaffolds; MSMC2's per-contig cost (a full mpileup/call/bamCaller
        # pass, each reserving several cpus/GB for up to 48h) can't.
        # Fallback for the edge case none of your actual references should
        # hit (all scaffold-level-or-better, so N50 should sit well above
        # 500000bp) -- if a species somehow has zero contigs that large,
        # take its single biggest contig instead of silently producing an
        # empty list and leaving jules_msmc2_run with no input files.
        """
        mkdir -p results/msmc2
        awk '$2>500000 {{print $1}}' {input.fai} > {output}
        if [ ! -s {output} ]; then
            sort -k2,2nr {input.fai} | head -1 | cut -f1 > {output}
        fi
        """

rule jules_msmc2_call:
    input:
        ref=expand(
            "{path_start}{ref}/{ref}.fasta",
            path_start=config["reference_genome_path"],
            ref=lookup(query="Run == '{srr}'", within=reads, cols="Species"),
        ),
        fai=expand(
            "{path_start}{ref}/{ref}.fasta.fai",
            path_start=config["reference_genome_path"],
            ref=lookup(query="Run == '{srr}'", within=reads, cols="Species"),
        ),
        bam="results/mark_reads/{srr}.bam",
        cov="results/coverages/{srr}.50k.coverage.txt",
        tools="resources/msmc-tools"
    output:
        vcf="results/msmc2/{srr}/{chrom}.vcf.gz",
        mask="results/msmc2/{srr}/{chrom}_mask.bed.gz"
    conda:
        "../envs/msmc2.yaml"
    params:
        q=config["q"],
        Q=config["Q"],
        mean_cov=lambda wc: get_avg_cov_value(wc.srr)
    shell:
        # msmc-tools' README documents `samtools mpileup -u`, but modern
        # samtools (pinned 1.22 in envs/msmc2.yaml) has removed BCF/VCF
        # output from `samtools mpileup` entirely -- it silently produces
        # plain-text pileup instead, which `bcftools call` then can't parse
        # ("Failed to read from standard input: unknown file type").
        # `bcftools mpileup -Ou` is the modern replacement (same -C/-q/-Q/-r/-f
        # flags), and matches what jules_bcftools_mpileup already uses
        # elsewhere in this repo.
        # bamCaller.py only instantiates its MaskGenerator (the only thing
        # that ever writes output.mask) upon seeing the FIRST non-header
        # line from bcftools call. A contig with zero reads -- entirely
        # possible on a real genome, not just this subset test -- produces
        # header-only input, so the mask file never gets created at all and
        # Snakemake fails with MissingOutputException even though the shell
        # pipeline itself exited 0. Pre-create a valid empty gzip file as a
        # fallback; bamCaller.py overwrites it normally whenever there IS
        # real coverage.
        """
        mkdir -p results/msmc2/{wildcards.srr}
        echo -n | gzip -c > {output.mask}
        bcftools mpileup -C50 -q {params.q} -Q {params.Q} -Ou -r {wildcards.chrom} -f {input.ref} {input.bam} \
            | bcftools call -c -V indels \
            | python {input.tools}/bamCaller.py {params.mean_cov} {output.mask} \
            | gzip -c > {output.vcf}
        """

rule jules_msmc2_contig:
    input:
        vcf="results/msmc2/{srr}/{chrom}.vcf.gz",
        mask="results/msmc2/{srr}/{chrom}_mask.bed.gz",
        tools="resources/msmc-tools"
    output:
        "results/msmc2/{srr}/{chrom}.multihetsep.txt"
    conda:
        "../envs/msmc2.yaml"
    shell:
        # No --chr and no pre-splitting needed here: jules_msmc2_call
        # already restricted mpileup to this one chromosome via -r, so the
        # VCF/mask this rule sees are single-chromosome from the start --
        # generate_multihetsep.py's monotonically-increasing position
        # counter (see msmc-tools source) is never at risk of a cross-contig
        # reset the way it would be from a whole-genome VCF/mask.
        """
        python {input.tools}/generate_multihetsep.py --mask {input.mask} {input.vcf} > {output}
        """

def jules_msmc2_inputs(wildcards):
    contigs_file = checkpoints.jules_msmc2_contigs.get(srr=wildcards.srr).output[0]
    with open(contigs_file) as f:
        chroms = [line.strip() for line in f if line.strip()]
    return expand(
        "results/msmc2/{srr}/{chrom}.multihetsep.txt",
        srr=wildcards.srr,
        chrom=chroms,
    )

rule jules_msmc2_run:
    input:
        jules_msmc2_inputs
    output:
        "results/msmc2/{srr}/msmc2.final.txt"
    conda:
        "../envs/msmc2.yaml"
    params:
        # No guessed fallback pattern -- only override -p if msmc2_p is
        # explicitly set in config; otherwise trust msmc2's own built-in
        # default time-segmentation rather than substituting a guess for it.
        p_flag=lambda wc: f"-p {config['msmc2_p']}" if config.get("msmc2_p") else "",
        prefix=lambda wc: f"results/msmc2/{wc.srr}/msmc2"
    threads: 4
    shell:
        # bioconda's msmc2 package installs its binary as msmc2_Linux, not
        # msmc2 -- confirmed via `ls .snakemake/conda/.../bin/ | grep msmc`.
        #
        # msmc2_Linux crashes (core.exception.ArrayIndexError@psmc_hmm.d:124,
        # "index [0] is out of bounds for array of length 0") if handed a
        # multihetsep.txt with zero segregating sites -- which is exactly
        # what a zero-coverage contig produces (see jules_msmc2_call). All
        # 74 declared inputs still have to exist for Snakemake, they just
        # aren't all worth feeding to msmc2 itself, so filter to non-empty
        # files at the shell level before invoking it.
        """
        FILES=""
        for f in {input}; do
            if [ -s "$f" ]; then
                FILES="$FILES $f"
            fi
        done
        msmc2_Linux -t {threads} {params.p_flag} -o {params.prefix} $FILES
        """
