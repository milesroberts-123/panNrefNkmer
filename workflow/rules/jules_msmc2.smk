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
        """
        mkdir -p results/msmc2
        cut -f1 {input.fai} > {output}
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
        """
        mkdir -p results/msmc2/{wildcards.srr}
        samtools mpileup -C50 -q {params.q} -Q {params.Q} -u -r {wildcards.chrom} -f {input.ref} {input.bam} \
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
        p=config.get("msmc2_p", "1*2+25*1+1*2+1*3"),
        prefix=lambda wc: f"results/msmc2/{wc.srr}/msmc2"
    threads: 4
    shell:
        """
        msmc2 -t {threads} -p {params.p} -o {params.prefix} {input}
        """
