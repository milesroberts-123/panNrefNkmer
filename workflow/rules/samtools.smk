# source:
# https://www.htslib.org/algorithms/duplicate.html
rule mark_dup:
    input:
        "results/bwa/{ID}_{ref}.bam"
    output:
        temp("results/mark_dup/{ID}_{ref}.bam")
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools collate -@ {threads} -O -u {input} | samtools fixmate -@ {threads} -m -u - - | samtools sort -@ {threads} -u - | samtools markdup -@ {threads} - {output}"

rule samtools_stats:
    input:
        "results/mark_dup/{ID}_{ref}.bam"
    output:
        "results/samtools_stats/{ID}_{ref}.txt"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools stats {input} | grep ^SN | cut -f 2- > {output}"

rule bam_index:
    input:
        "results/mark_dup/{ID}_{ref}.bam",
    output:
        temp("results/mark_dup/{ID}_{ref}.bam.bai"),
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index {input}"

rule bcftools_linref_merge:
    input:
        gz=expand("results/bcftools_linref/{ID}_{{ref}}_{{split}}.vcf.gz", ID=reads["BioSample"].unique()),
        tbi=expand("results/bcftools_linref/{ID}_{{ref}}_{{split}}.vcf.gz.tbi", ID=reads["BioSample"].unique())
    output:
        temp("results/bcftools_linref_merge/{ref}_{split}.vcf.gz")
    conda:
        "../envs/samtools.yaml"
    shell:
        "bcftools merge -Oz -o {output} {input.gz}"

rule bcftools_concat:
    input:
        expand("results/bcftools_filter/{{ref}}_{split}_sorted.vcf.gz", split = range(10, 10 + config["splits"]))
    output:
        gz="results/bcftools_concat/{ref}.vcf.gz",
        sorted="results/bcftools_concat/{ref}_sorted.vcf.gz",
        tbi="results/bcftools_concat/{ref}_sorted.vcf.gz.tbi"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        bcftools concat --allow-overlaps -Oz -o {output.gz} {input}

        bcftools sort -Oz -o {output.sorted} {output.gz}

        tabix {output.sorted}
        """

rule bcftools_stats:
    input:
        "results/bcftools_concat/{ref}_sorted.vcf.gz"
    output:
        "results/bcftools_stats/{ref}.txt"
    conda:
        "../envs/samtools.yaml"
    shell:
        "bcftools stats {input} > {output}"

rule bcftools_panref_merge:
    input:
        expand("results/bcftools_panref/{ID}_{{panref}}_{{linref}}.vcf", ID=reads[reads["Group"] == "ingroup"]["Run"].tolist())
    output:
        "results/bcftools_panref_merge/{ID}_{panref}_{linref}.vcf"
    conda:
        "../envs/samtools.yaml"
    shell:
        "bcftools merge {input} > {output}"

rule bcftools_linref_call:
    input:
        bam = "results/mark_dup/{ID}_{ref}.bam",
        bai = "results/mark_dup/{ID}_{ref}.bam.bai",
        ref = "../config/linear_genomes/sequence/{ref}.fa",
        sites = "results/split/split_sites_{ref}_{split}"
    output:
        gz=temp("results/bcftools_linref/{ID}_{ref}_{split}.vcf.gz"),
        tbi=temp("results/bcftools_linref/{ID}_{ref}_{split}.vcf.gz.tbi")
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        bcftools mpileup -f {input.ref} -R {input.sites} {input.bam} | bcftools call -f GQ -m -Oz -o {output.gz}
        
        tabix {output.gz}
        """

rule bcftools_panref_call:
    input:
        bam = "results/vg_surject/{ID}_{panref}_{linref}.bam",
        ref = "results/vg_paths/{panref}.fa",
    output:
        "results/bcftools_panref/{ID}_{panref}_{linref}.vcf"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools sort -O bam {input.bam} | bcftools mpileup -f {input.ref} -Ou - | bcftools call -mv > {output}"

rule bcftools_filter:
    input:
        "results/bcftools_linref_merge/{ref}_{split}.vcf.gz"
    output:
        invar_before = temp("results/bcftools_filter/{ref}_{split}_invar.vcf.gz"),
        var_before = temp("results/bcftools_filter/{ref}_{split}_var.vcf.gz"),
        invar_after = temp("results/bcftools_filter/{ref}_{split}_invar_filt.vcf.gz"),
        var_after = temp("results/bcftools_filter/{ref}_{split}_var_filt.vcf.gz"),
        invar_after_tbi = temp("results/bcftools_filter/{ref}_{split}_invar_filt.vcf.gz.tbi"),
        var_after_tbi = temp("results/bcftools_filter/{ref}_{split}_var_filt.vcf.gz.tbi"),
        unsorted = temp("results/bcftools_filter/{ref}_{split}_unsorted.vcf.gz"),
        sorted = temp("results/bcftools_filter/{ref}_{split}_sorted.vcf.gz")
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        # separate invariant vs variant (bi-allelic snp) sites
        echo Separating invariant and variant sites
        bcftools view --max-af 0 -Oz -o {output.invar_before} {input}
        bcftools view --types snps -m 2 -M 2 -Oz -o {output.var_before} {input}

        # filter variant sites
        echo Filtering variant sites...
        bcftools filter -i 'FMT/AD>=5 & FMT/GQ>=20 & F_MISSING<=0.8' -Oz -o {output.var_after} {output.var_before}

        # filter invariant sites
        echo Filtering invariant sites...
        bcftools filter -i 'FMT/AD>=5 & F_MISSING<=0.8' -Oz -o {output.invar_after} {output.invar_before}

        # index filtered sites
        echo Indexing filtered sites...
        tabix {output.invar_after}
        tabix {output.var_after}

        # recombine filtered sites
        echo Recombining sites after filtering...
        bcftools concat --allow-overlaps -Oz -o {output.unsorted} {output.var_after} {output.invar_after}
        bcftools sort -Oz -o {output.sorted} {output.unsorted}
        tabix {output.sorted}
        """

rule separate_chrom:
    input:
        expand("../config/linear_genomes/sequence/{ref}.fa", ref = linrefs)
    output:
        gz="results/separate_chrom/{chr}.fa.gz",
        tbi="results/separate_chrom/{chr}.fa.gz.tbi"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        echo Grab all instances of {wildcards.chr} 
        grep "{wildcards.chr}$" {input} > $(basename {output.gz} .gz)

        echo Compress fasta file
        bgzip $(basename {output.gz} .gz)

        echo Index fasta file
        tabix {output.gz}
        """
