# source:
# https://www.htslib.org/algorithms/duplicate.html
rule mark_dup:
    input:
        "bwa_results/{ID}_{ref}.bam"
    output:
        temp("mark_dup_results/{ID}_{ref}.bam")
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools collate -@ {threads} -O -u {input} | samtools fixmate -@ {threads} -m -u - - | samtools sort -@ {threads} -u - | samtools markdup -@ {threads} - {output}"

rule samtools_stats:
    input:
        "mark_dup_results/{ID}_{ref}.bam"
    output:
        "samtools_stats/{ID}_{ref}.txt"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools stats {input} | grep ^SN | cut -f 2- > {output}"

rule bam_index:
    input:
        "mark_dup_results/{ID}_{ref}.bam",
    output:
        temp("mark_dup_results/{ID}_{ref}.bam.bai"),
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index {input}"

rule bcftools_linref_merge:
    input:
        gz=expand("bcftools_linref_results/{ID}_{{ref}}_{{split}}.vcf.gz", ID=config["samples"] + config["outgroup"]),
        tbi=expand("bcftools_linref_results/{ID}_{{ref}}_{{split}}.vcf.gz.tbi", ID=config["samples"] + config["outgroup"])
    output:
        temp("bcftools_linref_merge_results/{ref}_{split}.vcf.gz")
    conda:
        "../envs/samtools.yaml"
    shell:
        "bcftools merge -Oz -o {output} {input.gz}"

rule bcftools_concat:
    input:
        expand("bcftools_filter_results/{{ref}}_{split}_sorted.vcf.gz", split = range(10, 10 + config["splits"]))
    output:
        gz="bcftools_concat_results/{ref}.vcf.gz",
        sorted="bcftools_concat_results/{ref}_sorted.vcf.gz",
        tbi="bcftools_concat_results/{ref}_sorted.vcf.gz.tbi"
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
        "bcftools_concat_results/{ref}_sorted.vcf.gz"
    output:
        "bcftools_stats_results/{ref}.txt"
    conda:
        "../envs/samtools.yaml"
    shell:
        "bcftools stats {input} > {output}"

rule bcftools_panref_merge:
    input:
        expand("bcftools_panref_results/{ID}_{{panref}}_{{linref}}.vcf", ID=config["samples"])
    output:
        "bcftools_panref_merge_results/{ID}_{panref}_{linref}.vcf"
    conda:
        "../envs/samtools.yaml"
    shell:
        "bcftools merge {input} > {output}"

rule bcftools_linref_call:
    input:
        bam = "mark_dup_results/{ID}_{ref}.bam",
        bai = "mark_dup_results/{ID}_{ref}.bam.bai",
        ref = "../config/linear_genomes/sequence/{ref}.fa",
        sites = "split_sites_{ref}_{split}"
    output:
        gz=temp("bcftools_linref_results/{ID}_{ref}_{split}.vcf.gz"),
        tbi=temp("bcftools_linref_results/{ID}_{ref}_{split}.vcf.gz.tbi")
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        bcftools mpileup -f {input.ref} -R {input.sites} {input.bam} | bcftools call -f GQ -m -Oz -o {output.gz}
        
        tabix {output.gz}
        """

rule bcftools_panref_call:
    input:
        bam = "vg_surject_results/{ID}_{panref}_{linref}.bam",
        ref = "vg_paths_results/{panref}.fa",
    output:
        "bcftools_panref_results/{ID}_{panref}_{linref}.vcf"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools sort -O SAM {input.bam} | bcftools mpileup -f {input.ref} -Ou - | bcftools call -mv > {output}"

rule bcftools_filter:
    input:
        "bcftools_linref_merge_results/{ref}_{split}.vcf.gz"
    output:
        invar_before = temp("bcftools_filter_results/{ref}_{split}_invar.vcf.gz"),
        var_before = temp("bcftools_filter_results/{ref}_{split}_var.vcf.gz"),
        invar_after = temp("bcftools_filter_results/{ref}_{split}_invar_filt.vcf.gz"),
        var_after = temp("bcftools_filter_results/{ref}_{split}_var_filt.vcf.gz"),
        invar_after_tbi = temp("bcftools_filter_results/{ref}_{split}_invar_filt.vcf.gz.tbi"),
        var_after_tbi = temp("bcftools_filter_results/{ref}_{split}_var_filt.vcf.gz.tbi"),
        unsorted = temp("bcftools_filter_results/{ref}_{split}_unsorted.vcf.gz"),
        sorted = temp("bcftools_filter_results/{ref}_{split}_sorted.vcf.gz")
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
        expand("../config/linear_genomes/sequence/{ref}.fa", ref = config["linrefs"])
    output:
        gz="{chr}.fa.gz",
        tbi="{chr}.fa.gz.tbi"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        echo Grab all instances of {wildcards.chr} 
        grep "{wildcards.chr}$" {input} > $(basename {output} .gz)

        echo Compress fasta file
        bgzip $(basename {output} .gz)

        echo Index fasta file
        tabix {output.gz}
        """
