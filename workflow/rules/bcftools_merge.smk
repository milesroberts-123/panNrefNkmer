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

