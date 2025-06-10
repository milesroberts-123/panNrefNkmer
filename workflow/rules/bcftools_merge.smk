rule bcftools_linref_merge:
    input:
        expand("bcftools_linref_results/{ID}_{{ref}}_{{chrom}}.vcf", ID=config["samples"])
    output:
        "bcftools_linref_merge_results/{ref}_{chrom}.vcf"
    conda:
        "../envs/samtools.yaml"
    shell:
        "bcftools merge {input} > {output}"

rule bcftools_concat:
    input:
        expand("bcftools_linref_merge_results/{{ref}}_{chrom}.vcf", chrom = config["chroms"])
    output:
        "bcftools_concat_results/{ref}.vcf"
    shell:
        "bcftools concat {input} > {output}"

rule bcftools_panref_merge:
    input:
        expand("bcftools_panref_results/{ID}_{{panref}}_{{linref}}.vcf", ID=config["samples"])
    output:
        "bcftools_panref_merge_results/{ID}_{panref}_{linref}.vcf"
    conda:
        "../envs/samtools.yaml"
    shell:
        "bcftools merge {input} > {output}"

