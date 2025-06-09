rule bcftools_linref_merge:
    input:
        expand("bcftools_linref_results/{ID}_{{ref}}.vcf", ID=config["samples"])
    output:
        "bcftools_linref_merge_results/{ref}.vcf"
    conda:
        "../envs/samtools.yaml"
    shell:
        "bcftools merge {input} > {output}"

rule bcftools_panref_merge:
    input:
        expand("bcftools_panref_results/{ID}_{{panref}}_{{linref}}.vcf", ID=config["samples"])
    output:
        "bcftools_panref_merge_results/{ID}_{panref}_{linref}.vcf"
    conda:
        "../envs/samtools.yaml"
    shell:
        "bcftools merge {input} > {output}"

