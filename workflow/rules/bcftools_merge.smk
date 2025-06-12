rule bcftools_linref_merge:
    input:
        gz=expand("bcftools_linref_results/{ID}_{{ref}}_{{split}}.vcf.gz", ID=config["samples"]),
        tbi=expand("bcftools_linref_results/{ID}_{{ref}}_{{split}}.vcf.tbi", ID=config["samples"])
    output:
        "bcftools_linref_merge_results/{ref}_{split}.vcf.gz"
    conda:
        "../envs/samtools.yaml"
    shell:
        "bcftools merge -Oz -o {output} {input.gz}"

rule bcftools_concat:
    input:
        expand("bcftools_linref_merge_results/{{ref}}_{split}.vcf.gz", split = range(10, 10 + config["splits"]))
    output:
        "bcftools_concat_results/{ref}.vcf.gz"
    conda:
        "../envs/samtools.yaml"
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

