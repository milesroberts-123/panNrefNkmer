rule bcftools_linref_call:
    input:
        bam = "mark_dup_results/{ID}_{ref}.bam",
        ref = "../config/linear_genomes/sequence/{ref}.fa",
        sites = "split_sites_{ref}_{split}"
    output:
        "bcftools_linref_results/{ID}_{ref}_{split}.vcf"
    conda:
        "../envs/samtools.yaml"
    shell:
        "bcftools mpileup -f {input.ref} -R {input.sites} {input.bam} | bcftools call -f GQ -m -Oz -o {output}"

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

