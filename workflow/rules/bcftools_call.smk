rule bcftools_linref_call:
    input:
        bam = "bwa_results/{ID}_{ref}.bam",
        ref = "../config/linear_genomes/sequence/{ref}.fa",
        sites = "degenotate_results/{ref}/degeneracy-all-sites.bed"
    output:
        "bcftools_linref_results/{ID}_{ref}.vcf"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools sort -O SAM {input.bam} | bcftools mpileup -f {input.ref} -Ou - | bcftools call -mv > {output}"

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

