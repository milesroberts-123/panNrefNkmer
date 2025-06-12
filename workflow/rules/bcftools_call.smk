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

rule bam_index:
    input:
        "mark_dup_results/{ID}_{ref}.bam",
    output:
        temp("mark_dup_results/{ID}_{ref}.bam.bai"),
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index {input}"

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

