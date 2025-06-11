rule salmon_quant:
    input:
        idx = "salmon_index_{ref}.done",
        read1 = "fastp_results/trimmed_paired_R1_{ID}.fastq.gz",
        read2 = "fastp_results/trimmed_paired_R2_{ID}.fastq.gz"
    output:
        "salmon_quant_results/{ref}_{ID}/quant.sf"
    params:
        prefix = "salmon_quant_results/{ref}_{ID}",
        index = "salmon_index_{ref}"
    conda:
        "../envs/salmon.yaml"
    shell:
        """
        salmon quant -i {params.index} \
            -l A \
            -1 {input.read1} \
            -2 {input.read2} \
            -p {threads} \
            --validateMappings \
            -o {params.prefix}
        """
