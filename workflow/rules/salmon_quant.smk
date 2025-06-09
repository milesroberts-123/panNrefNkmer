rule salmon_quant:
    input:
        idx = "salmon_index_{ref}",
        read1 = "fastp_results/trimmed_paired_R1_{ID}.fastq.gz",
        read2 = "fastp_results/trimmed_paired_R2_{ID}.fastq.gz"
    output:
        "qaunts/{ID}_{ref}_quant.sf"
    params:
        prefix = "qaunts/{ID}_{ref}_quant",
        index = "salmon_index_{ref}"
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
