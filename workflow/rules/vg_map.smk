rule vg_map:
    input:
        read1 = "fastp_results/trimmed_paired_R1_{ID}.fastq.gz",
        read2 = "fastp_results/trimmed_paired_R1_{ID}.fastq.gz",
        xg = "{ref}.xg",
        gcsa = "{ref}.gcsa",
    output:
        "{ID}_{ref}.gam"
    conda:
        "../envs/vg.yaml"
    shell:
        "vg map -f {input.read1} -f {input.read2} -x {input.xg} -g {input.gcsa} > {output}"
