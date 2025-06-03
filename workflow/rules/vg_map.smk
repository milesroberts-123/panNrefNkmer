rule vg_map:
    input:
        read1 = "{ID}_R1.fastq",
        read2 = "{ID}_R2.fastq",
        xg = "{ref}.xg",
        gcsa = "{ref}.gcsa",
    output:
        "{ID}_{ref}.gam"
    conda:
        "../envs/vg.yaml"
    shell:
        "vg map -f {input.read1} -f {input.read2} -x {input.xg} -g {input.gcsa} > {output}"
