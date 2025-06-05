rule vg_giraffe:
    input:
        read1 = "fastp_results/trimmed_paired_R1_{ID}.fastq.gz",
        read2 = "fastp_results/trimmed_paired_R1_{ID}.fastq.gz",
        dist = "../config/pangenomes/{ref}.dist",
        gbz = "../config/pangenomes/{ref}.giraffe.gbz",
        min = "../config/pangenomes/{ref}.shortread.withzip.min",
        zip = "../config/pangenomes/{ref}.shortread.zipcodes"
    output:
        "vg_giraffe_results/{ID}_{ref}.gam"
    conda:
        "../envs/vg.yaml"
    benchmark:
        "benchmarks/vg_giraffe/{ref}/{ID}.log"
    shell:
        """
        vg giraffe -p -t {threads} -Z {input.gbz} -m {input.min} -z {input.zip} -d {input.dist} -f {input.read1} -f {input.read2} > {output}
        """
