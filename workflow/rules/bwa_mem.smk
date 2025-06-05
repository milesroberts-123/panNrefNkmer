rule bwa_mem:
    input:
        amb="../config/linear_genomes/{ref}.fa.amb",
        ann="../config/linear_genomes/{ref}.fa.ann",
        bwt="../config/linear_genomes/{ref}.fa.bwt",
        pac="../config/linear_genomes/{ref}.fa.pac",
        sa="../config/linear_genomes/{ref}.fa.sa",
        reffasta="../config/linear_genomes/{ref}.fa",
        read1="fastp_results/trimmed_paired_R1_{ID}.fastq.gz",
        read2="fastp_results/trimmed_paired_R2_{ID}.fastq.gz"
    output:
        temp("bwa_results/{ID}_{ref}.bam")
    conda:
        "../envs/bwa.yaml"
    log:
        "logs/bwa_mem/{ID}_{ref}.log",
    benchmark:
        "benchmarks/bwa_mem/{ID}_{ref}.bench",        
    shell:
        """
        # align reads to reference
        bwa mem -R '@RG\\tID:{wildcards.ID}\\tSM:{wildcards.ID}' -t {threads} {input.reffasta} {input.read1} {input.read2} | samtools sort -O bam > {output}
        """
