rule bwa_mem:
    input:
        amb="reference_genomes/{ref}.fasta.amb",
        ann="reference_genomes/{ref}.fasta.ann",
        bwt="reference_genomes/{ref}.fasta.bwt",
        pac="reference_genomes/{ref}.fasta.pac",
        sa="reference_genomes/{ref}.fasta.sa",
        reffasta="reference_genomes/{ref}.fasta",
        read1="fastp_results/trimmed_paired_R1_{ID}.fastq",
        read2="fastp_results/trimmed_paired_R2_{ID}.fastq",
    output:
        temp("bwa_results/{ID}_{ref}.bam")
    conda:
        "../envs/bwa.yaml"
    log:
        "logs/bwa_mem/{ID}_{PID}.log",
    shell:
        """
        # align reads to reference
        bwa mem -R '@RG\\tID:{wildcards.SID}\\tSM:{wildcards.ID}' -t {threads} {input.reffasta} {input.read1} {input.read2} | samtools sort -O bam > {output}
        """
