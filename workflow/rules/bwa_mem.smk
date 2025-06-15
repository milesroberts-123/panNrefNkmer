rule bwa_mem:
    input:
        amb="../config/linear_genomes/sequence/{ref}.fa.amb",
        ann="../config/linear_genomes/sequence/{ref}.fa.ann",
        bwt="../config/linear_genomes/sequence/{ref}.fa.bwt",
        pac="../config/linear_genomes/sequence/{ref}.fa.pac",
        sa="../config/linear_genomes/sequence/{ref}.fa.sa",
        ref="../config/linear_genomes/sequence/{ref}.fa",
        read1="change_headers_results/paired_R1_{ID}.fastq.gz",
        read2="change_headers_results/paired_R2_{ID}.fastq.gz",
        readu="change_headers_results/unpaired_{ID}.fastq.gz"
    output:
        paired=temp("bwa_results/paired_{ID}_{ref}.bam"),
        paired_bai=temp("bwa_results/paired_{ID}_{ref}.bam.bai"),
        unpaired=temp("bwa_results/unpaired_{ID}_{ref}.bam"),
        unpaired_bai=temp("bwa_results/unpaired_{ID}_{ref}.bam.bai"),
        final=temp("bwa_results/{ID}_{ref}.bam"),
        bai=temp("bwa_results/{ID}_{ref}.bam.bai")
    conda:
        "../envs/bwa.yaml"
    log:
        "logs/bwa_mem/{ID}_{ref}.log",
    benchmark:
        "benchmarks/bwa_mem/{ID}_{ref}.bench",        
    shell:
        """
        # align reads to reference
        echo Aligning reads...
        bwa mem -R '@RG\\tID:{wildcards.ID}\\tSM:{wildcards.ID}' -t {threads} {input.ref} {input.read1} {input.read2} | samtools sort -O bam > {output.paired}
        bwa mem -R '@RG\\tID:{wildcards.ID}\\tSM:{wildcards.ID}' -t {threads} {input.ref} {input.readu} | samtools sort -O bam > {output.unpaired}

        # index reads
        echo Inexing alignments...
        samtools index {output.paired}
        samtools index {output.unpaired}

        # merge alignments
        echo Merging paired and unpaired alignments...
        samtools merge {output.paired} {output.unpaired} > {output.final}

        # index final output
        echo Indexing final output...
        samtools index {output.final}
        """
