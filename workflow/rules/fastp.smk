rule fastp:
    input:
        read1="raw_reads/{ID}_1.fastq",
        read2="raw_reads/{ID}_2.fastq",
    output:
        dpread1=temp("fastp_results/dedup_paired_R1_{ID}.fastq.gz"),
        dpread2=temp("fastp_results/dedup_paired_R2_{ID}.fastq.gz"),
        duread1=temp("fastp_results/dedup_unpaired_R1_{ID}.fastq.gz"),
        duread2=temp("fastp_results/dedup_unpaired_R2_{ID}.fastq.gz"),
        pread1=temp("fastp_results/trimmed_paired_R1_{ID}.fastq.gz"),
        pread2=temp("fastp_results/trimmed_paired_R2_{ID}.fastq.gz"),
        uread1=temp("fastp_results/trimmed_unpaired_R1_{ID}.fastq.gz"),
        uread2=temp("fastp_results/trimmed_unpaired_R2_{ID}.fastq.gz"),
        jsonR1R2=temp("fastp_results/{ID}_R1R2.json"),
        jsonU1=temp("fastp_results/{ID}_U1.json"),
        jsonU2=temp("fastp_results/{ID}_U2.json"),
    conda:
        "../envs/fastp.yaml"
    log:
        "logs/fastp/{ID}.log",
    params:
        unqualLimit=config["unqualLimit"],
        k=config["k"],
        qualThresh=config["qualThresh"],
        windowLength=config["windowLength"],
        nBaseLimit=config["nBaseLimit"]
    shell:
        """
        # remove duplicates, do read correction, drop low quality reads
        fastp --thread {threads} --n_base_limit {params.nBaseLimit} -u {params.unqualLimit} -q {params.qualThresh} --dedup --correction -i {input.read1} -I {input.read2} -o {output.dpread1} -O {output.dpread2} --unpaired1 {output.duread1} --unpaired2 {output.duread2} &>> {log}

        # trim low quality bases from remaining reads
        fastp --thread {threads} -Q -l {params.k} --cut_tail --cut_tail_window_size {params.windowLength} --cut_tail_mean_quality {params.qualThresh} --json {output.jsonR1R2} -i {output.dpread1} -I {output.dpread2} -o {output.pread1} -O {output.pread2} &>> {log}
        fastp --thread {threads} -Q -l {params.k} --cut_tail --cut_tail_window_size {params.windowLength} --cut_tail_mean_quality {params.qualThresh} --json {output.jsonU1} -i {output.duread1} -o {output.uread1} &>> {log}
        fastp --thread {threads} -Q -l {params.k} --cut_tail --cut_tail_window_size {params.windowLength} --cut_tail_mean_quality {params.qualThresh} --json {output.jsonU2} -i {output.duread2} -o {output.uread2} &>> {log}
        """

rule qc_post_rm_contam:
    input:
        "no_contam_reads/{ID}.fastq"    
    output:
        temp("fastp_results/no_contam_{ID}.json")
    shell:
        "fastp -A -Q -L -G -i {input} --json {output}"
