rule fastp:
    input:
        read1="raw_reads/{ID}_1.fastq",
        read2="raw_reads/{ID}_2.fastq",
    output:
        pread1=temp("fastp_results/trimmed_paired_R1_{ID}.fastq.gz"),
        pread2=temp("fastp_results/trimmed_paired_R2_{ID}.fastq.gz"),
        uread1=temp("fastp_results/trimmed_unpaired_R1_{ID}.fastq.gz"),
        uread2=temp("fastp_results/trimmed_unpaired_R2_{ID}.fastq.gz"),
        json=temp("fastp_results/{ID}.json"),
    conda:
        "../envs/fastp.yaml"
    params:
        unqualLimit=config["unqualLimit"],
        k=config["k"],
        qualThresh=config["qualThresh"],
        windowLength=config["windowLength"],
        nBaseLimit=config["nBaseLimit"]
    shell:
        """
        fastp --thread {threads} --n_base_limit {params.nBaseLimit} -u {params.unqualLimit} -q {params.qualThresh} -l {params.k} --cut_tail --cut_tail_window_size {params.windowLength} --cut_tail_mean_quality {params.qualThresh} --dedup --correction -i {input.read1} -I {input.read2} -o {output.pread1} -O {output.pread2} --unpaired1 {output.uread1} --unpaired2 {output.uread2} --json {output.json}
        """

rule qc_post_rm_contam:
    input:
        "no_contam_reads/{ID}_stage1.fastq"    
    output:
        temp("fastp_results/no_contam_{ID}.json")
    conda:
        "../envs/fastp.yaml"
    shell:
        "fastp -A -Q -L -G -i {input} --json {output}"
