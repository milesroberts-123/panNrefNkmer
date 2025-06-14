rule kmc_rm_contam:
    input:
        pread1="fastp_results/trimmed_paired_R1_{ID}.fastq.gz",
        pread2="fastp_results/trimmed_paired_R2_{ID}.fastq.gz",
        uread1="fastp_results/trimmed_unpaired_R1_{ID}.fastq.gz",
        uread2="fastp_results/trimmed_unpaired_R2_{ID}.fastq.gz",
        pre="contam.kmc_pre",
        suf="contam.kmc_suf"
    output:
        filt=temp("no_contam_reads/{ID}.fastq"),
        list=temp("input_{ID}.txt")
    params:
        contamMatchLimit = config["contam_match_limit"]
    conda:
        "../envs/kmc.yaml"
    shell:
        """
        if [ ! -d "no_contam_reads" ]; then
            mkdir no_contam_reads/
        fi

        # create file list
        echo {input.pread1} {input.pread2} {input.uread1} {input.uread2} | tr ' ' '\n' > {output.list}

        # filter reads for contamination
        kmc_tools -t{threads} filter contam @{output.list} -ci0.0 -cx{params.contamMatchLimit} {output.filt}

        # compress reads
        # gzip no_contam_reads/{wildcards.ID}.fastq
        """
