rule kmc_rm_contam:
    input:
        pread1="biosample_results/{ID}_paired_R1.fastq.gz",
        pread2="biosample_results/{ID}_paired_R2.fastq.gz",
        uread1="biosample_results/{ID}_unpaired_R1.fastq.gz",
        uread2="biosample_results/{ID}_unpaired_R2.fastq.gz",
        pre="contam.kmc_pre",
        suf="contam.kmc_suf"
    output:
        filt1=temp("no_contam_reads/{ID}_stage1.fastq"),
        list=temp("input_{ID}.txt")
    params:
        contamMatchLimitCount = config["contam_match_limit_count"],
        #contamMatchLimitPercent = config["contam_match_limit_percent"]
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
        kmc_tools -t{threads} filter contam @{output.list} -ci0 -cx{params.contamMatchLimitCount} {output.filt1}
        """
