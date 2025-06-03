rule kmc:
    input:
        pread1="fastp_results/trimmed_paired_R1_{ID}.fastq",
        pread2="fastp_results/trimmed_paired_R2_{ID}.fastq",
        uread1="fastp_results/trimmed_unpaired_R1_{ID}.fastq",
        uread2="fastp_results/trimmed_unpaired_R2_{ID}.fastq",
    output:
        counts=temp("kmc_results/kmer_counts_{ID}.txt"),
        tmp_R1_pre=temp("tmp_R1_{ID}.kmc_pre"),
        tmp_R1_suf=temp("tmp_R1_{ID}.kmc_suf"),
        tmp_R2_pre=temp("tmp_R2_{ID}.kmc_pre"),
        tmp_R2_suf=temp("tmp_R2_{ID}.kmc_suf"),
        tmp_u_R1_pre=temp("tmp_u_R1_{ID}.kmc_pre"),
        tmp_u_R1_suf=temp("tmp_u_R1_{ID}.kmc_suf"),
        tmp_u_R2_pre=temp("tmp_u_R2_{ID}.kmc_pre"),
        tmp_u_R2_suf=temp("tmp_u_R2_{ID}.kmc_suf"),
        union_R1_R2_pre=temp("union_R1_R2_{ID}.kmc_pre"),
        union_R1_R2_suf=temp("union_R1_R2_{ID}.kmc_suf"),
        union_R1_R2_u1_pre=temp("union_R1_R2_u1_{ID}.kmc_pre"),
        union_R1_R2_u1_suf=temp("union_R1_R2_u1_{ID}.kmc_suf"),
        union_R1_R2_u1_u2_pre=temp("union_R1_R2_u1_u2_{ID}.kmc_pre"),
        union_R1_R2_u1_u2_suf=temp("union_R1_R2_u1_u2_{ID}.kmc_suf"),
    conda:
        "../envs/kmc.yaml"
    log:
        "logs/kmc/{ID}.log",
    benchmark:
        "benchmarks/kmc/{ID}.bench"
    params:
        mincount=config["mincount"],
        maxcount=config["maxcount"],
        k=config["k"],
    shell:
        """
        # create directory
        if [ -d "tmp_kmc_{wildcards.SID}_{wildcards.PID}" ]; then
            rm -r tmp_kmc_{wildcards.SID}_{wildcards.PID}
        fi

        mkdir tmp_kmc_{wildcards.SID}_{wildcards.PID}

        # count k-mers
        kmc -t{threads} -ci{params.mincount} -cs{params.maxcount} -k{params.k} {input.pread1} tmp_R1_{wildcards.SID}_{wildcards.PID} tmp_kmc_{wildcards.SID}_{wildcards.PID} &>> {log}
        kmc -t{threads} -ci{params.mincount} -cs{params.maxcount} -k{params.k} {input.pread2} tmp_R2_{wildcards.SID}_{wildcards.PID} tmp_kmc_{wildcards.SID}_{wildcards.PID} &>> {log}
        kmc -t{threads} -ci{params.mincount} -cs{params.maxcount} -k{params.k} {input.uread1} tmp_u_R1_{wildcards.SID}_{wildcards.PID} tmp_kmc_{wildcards.SID}_{wildcards.PID} &>> {log}
        kmc -t{threads} -ci{params.mincount} -cs{params.maxcount} -k{params.k} {input.uread2} tmp_u_R2_{wildcards.SID}_{wildcards.PID} tmp_kmc_{wildcards.SID}_{wildcards.PID} &>> {log}

        # combine k-mer counts into one database
        kmc_tools simple tmp_R1_{wildcards.SID}_{wildcards.PID} tmp_R2_{wildcards.SID}_{wildcards.PID} union union_R1_R2_{wildcards.SID}_{wildcards.PID} &>> {log}
        kmc_tools simple union_R1_R2_{wildcards.SID}_{wildcards.PID} tmp_u_R1_{wildcards.SID}_{wildcards.PID} union union_R1_R2_u1_{wildcards.SID}_{wildcards.PID} &>> {log}
        kmc_tools simple union_R1_R2_u1_{wildcards.SID}_{wildcards.PID} tmp_u_R2_{wildcards.SID}_{wildcards.PID} union union_R1_R2_u1_u2_{wildcards.SID}_{wildcards.PID} &>> {log}

        # dump all k-mers to text file
        kmc_tools transform union_R1_R2_u1_u2_{wildcards.SID}_{wildcards.PID} dump {output.counts} &>> {log}

        # delete tmp directories
        rm -r tmp_kmc_{wildcards.SID}_{wildcards.PID}
        """
