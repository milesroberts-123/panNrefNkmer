rule kmc:
    input:
        "results/no_contam_reads/{ID}_stage1.fastq"
    output:
        counts=temp("results/kmc/{ID}.txt"),
        tmp_pre=temp("results/kmc_db_{ID}.kmc_pre"),
        tmp_suf=temp("results/kmc_db_{ID}.kmc_suf")
    conda:
        "../envs/kmc.yaml"
    log:
        "logs/kmc/{ID}.log",
    params:
        mincount=config["mincount"],
        maxcount=config["maxcount"],
        k=config["k"],
    shell:
        """
        # create directory
        if [ -d "tmp_kmc_{wildcards.ID}" ]; then
            rm -r tmp_kmc_{wildcards.ID}
        fi

        mkdir tmp_kmc_{wildcards.ID}

        # count k-mers
        kmc -m15 -t{threads} -ci{params.mincount} -cs{params.maxcount} -k{params.k} {input} results/kmc_db_{wildcards.ID} tmp_kmc_{wildcards.ID} &>> {log}

        # dump all k-mers to text file
        kmc_tools transform results/kmc_db_{wildcards.ID} dump {output.counts} &>> {log}

        # delete tmp directories
        rm -r tmp_kmc_{wildcards.ID}
        """

rule kmc_ref_db:
    input:
        "../config/reference_genomes/{species}/{species}.fasta"
    output:
        pre=temp("results/kmc_ref_dbs/{species}.kmc_pre"),
        suf=temp("results/kmc_ref_dbs/{species}.kmc_suf")
    conda:
        "../envs/kmc.yaml"
    params:
        k = config["k"]
    shell:
        """
        if [ ! -d "kmc_tmp_{wildcards.species}" ]; then
            mkdir kmc_tmp_{wildcards.species}/
        fi

        # build kmc db
        kmc -k{params.k} -m23 -t{threads} -ci1 -cs2 -fm {input} results/kmc_ref_dbs/{wildcards.species} kmc_tmp_{wildcards.species}/

        # clean up
        rm -r kmc_tmp_{wildcards.species}/
        """

rule kmc_rm_contam:
    input:
        pread1="results/biosample/{ID}_paired_R1.fastq.gz",
        pread2="results/biosample/{ID}_paired_R2.fastq.gz",
        uread1="results/biosample/{ID}_unpaired_R1.fastq.gz",
        uread2="results/biosample/{ID}_unpaired_R2.fastq.gz",
        refdb=expand("results/kmc_ref_dbs/{species}.{ext}", ext=["kmc_pre", "kmc_suf"], species=lookup(query="BioSample == '{ID}'", within = reads, cols="Species")),
    output:
        filt1=temp("results/no_contam_reads/{ID}_stage1.fastq"),
        list=temp("results/input_{ID}.txt")
    params:
        contamMatchLimitCount = config["contam_match_limit_count"],
    conda:
        "../envs/kmc.yaml"
    shell:
        """
        if [ ! -d "results/no_contam_reads" ]; then
            mkdir -p results/no_contam_reads/
        fi

        # create file list
        echo {input.pread1} {input.pread2} {input.uread1} {input.uread2} | tr ' ' '\n' > {output.list}

        # filter reads for contamination
        kmc_tools -t{threads} filter {input.refdb[0]} @{output.list} -ci{params.contamMatchLimitCount} -cx1000000 {output.filt1}
        """
