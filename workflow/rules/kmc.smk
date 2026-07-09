rule kmc:
    input:
        "no_contam_reads/{ID}_stage1.fastq"
    output:
        counts=temp("kmc_results/{ID}.txt"),
        tmp_pre=temp("kmc_db_{ID}.kmc_pre"),
        tmp_suf=temp("kmc_db_{ID}.kmc_suf")
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
        kmc -m15 -t{threads} -ci{params.mincount} -cs{params.maxcount} -k{params.k} {input} kmc_db_{wildcards.ID} tmp_kmc_{wildcards.ID} &>> {log}

        # dump all k-mers to text file
        kmc_tools transform kmc_db_{wildcards.ID} dump {output.counts} &>> {log}

        # delete tmp directories
        rm -r tmp_kmc_{wildcards.ID}
        """

#rule kmc_contam_db:
#    input:
#        expand("../config/contaminants/{contam}.fa" , contam = config["custom_contams"]),
#        "contam.fa"
#    output:
#        list="contam.list",
#        pre="contam.kmc_pre",
#        suf="contam.kmc_suf"
#    conda:
#        "../envs/kmc.yaml"
#    params:
#        k = config["k"]
#    shell:
#        """
#        if [ ! -d "kmc_tmp_dir" ]; then
#            mkdir kmc_tmp_dir/
#        fi
#
#        # create file list
#        echo {input} | tr ' ' '\n' > {output.list}
#
#        # build kmc db
#        kmc -k{params.k} -m23 -t{threads} -ci1 -cs2 -fm @{output.list} contam kmc_tmp_dir/
#
#        # clean up
#        rm -r kmc_tmp_dir/
#        """

rule kmc_ref_db:
    input:
        "../config/reference_genomes/{species}/{species}.fasta"
    output:
        pre=temp("kmc_ref_dbs/{species}.kmc_pre"),
        suf=temp("kmc_ref_dbs/{species}.kmc_suf")
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
        kmc -k{params.k} -m23 -t{threads} -ci1 -cs2 -fm {input} kmc_ref_dbs/{wildcards.species} kmc_tmp_{wildcards.species}/

        # clean up
        rm -r kmc_tmp_{wildcards.species}/
        """

rule kmc_rm_contam:
    input:
        pread1="biosample_results/{ID}_paired_R1.fastq.gz",
        pread2="biosample_results/{ID}_paired_R2.fastq.gz",
        uread1="biosample_results/{ID}_unpaired_R1.fastq.gz",
        uread2="biosample_results/{ID}_unpaired_R2.fastq.gz",
        refdb=expand("kmc_ref_dbs/{species}.{ext}", ext=["kmc_pre", "kmc_suf"], species=lookup(query="BioSample == '{ID}'", within = reads, cols="Species")),
    output:
        filt1=temp("no_contam_reads/{ID}_stage1.fastq"),
        list=temp("input_{ID}.txt")
    params:
        contamMatchLimitCount = config["contam_match_limit_count"],
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
        kmc_tools -t{threads} filter $(basename {input.refdb[0]} .kmc_pre) @{output.list} -ci{params.contamMatchLimitCount} -cx1000000 {output.filt1}
        """
