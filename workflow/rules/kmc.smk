rule kmc:
    input:
        "no_contam_reads/{ID}_stage2.fastq"
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
