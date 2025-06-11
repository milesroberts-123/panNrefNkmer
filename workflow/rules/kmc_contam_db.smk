rule kmc_contam_db:
    input:
        expand("../config/contaminants/{contam}.fa" , contam = config["custom_contams"]),
        "contam.fa"
    output:
        list="contam.list",
        pre="contam.kmc_pre",
        suf="contam.kmc_suf"
    conda:
        "../envs/kmc.yaml"
    params:
        k = config["k"]
    shell:
        """
        if [ ! -d "kmc_tmp_dir" ]; then
            mkdir kmc_tmp_dir/
        fi

        # create file list
        echo {input} | tr ' ' '\n' > {output.list}

        # build kmc db
        kmc -k{params.k} -m7 -t{threads} -fm @{output.list} contam kmc_tmp_dir/

        # clean up
        rm -r kmc_tmp_dir/
        """
