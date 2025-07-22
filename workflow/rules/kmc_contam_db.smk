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
        kmc -k{params.k} -m23 -t{threads} -ci1 -cs2 -fm @{output.list} contam kmc_tmp_dir/

        # clean up
        rm -r kmc_tmp_dir/
        """
