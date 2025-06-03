rule kmc_contam_db:
    input:
        expand("../config/contaminants/{contam}.fa" , contam = config["contams"])
    output:
        "contam.kmc_pre",
        "contam.kmc_suf"
    conda:
        "../envs/kmc.yaml"
    params:
        k = config["k"]
    shell:
        """
        # create file list
        echo {input} | tr ' ' '\n' > {output.list}

        # create tmp dir
        mkdir kmc_tmp_dir

        # build kmc db
        kmc -k{params.k} -fq @{output.list} contam kmc_tmp_dir/

        # clean up
        rm -r kmc_tmp_dir/
        """
