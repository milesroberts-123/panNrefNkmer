rule kmc_contam_db:
    input:
        #expand("../config/contaminants/{contam}.fasta.gz" , contam = config["contams"])
        "contam.fa"
    output:
        "contam.kmc_pre",
        "contam.kmc_suf"
    conda:
        "../envs/kmc.yaml"
    params:
        k = config["k"]
    shell:
        """
        # create tmp dir
        mkdir kmc_tmp_dir

        # build kmc db
        kmc -k{params.k} -fm {input} contam kmc_tmp_dir/

        # clean up
        rm -r kmc_tmp_dir/
        """
