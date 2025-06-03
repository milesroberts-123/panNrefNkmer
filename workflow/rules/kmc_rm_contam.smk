rule kmc_contam_db:
    input:
        "raw_reads/{ID}_1.fastq.gz",
        "raw_reads/{ID}_2.fastq.gz",
        "contam.kmc_pre",
        "contam.kmc_suf"
    output:
        "no_contam_reads/{ID}_1.fastq.gz",
        "no_contam_reads/{ID}_2.fastq.gz"
    conda:
        "../envs/kmc.yaml"
    params:
        k = config["k"]
    shell:
        """
        kmc_tools contam kmc_db @input_files.txt -cx50 filtered.fastq

        # create file list
        echo {input} | tr ' ' '\n' > {output.list}

        # create tmp dir
        mkdir kmc_tmp_dir

        # build kmc db
        kmc -k{params.k} -fq @{output.list} contam kmc_tmp_dir/

        # clean up
        rm -r kmc_tmp_dir/
        """
