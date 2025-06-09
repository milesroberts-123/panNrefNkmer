rule datasets:
    output:
        "contam.fa"
    params:
        contams = config["ncbi_contams"]
    conda:
        "../envs/datasets.yaml"
    shell:
        r"""
        # get all contaminating genomes
        datasets download genome taxon {params.contams} --reference --dehydrated --filename contam.zip

        # unpack metadata
        unzip contam.zip -d contam_genomes

        # download genomes based on metadata
        datasets rehydrate --directory contam_genomes/

        # search directory for all genomes and copy them
        find contam_genomes/ncbi_dataset/data/ -type f -name '*.fna' -exec cat {{}} \; > {output}

        # clean up
        #rm -r contam_genomes/
        """
