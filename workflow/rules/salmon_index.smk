rule salmon_index:
    input:
        cds = "degenotate_results/{ref}/cds-nt-longest.fa",
        genome = "../config/linear_genomes/sequence/{ref}.fa",
    output:
        decoy_list = temp("seqkit_results/{ref}_decoys.txt"),
        decoy = temp("{ref}_decoys.fa"),
        done = temp(touch("salmon_index_{ref}.done"))
    params:
        prefix = "salmon_index_{ref}",
        k = config["k"]
    conda:
        "../envs/salmon.yaml"
    shell:
        """
        # create decoy transcriptome file
        cat {input.cds} {input.genome} > {output.decoy}

        # make decoy list
        seqkit seq --name {input.genome} > {output.decoy_list}

        # make salmon index
        salmon index --threads {threads} -t {output.decoy} -i {params.prefix} --decoys {output.decoy_list} -k {params.k}
        """
