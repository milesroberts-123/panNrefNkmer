rule salmon_quantmerge:
    input:
        expand("results/salmon_quant/{{ref}}_{ID}/quant.sf", ID = config["rna"])
    output:
        "results/salmon_quantmerge/{ref}.txt"
    params:
        samnames = config["rna"]
    conda:
        "../envs/salmon.yaml"
    shell:
        """
        salmon quantmerge --quants {input} --names {params.samnames} -o {output}
        """

rule salmon_quant:
    input:
        idx = "results/salmon_index_{ref}.done",
        read1 = "results/fastp/trimmed_paired_R1_{ID}.fastq.gz",
        read2 = "results/fastp/trimmed_paired_R2_{ID}.fastq.gz"
    output:
        temp("results/salmon_quant/{ref}_{ID}/quant.sf")
    params:
        prefix = "results/salmon_quant/{ref}_{ID}",
        index = "results/salmon_index_{ref}"
    conda:
        "../envs/salmon.yaml"
    shell:
        """
        salmon quant -i {params.index} \
            -l A \
            -1 {input.read1} \
            -2 {input.read2} \
            -p {threads} \
            --validateMappings \
            -o {params.prefix}
        """

rule salmon_index:
    input:
        cds = "results/degenotate/{ref}/cds-nt-longest.fa",
        genome = "../config/linear_genomes/sequence/{ref}.fa",
    output:
        decoy_list = temp("results/seqkit/{ref}_decoys.txt"),
        decoy = temp("results/{ref}_decoys.fa"),
        done = temp(touch("results/salmon_index_{ref}.done"))
    params:
        prefix = "results/salmon_index_{ref}",
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
