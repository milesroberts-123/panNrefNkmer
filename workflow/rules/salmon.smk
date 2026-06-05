rule salmon_quantmerge:
    input:
        expand("salmon_quant_results/{{ref}}_{ID}/quant.sf", ID = config["rna"])
    output:
        "salmon_quantmerge_results/{ref}.txt"
    params:
        prefix = expand("salmon_quant_results/{{ref}}_{ID}", ID = config["rna"]),
        samnames = config["rna"]
    conda:
        "../envs/salmon.yaml"
    shell:
        """
        salmon quantmerge --quants {params.prefix} --names {params.samnames} -o {output}
        """

rule salmon_quant:
    input:
        idx = "salmon_index_{ref}.done",
        read1 = "fastp_results/trimmed_paired_R1_{ID}.fastq.gz",
        read2 = "fastp_results/trimmed_paired_R2_{ID}.fastq.gz"
    output:
        temp("salmon_quant_results/{ref}_{ID}/quant.sf")
    params:
        prefix = "salmon_quant_results/{ref}_{ID}",
        index = "salmon_index_{ref}"
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

#rule alevin:
#    input:
#        cb =,
#        reads =,
#        tgmap =,
#        index =
#    output:
#
#    params:
#        index = "",
#        outdir = ""
#    conda:
#        "../envs/salmon.yaml"
#    shell:
#        "salmon alevin -l ISR -1 cb.fastq.gz -2 reads.fastq.gz --chromium  -i salmon_index_directory -p {threads} -o alevin_output --tgMap txp2gene.tsv"
