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
