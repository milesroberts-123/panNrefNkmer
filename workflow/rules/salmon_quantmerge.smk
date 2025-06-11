rule salmon_quantmerge:
    input:
        expand("qaunts/{ID}_{{ref}}_quant.sf", ID = config["rna"])
    output:
        "{ref}_salmon_quant.txt"
    params:
        prefix = expand("qaunts/{ID}_{{ref}}_quant", ID = config["rna"]),
    conda:
        "../envs/salmon.yaml"
    shell:
        """
        salmon quantmerge --quants {params.prefix} -c tpm -o {output}
        """
