# https://unix.stackexchange.com/questions/16443/combine-text-files-column-wise

rule cbind:
    input:
        expand("cbf_results/{ID}.txt", ID=config["samples"])
    output:
        "cbf_table.txt"
    shell:
        r"""
        paste -d' ' {input} > {output}
        """
