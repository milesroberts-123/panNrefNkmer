# https://unix.stackexchange.com/questions/16443/combine-text-files-column-wise

rule cbind:
    input:
        expand("cbf_results/{ID}.txt", ID=lookup(query="Species == '{species}'", within = reads, cols="BioSample"))
    output:
        "cbf_table_{species}.txt"
    shell:
        r"""
        paste -d' ' {input} > {output}
        """
