# https://unix.stackexchange.com/questions/16443/combine-text-files-column-wise

rule cbind:
    input:
        expand("cbf_results/{ID}.txt", ID=config["samples"])
    output:
        "cbf_table.txt"
    shell:
        "paste {input} | column -s $'\t' -t > {output}"
