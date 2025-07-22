rule alevin:
    input:
        cb =,
        reads =,
        tgmap =,
        index =
    output:

    params:
        index = "",
        outdir = ""
    conda:
        "../envs/salmon.yaml"
    shell:
        "salmon alevin -l ISR -1 cb.fastq.gz -2 reads.fastq.gz --chromium  -i salmon_index_directory -p {threads} -o alevin_output --tgMap txp2gene.tsv"
