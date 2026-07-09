rule divergence_estimation:
    input:
        "results/separate_chrom/{chr}.fa.gz",
        "results/separate_chrom/{chr}.fa.gz.tbi"
    output:
        "results/mash/{chr}_triangle.txt"
    params:
        k = config["k"]
    conda:
        "../envs/mash.yaml"
    shell:
        "mash triangle -k {params.k} {input} > {output}"
