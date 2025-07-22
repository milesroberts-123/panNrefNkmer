rule divergence_estimation:
    input:
        "{chr}.fa.gz",
        "{chr}.fa.gz.tbi"
    output:
        "{chr}_triangle.txt"
    params:
        k = config["k"]
    conda:
        "../envs/mash.yaml"
    shell:
        "mash triangle -k {params.k} {input} > {output}"
