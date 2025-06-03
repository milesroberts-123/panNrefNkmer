rule vg_index:
    input:
        "{ref}.vg"
    output:
        gcsa = "{ref}.gcsa",
        xg = "{ref}.xg"
    params:
        k = config["k"]
    log:
        "logs/vg_index/{ref}.log"
    conda:
        "../envs/vg.yaml"
    shell:
        "vg index -x {output.xg} -g {output.gcsa} -k {params.k} {input} &> {log}"
