rule vg_index:
    input:
        "../config/pangenomes/{ref}.gfa"
    output:
        dist = "../config/pangenomes/{ref}.dist",
        gbz = "../config/pangenomes/{ref}.giraffe.gbz",
        min = "../config/pangenomes/{ref}.shortread.withzip.min",
        zip = "../config/pangenomes/{ref}.shortread.zipcodes"
    log:
        "logs/vg_index/{ref}.log"
    conda:
        "../envs/vg.yaml"
    shell:
        "vg autoindex -w giraffe -g {input} -p {wildcards.ref}"
