rule vg_surject:
    input:
        gbz = "../config/pangenomes/{panref}.giraffe.gbz",
        gam = "vg_filter_results/{ID}_{panref}.gam",
        paths = "../config/pangenomes/paths/{linref}.txt"
    output:
        temp("vg_surject_results/{ID}_{panref}_{linref}.bam")
    conda:
        "../envs/vg.yaml"
    shell:
        "vg surject -x {input.gbz} -F {input.paths} -b {input.gam} > {output}"
