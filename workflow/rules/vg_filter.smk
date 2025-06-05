rule vg_filter:
    input:
        gbz = "../config/pangenomes/{ref}.giraffe.gbz",
        gam = "vg_giraffe_results/{ID}_{ref}.gam"
    output:
        "vg_filter_results/{ID}_{ref}.gam"
    conda:
        "../envs/vg.yaml"
    shell:
        "vg filter {input.gam} -r 0.90 -fu -m 1 -q 15 -D 999 -x {input.gbz} > {output}"
