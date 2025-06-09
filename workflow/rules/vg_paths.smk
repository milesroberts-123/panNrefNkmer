rule vg_paths:
    input:
        "../config/pangenomes/{ref}.giraffe.gbz",
    output:
        "vg_paths_results/{ref}.fa"
    conda:
        "../envs/vg.yaml"
    shell:
        "vg paths -F -x {input} > {output}"
