rule vg_surject:
    input:
        xg = "{ref}.xg",
        gam = "{ID}_{ref}_filtered.gam"
    output:
        "vg_results/{ID}_{ref}.bam"
    conda:
        "../envs/vg.yaml"
    shell:
        "vg surject -x {input.xg} -b {input.gam} > {output}"
