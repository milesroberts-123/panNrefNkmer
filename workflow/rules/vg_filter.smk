rule vg_filter:
    input:
        xg = "{ref}.xg",
        gam = "{ID}_{ref}.gam"
    output:
        "{ID}_{ref}_filtered.gam"
    conda:
        "../envs/vg.yaml"
    shell:
        "vg filter {input.gam} -r 0.90 -fu -m 1 -q 15 -D 999 -x {input.xg} > {output}"
