rule vg_stats:
    input:
        "vg_filter_results/{ID}_{ref}.gam"
    output:
        "vg_stats_results/{ID}_{ref}.txt"
    conda:
        "../envs/vg.yaml"
    shell:
        """
        vg stats -a {input}  {output}
        """
