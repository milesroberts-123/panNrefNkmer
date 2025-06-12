rule split_sites:
    input:
        "degenotate_results/{ref}/degeneracy-all-sites.bed"
    output:
        tmp=temp("{ref}_4.bed"),
        split_files = expand("split_sites_{{ref}}_{split}", split = range(10, 10 + config["splits"]))
    params:
        prefix = "split_sites_{ref}_",
        splits = config["splits"]
    shell:
        """
        awk '(($5 == 4 || $5 == 0))' {input} > {output.tmp}

        split -n l/{params.splits} --numeric-suffixes=10 {output.tmp} {params.prefix}
        """
