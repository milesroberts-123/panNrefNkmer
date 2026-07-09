rule split_sites:
    input:
        "results/degenotate/{ref}/degeneracy-all-sites.bed"
    output:
        tmp=temp("results/split/{ref}_04.bed"),
        split_files = temp(expand("results/split/split_sites_{{ref}}_{split}", split = range(10, 10 + config["splits"])))
    params:
        prefix = "results/split/split_sites_{ref}_",
        splits = config["splits"]
    shell:
        """
        awk '(($5 == 4 || $5 == 0))' {input} > {output.tmp}

        split -n l/{params.splits} --numeric-suffixes=10 {output.tmp} {params.prefix}
        """
