rule split_sites:
    input:
        "degenotate_results/{ref}/degeneracy-all-sites.bed"
    output:
        "split_sites_results/{ref}_{chrom}.bed"
    shell:
        """
        awk '(($1 == "{wildcards.chrom}"))' {input} > {output}
        """
