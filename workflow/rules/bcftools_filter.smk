rule bcftools_filter:
    input:
        "bcftools_linref_merge_results/{ref}_{split}.vcf.gz"
    output:
        invar_before = temp("bcftools_filter_results/{ref}_{split}_invar.vcf.gz"),
        var_before = temp("bcftools_filter_results/{ref}_{split}_var.vcf.gz"),
        invar_after = temp("bcftools_filter_results/{ref}_{split}_invar_filt.vcf.gz"),
        var_after = temp("bcftools_filter_results/{ref}_{split}_var_filt.vcf.gz"),
        invar_after_tbi = temp("bcftools_filter_results/{ref}_{split}_invar_filt.vcf.gz.tbi"),
        var_after_tbi = temp("bcftools_filter_results/{ref}_{split}_var_filt.vcf.gz.tbi"),
        unsorted = temp("bcftools_filter_results/{ref}_{split}_unsorted.vcf.gz"),
        sorted = temp("bcftools_filter_results/{ref}_{split}_sorted.vcf.gz")
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        # separate invariant vs variant (bi-allelic snp) sites
        echo Separating invariant and variant sites
        bcftools view --max-af 0 -Oz -o {output.invar_before} {input}
        bcftools view --types snps -m 2 -M 2 -Oz -o {output.var_before} {input}

        # filter variant sites
        echo Filtering variant sites...
        bcftools filter -i 'FMT/AD>=5 & FMT/GQ>=20 & F_MISSING<=0.8' -Oz -o {output.var_after} {output.var_before}

        # filter invariant sites
        echo Filtering invariant sites...
        bcftools filter -i 'FMT/AD>=5 & F_MISSING<=0.8' -Oz -o {output.invar_after} {output.invar_before}

        # index filtered sites
        echo Indexing filtered sites...
        tabix {output.invar_after}
        tabix {output.var_after}

        # recombine filtered sites
        echo Recombining sites after filtering...
        bcftools concat --allow-overlaps -Oz -o {output.unsorted} {output.var_after} {output.invar_after}
        bcftools sort -Oz -o {output.sorted} {output.unsorted}
        tabix {output.sorted}
        """
