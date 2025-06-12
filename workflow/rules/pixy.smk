rule pixy:
    input:
        vcf = "bcftools_concat_results/{ref}_sorted.vcf.gz",
        tbi = "bcftools_concat_results/{ref}_sorted.vcf.gz.tbi",
        gff = "../config/linear_genomes/annotation/{ref}.gff",
    output:
        bed = temp("pixy_results/{ref}_genes.bed"),
        ingroup = "ingroup_{ref}.txt",
        outgroup = "outgroup_{ref}.txt",
        populations = "populations_{ref}.txt",
        pixy = "pixy_results/{ref}.txt"
    params:
        ingroup = config["samples"],
        outgroup = config["outgroup"]
    conda:
        "../envs/pixy.yaml"
    shell:
        r"""
        # get bed intervals by gene
        awk '(($3 == "gene"))' {input.gff} > {output.bed}
        
        echo {params.ingroup} | sed 's: :\n:g' | sed 's:$:\tingroup:g' > {output.ingroup}
        echo {params.outgroup} | sed 's: :\n:g' | sed 's:$:\toutgroup:g' > {output.outgroup}

        cat {output.ingroup} {output.outgroup} > {output.populations}

        # calculate statistics by gene
        pixy --populations {output.populations} --vcf {input.vcf} --bed_file {output.bed} --stats pi fst dxy watterson_theta tajima_d --output_folder pixy_results --output_prefix {wildcards.ref}_
        """
