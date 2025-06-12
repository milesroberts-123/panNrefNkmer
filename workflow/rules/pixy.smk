rule pixy:
    input:
        vcf = "bcftools_concat_results/{ref}_sorted.vcf.gz",
        tbi = "bcftools_concat_results/{ref}_sorted.vcf.gz.tbi",
        gff = "../config/linear_genomes/annotation/{ref}.gff",
    output:
        bed = temp("pixy_results/{ref}_genes.bed"),
        populations = "populations_{ref}.txt",
        pixy = "pixy_results/{ref}.txt"
    params:
        ingroup = config["samples"],
    conda:
        "../envs/pixy.yaml"
    shell:
        r"""
        # get bed intervals by gene
        awk '(($3 == "gene"))' {input.gff} | cut -f1,4,5 > {output.bed}
        
        echo {params.ingroup} | sed 's: :\n:g' | sed 's:$:\tingroup:g' > {output.populations}

        # calculate statistics by gene
        pixy --populations {output.populations} --vcf {input.vcf} --bed_file {output.bed} --stats pi fst dxy watterson_theta tajima_d --output_folder pixy_results --output_prefix {wildcards.ref}_
        """
