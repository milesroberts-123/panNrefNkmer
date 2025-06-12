rule pixy:
    input:
        vcf = "bcftools_concat_results/{ref}.vcf",
        gff = "../config/linear_genomes/annotation/{ref}.gff",
    output:
        bed = temp("pixy_results/{ref}_genes.bed"),
        pixy = "pixy_results/{ref}.txt"
    conda:
        "../envs/pixy.yaml"
    shell:
        """
        # get bed intervals by gene
        awk '(($3 == "gene"))' {input.gff} > {output.bed}
        
        # calculate statistics by gene
        pixy --vcf {input.vcf} --windows {output.bed} --stats pi watterson_theta tajima_d --output_folder pixy_results --output_prefix {wildcards.ref}_
        """
