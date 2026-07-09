rule multiqc:
    input:
        expand("results/fastp/{ID}.json", ID = reads["Run"].tolist()),
        expand("results/fastp/no_contam_{ID}.json", ID = reads[reads["Group"] == "ingroup"]["Run"].tolist())
    output:
        "results/multiqc/multiqc_report.html"
    conda:
        "../envs/multiqc.yaml"
    shell:
        "multiqc results/fastp/ -o results/multiqc"
