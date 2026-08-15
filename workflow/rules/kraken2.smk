rule kraken2_rm_microbial:
    input:
        pread1="results/biosample/{ID}_paired_R1.fastq.gz",
        pread2="results/biosample/{ID}_paired_R2.fastq.gz"
    output:
        unclass1=temp("results/kraken2_unclassified/{ID}_unclassified_1.fastq"),
        unclass2=temp("results/kraken2_unclassified/{ID}_unclassified_2.fastq"),
        report="results/kraken2_reports/{ID}.txt"
    params:
        db=config["kraken_db_path"],
        confidence=config["kraken_confidence"]
    conda:
        "../envs/kraken2.yaml"
    shell:
        """
        mkdir -p results/kraken2_unclassified results/kraken2_reports

        kraken2 --db {params.db} --paired --threads {threads} \
            --confidence {params.confidence} \
            --unclassified-out results/kraken2_unclassified/{wildcards.ID}_unclassified#.fastq \
            --report {output.report} \
            {input.pread1} {input.pread2}
        """
