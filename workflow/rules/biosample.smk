rule biosample:
    input:
        pread1=expand("fastp_results/trimmed_paired_R1_{run}.fastq.gz", run = lookup(query="BioSample == '{bio}'", within=reads, cols="Run")),
        pread2=expand("fastp_results/trimmed_paired_R2_{run}.fastq.gz", run = lookup(query="BioSample == '{bio}'", within=reads, cols="Run")),
        uread1=expand("fastp_results/trimmed_unpaired_R1_{run}.fastq.gz", run = lookup(query="BioSample == '{bio}'", within=reads, cols="Run")),
        uread2=expand("fastp_results/trimmed_unpaired_R2_{run}.fastq.gz", run = lookup(query="BioSample == '{bio}'", within=reads, cols="Run")),
    output:
        pread1="biosample_results/{bio}_paired_R1.fastq.gz",
        pread2="biosample_results/{bio}_paired_R2.fastq.gz",
        uread1="biosample_results/{bio}_unpaired_R1.fastq.gz",
        uread2="biosample_results/{bio}_unpaired_R2.fastq.gz",
    shell:
        """
        zcat {input.pread1} | gzip > {output.pread1}
        zcat {input.pread2} | gzip > {output.pread2}
        zcat {input.uread1} | gzip > {output.uread1}
        zcat {input.uread2} | gzip > {output.uread2}
        """
