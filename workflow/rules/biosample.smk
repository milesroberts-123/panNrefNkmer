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
    params:
        no_inputs=lambda wildcards, input: len(input)
    shell:
        """
        if [ {params.no_inputs} -gt 4 ]; then
            zcat {input.pread1} | gzip > {output.pread1}
            zcat {input.pread2} | gzip > {output.pread2}
            zcat {input.uread1} | gzip > {output.uread1}
            zcat {input.uread2} | gzip > {output.uread2}
        else
            mv {input.pread1} {output.pread1}
            mv {input.pread2} {output.pread2}
            mv {input.uread1} {output.uread1}
            mv {input.uread2} {output.uread2}
        fi
        """
