# Need to make sure paired reads have the same names so that BWA works
# Sometimes reads on SRA don't, even though they are obviously paired
# https://www.biostars.org/p/68477/
rule change_headers:
    input:
        pread1="results/fastp/trimmed_paired_R1_{ID}.fastq.gz",
        pread2="results/fastp/trimmed_paired_R2_{ID}.fastq.gz",
        uread1="results/fastp/trimmed_unpaired_R1_{ID}.fastq.gz",
        uread2="results/fastp/trimmed_unpaired_R2_{ID}.fastq.gz",
    output:
        pread1=temp("results/change_headers/paired_R1_{ID}.fastq.gz"),
        pread2=temp("results/change_headers/paired_R2_{ID}.fastq.gz"),
        cat=temp("results/change_headers/unpaired_{ID}.fastq.gz")
    shell:
        """
        if [ ! -d "results/change_headers" ]; then
            mkdir -p results/change_headers
        fi

        # change headers
        zcat {input.pread1} |  awk '{{if(NR%4==1) $0=sprintf("@1_%d",(1+i++)); print;}}' | gzip -c > {output.pread1}

        zcat {input.pread2} |  awk '{{if(NR%4==1) $0=sprintf("@1_%d",(1+i++)); print;}}' | gzip -c > {output.pread2}

        # combine unpaired reads into one file
        zcat {input.uread1} {input.uread2} | gzip -c > {output.cat}
        """
