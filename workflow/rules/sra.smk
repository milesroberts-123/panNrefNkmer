rule sra:
    output:
        r1=temp("raw_reads/{ID}_1.fastq"),
        r2=temp("raw_reads/{ID}_2.fastq")
    conda:
        "../envs/sra.yaml"
    shell:
        """
        # create directory
        if [ ! -d "raw_reads" ]; then
            mkdir raw_reads
        fi

        # download data
        fasterq-dump --threads {threads} --split-files --skip-technical {wildcards.ID} 
        
        # move data to folder
        mv {wildcards.ID}_1.fastq raw_reads/
        mv {wildcards.ID}_2.fastq raw_reads/
        """
