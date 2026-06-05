rule jules_bwa_index:
    input:
        config["reference_genome_path"] + "{ref}/{ref}.fasta",
    output:
        amb=temp(config["reference_genome_path"] + "{ref}/{ref}.fasta.amb"),
        ann=temp(config["reference_genome_path"] + "{ref}/{ref}.fasta.ann"),
        bwt=temp(config["reference_genome_path"] + "{ref}/{ref}.fasta.bwt"),
        pac=temp(config["reference_genome_path"] + "{ref}/{ref}.fasta.pac"),
        sa=temp(config["reference_genome_path"] + "{ref}/{ref}.fasta.sa"),
    conda:
        "../envs/bwa.yaml"
    shell:
        """
        bwa index {input}
        """

rule jules_samtools_faidx:
    input:
        config["reference_genome_path"] + "{ref}/{ref}.fasta"    
    output:
        temp(config["reference_genome_path"] + "{ref}/{ref}.fasta.fai")
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        samtools faidx {ref}
        """

rule jules_fastq_dump:
    output:
        r1=temp("raw_reads/{ID}_1.fastq.gz"),
        r2=temp("raw_reads/{ID}_2.fastq.gz")
    conda:
        "../envs/sra.yaml"
    shell:
        """
        # create directory
        if [ ! -d "raw_reads" ]; then
            mkdir raw_reads
        fi

        prefetch {wildcards.ID}

        # download data
        fastq-dump --progress --temp /tmp --gzip --outdir ./raw_reads --split-files --skip-technical ./{wildcards.ID}
        """

rule jules_trimmomatic:
    input:
        r1="raw_reads/{srr}_1.fastq.gz",
        r2="raw_reads/{srr}_2.fastq.gz",
        adapters=config["adapter_path"]
    output:
        r1=temp("trim_reads/{srr}_r1.fastq.gz"),
        r2=temp("trim_reads/{srr}_r2.fastq.gz"),
        u1=temp("trim_reads/{srr}_u1.fastq.gz"),
        u2=temp("trim_reads/{srr}_u2.fastq.gz")
    conda:
        "../envs/trimmomatic.yaml"
    shell:
        """
        trimmomatic PE -threads {threads} {input.r1} {input.r2} {output.r1} {output.u1} {output.r2} {output.u2} ILLUMINACLIP:{input.adapters}:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36
        """

rule jules_bwa_mem:
    input:
        r1="trim_reads/{srr}_r1.fastq.gz",
        r2="trim_reads/{srr}_r2.fastq.gz",
        ref=expand("{path_start}{ref}/{ref}.fasta", path_start = config["reference_genome_path"], ref = lookup(query = "Run == '{srr}'", within = reads, cols = "Species"))
    output:
        temp("align_reads/{srr}.bam")
    conda:
        "../envs/bwa.yaml"
    shell:
        """
        bwa mem -t {threads} -M -R '@RG\\tID:{wildcards.srr}\\tSM:{wildcards.srr}\\tPL:ILLUMINA\\tLB:lib1' {input.ref} {input.r1} {input.r2} | samtools sort -@ {threads} -o {output}
        """

rule jules_picard_mark_dup:
    input:
        "align_reads/{srr}.bam",
    output:
        bam="mark_reads/{srr}.bam",
        metrics=temp("picard_metrics/{srr}.txt"),
        tmp_dir=temp(directory("tmp_picard_{srr}"))
    conda:
        "../envs/picard.yaml"
    shell:
        """
        picard MarkDuplicates I={input} O={output.bam} M={output.metrics} VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true TMP_DIR={output.tmp_dir} 
        """

rule jules_samtools_index:
    input:
        "mark_reads/{srr}.bam"
    output:
        "mark_reads/{srr}.bam.bai"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        samtools index -@ {threads} {input}
        """

rule jules_sample_coverage:
    input:
        "mark_reads/{srr}.bam"
    output:
        "coverages/{srr}.50k.coverage.txt"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        samtools depth -a {input} | awk '{{sum+=$3}} END {{print "Average =", sum/NR}}' > {output}
        """

rule jules_bcftools_mpileup:
    input:
        ref=expand("{path_start}{ref}/{ref}.fasta", path_start = config["reference_genome_path"], ref = lookup(query = "Run == '{srr}'", within = reads, cols = "Species")),
        bam="mark_reads/{srr}.bam"
    output:
        vcf="bcfvcfs/{srr}_raw.vcf.gz",
        tbi="bcfvcfs/{srr}_raw.vcf.gz.tbi"
    conda:
        "../envs/bcftools.yaml"
    params:
        q=config["q"],
        Q=config["Q"]
    shell:
        """
        bcftools mpileup -f {input.ref} -a "FORMAT/AD,FORMAT/DP" -q {params.q} -Q {params.Q} {input.bam} | bcftools call -mv -Oz -o {output.vcf}

        bcftools index {output.vcf}
        """

rule jules_bcftools_filter:
    input:
        "bcfvcfs/{srr}_raw.vcf.gz"
    output:
        "bcfvcfs/{srr}_filtered.vcf.gz"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        bcftools filter -i 'QUAL>20 && FORMAT/DP>=4 && INFO/DP>=4 && INFO/MQ>30' {input} | bcftools view -Oz -o {output}
        """

rule jules_bcftools_roh:
    input:
        "bcfvcfs/{srr}_filtered.vcf.gz"
    output:
        "roh/{srr}_ROH.txt"
    conda:
        "../envs/bcftools.yaml"
    params:
        G=config["G"],
        AFdflt=config["AFdflt"]
    shell:
        """
        bcftools roh -G{params.G} --AF-dflt {params.AFdflt} -o {output} {input}
        """

rule jules_psmc_50k_bed:
    input: 
        expand("{path_start}{{ref}}/{{ref}}.fasta.fai", path_start = config["reference_genome_path"])
    output:
        "{ref}.50k.bed"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        zless {input} | awk '$2>50000 {{print $1, "0", $2}}' > {output}
        """

rule jules_psmc_subset_bam:
    input: 
        bam="mark_reads/{srr}.bam",
        bed=expand("{species}.50k.bed", species=lookup(query="Run == '{srr}'", within=reads, cols="Species"))
    output:
        "psmc/{srr}.50k.bam"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools view -@ {threads} -bh -L {input.bed} \
        -o {output} {input.bam}
        """

rule jules_psmc_gen_consensus: 
    input:
        ref=expand("{path_start}{ref}/{ref}.fasta", path_start = config["reference_genome_path"], ref = lookup(query = "Run == '{srr}'", within = reads, cols = "Species")),
        bam="psmc/{srr}.50k.bam",
        cov="coverages/{srr}.50k.coverage.txt"
    output: 
        "psmc/{srr}.con.fq.gz"
    conda:  
        "../envs/psmc_legacy.yaml"
    params:
        mincov=lambda wc: int(open("coverages/{}.50k.coverage.txt".format(wc.srr)).read().split()[2]) // 3,
        maxcov=lambda wc: int(open("coverages/{}.50k.coverage.txt".format(wc.srr)).read().split()[2]) * 2
    shell:
        """
        samtools mpileup -C50 -uf {input.ref} {input.bam} | \
            bcftools call -c - | \
            vcfutils.pl vcf2fq -d {params.mincov} -D {params.maxcov} | \
            gzip > {output}
        """

rule jules_psmc_gen_input:
    input: 
        "psmc/{srr}.con.fq.gz"
    output:
        "psmc/{srr}.psmcfa"
    conda: 
        "../envs/psmc_legacy.yaml"
    params:
        psmc_path=config["psmc_path"]
    shell:
        """
        {params.psmc_path}/fq2psmcfa -q20 {input} \
        > {output}
        """

rule jules_psmc_run_psmc: 
    input: 
         "psmc/{srr}.psmcfa"
    output: 
         "psmc/{srr}.psmc"
    conda:
         "../envs/psmc_legacy.yaml"
    params:
        psmc_path=config["psmc_path"]
    shell: 
        """
        {params.psmc_path}/psmc -N25 -t15 -r5 -p "1+1+1+1+25*2+4+6" \
        -o {output} {input}
        """
