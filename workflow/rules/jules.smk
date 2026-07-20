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
        samtools faidx {input}
        """

rule jules_fastq_dump:
    output:
        r1=temp("results/raw_reads/{ID}_1.fastq.gz"),
        r2=temp("results/raw_reads/{ID}_2.fastq.gz")
    conda:
        "../envs/sra.yaml"
    shell:
        """
        # create directory
        if [ ! -d "results/raw_reads" ]; then
            mkdir -p results/raw_reads
        fi

        prefetch {wildcards.ID}

        # download data
        fastq-dump --progress --temp /tmp --gzip --outdir ./results/raw_reads --split-files --skip-technical ./{wildcards.ID}
        """

rule jules_trimmomatic:
    input:
        r1="results/raw_reads/{srr}_1.fastq.gz",
        r2="results/raw_reads/{srr}_2.fastq.gz",
        adapters=config["adapter_path"]
    output:
        r1=temp("results/trim_reads/{srr}_r1.fastq.gz"),
        r2=temp("results/trim_reads/{srr}_r2.fastq.gz"),
        u1=temp("results/trim_reads/{srr}_u1.fastq.gz"),
        u2=temp("results/trim_reads/{srr}_u2.fastq.gz")
    conda:
        "../envs/trimmomatic.yaml"
    shell:
        """
        trimmomatic PE -threads {threads} {input.r1} {input.r2} {output.r1} {output.u1} {output.r2} {output.u2} ILLUMINACLIP:{input.adapters}:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36
        """

rule jules_bwa_mem:
    input:
        r1="results/trim_reads/{srr}_r1.fastq.gz",
        r2="results/trim_reads/{srr}_r2.fastq.gz",
        ref=expand("{path_start}{ref}/{ref}.fasta", path_start = config["reference_genome_path"], ref = lookup(query = "Run == '{srr}'", within = reads, cols = "Species"))
    output:
        temp("results/align_reads/{srr}.bam")
    conda:
        "../envs/bwa.yaml"
    shell:
        """
        bwa mem -t {threads} -M -R '@RG\\tID:{wildcards.srr}\\tSM:{wildcards.srr}\\tPL:ILLUMINA\\tLB:lib1' {input.ref} {input.r1} {input.r2} | samtools sort -@ {threads} -o {output}
        """

rule jules_picard_mark_dup:
    input:
        "results/align_reads/{srr}.bam",
    output:
        bam="results/mark_reads/{srr}.bam",
        metrics=temp("results/picard_metrics/{srr}.txt"),
        tmp_dir=temp(directory("results/tmp_picard_{srr}"))
    conda:
        "../envs/picard.yaml"
    shell:
        """
        picard MarkDuplicates I={input} O={output.bam} M={output.metrics} VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true TMP_DIR={output.tmp_dir} 
        """

rule jules_samtools_index:
    input:
        "results/mark_reads/{srr}.bam"
    output:
        "results/mark_reads/{srr}.bam.bai"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        samtools index -@ {threads} {input}
        """

rule jules_sample_coverage:
    input:
        "results/mark_reads/{srr}.bam"
    output:
        "results/coverages/{srr}.50k.coverage.txt"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        samtools depth -a {input} | awk '{{sum+=$3}} END {{print "Average =", sum/NR}}' > {output}
        """

rule jules_bcftools_mpileup:
    input:
        ref=expand("{path_start}{ref}/{ref}.fasta", path_start = config["reference_genome_path"], ref = lookup(query = "Run == '{srr}'", within = reads, cols = "Species")),
        bam="results/mark_reads/{srr}.bam"
    output:
        vcf="results/bcfvcfs/{srr}_raw.vcf.gz",
        tbi="results/bcfvcfs/{srr}_raw.vcf.gz.tbi"
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
        vcf="results/bcfvcfs/{srr}_raw.vcf.gz",
        cov="results/coverages/{srr}.50k.coverage.txt"
    output:
        "results/bcfvcfs/{srr}_filtered.vcf.gz"
    conda:
        "../envs/bcftools.yaml"
    params:
        mincov=lambda wc: int(float(open("coverages/{}.50k.coverage.txt".format(wc.srr)).read().split()[2])) // 3,
        maxcov=lambda wc: int(float(open("coverages/{}.50k.coverage.txt".format(wc.srr)).read().split()[2])) * 2
    shell:
        """
        bcftools filter -i 'QUAL>=30 && FORMAT/DP>={params.mincov} && FORMAT/DP<={params.maxcov} && INFO/DP>={params.mincov} && INFO/MQ>=30 && FORMAT/SP<60' {input.vcf} | bcftools view -v snps -m2 -M2 -Oz -o {output}
        """

rule jules_bcftools_roh:
    input:
        "results/bcfvcfs/{srr}_filtered.vcf.gz"
    output:
        "results/roh/{srr}_ROH.txt"
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
        "results/psmc_bed/{ref}.50k.bed"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        cat {input} | awk '$2>50000 {{print $1, "0", $2}}' > {output}
        """

rule jules_psmc_subset_bam:
    input: 
        bam="results/mark_reads/{srr}.bam",
        bed=expand("results/psmc_bed/{species}.50k.bed", species=lookup(query="Run == '{srr}'", within=reads, cols="Species"))
    output:
        "results/psmc/{srr}.50k.bam"
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
        bam="results/psmc/{srr}.50k.bam",
        cov="results/coverages/{srr}.50k.coverage.txt"
    output: 
        "results/psmc/{srr}.con.fq.gz"
    conda:  
        "../envs/psmc_legacy.yaml"
    params:
        mincov=lambda wc: int(open("results/coverages/{}.50k.coverage.txt".format(wc.srr)).read().split()[2]) // 3,
        maxcov=lambda wc: int(open("results/coverages/{}.50k.coverage.txt".format(wc.srr)).read().split()[2]) * 2
    shell:
        """
        samtools mpileup -C50 -uf {input.ref} {input.bam} | \
            bcftools call -c - | \
            vcfutils.pl vcf2fq -d {params.mincov} -D {params.maxcov} | \
            gzip > {output}
        """

rule jules_psmc_gen_consensus:
    input:
        ref=expand("{path_start}{ref}/{ref}.fasta", path_start = config["reference_genome_path"], ref = lookup(query = "Run == '{srr}'", within = reads, cols = "Species")),
        bam="results/psmc/{srr}.50k.bam",
        cov="results/coverages/{srr}.50k.coverage.txt"
    output:
        "results/psmc/{srr}.con.fq.gz"
    conda:
        "../envs/psmc_legacy.yaml"
    params:
        mincov=lambda wc: int(float(open("results/coverages/{}.50k.coverage.txt".format(wc.srr)).read().split()[2])) // 3,
        maxcov=lambda wc: int(float(open("results/coverages/{}.50k.coverage.txt".format(wc.srr)).read().split()[2])) * 2
    shell:
        """
        samtools mpileup -C50 -uf {input.ref} {input.bam} | \
            bcftools call -c - | \
            vcfutils.pl vcf2fq -d {params.mincov} -D {params.maxcov} | \
            gzip > {output}
        """

rule jules_psmc_gen_input:
    input: 
        "results/psmc/{srr}.con.fq.gz"
    output:
        "results/psmc/{srr}.psmcfa"
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
         "results/psmc/{srr}.psmcfa"
    output: 
         "results/psmc/{srr}.psmc"
    conda:
         "../envs/psmc_legacy.yaml"
    params:
        psmc_path=config["psmc_path"]
    shell: 
        """
        {params.psmc_path}/psmc -N25 -t15 -r5 -p "1+1+1+1+25*2+4+6" \
        -o {output} {input}
        """
