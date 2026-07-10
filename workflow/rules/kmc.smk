rule kmc:
    input:
        #pread1="results/biosample/{ID}_paired_R1.fastq.gz",
        #pread2="results/biosample/{ID}_paired_R2.fastq.gz",
        #uread1="results/biosample/{ID}_unpaired_R1.fastq.gz",
        #uread2="results/biosample/{ID}_unpaired_R2.fastq.gz"
        "results/no_contam_reads/{ID}_stage1.fastq"
    output:
        counts=temp("results/kmc/{ID}.txt"),
        tmp_pre=temp("results/kmc_db_{ID}.kmc_pre"),
        tmp_suf=temp("results/kmc_db_{ID}.kmc_suf"),
        #list=temp("results/input_{ID}.txt")
    conda:
        "../envs/kmc.yaml"
    params:
        mincount=config["mincount"],
        maxcount=config["maxcount"],
        k=config["k"],
    shell:
        """
        # create directory
        if [ -d "tmp_kmc_{wildcards.ID}" ]; then
            rm -r tmp_kmc_{wildcards.ID}
        fi

        mkdir tmp_kmc_{wildcards.ID}

        # count k-mers
        kmc -m15 -t{threads} -ci{params.mincount} -cs{params.maxcount} -k{params.k} {input} results/kmc_db_{wildcards.ID} tmp_kmc_{wildcards.ID}

        # dump all k-mers to text file
        kmc_tools transform results/kmc_db_{wildcards.ID} dump {output.counts}

        # delete tmp directories
        rm -r tmp_kmc_{wildcards.ID}
        """

rule kmc_ref_db:
    input:
        "../config/reference_genomes/{species}/{species}.fasta"
    output:
        pre=temp("results/kmc_ref_dbs/{species}.kmc_pre"),
        suf=temp("results/kmc_ref_dbs/{species}.kmc_suf")
    conda:
        "../envs/kmc.yaml"
    params:
        k = config["k"]
    shell:
        """
        if [ ! -d "kmc_tmp_{wildcards.species}" ]; then
            mkdir kmc_tmp_{wildcards.species}/
        fi

        # build kmc db
        kmc -k{params.k} -m23 -t{threads} -ci1 -cs2 -fm {input} results/kmc_ref_dbs/{wildcards.species} kmc_tmp_{wildcards.species}/

        # clean up
        rm -r kmc_tmp_{wildcards.species}/
        """

rule kmc_rm_contam:
    input:
        pread1="results/biosample/{ID}_paired_R1.fastq.gz",
        pread2="results/biosample/{ID}_paired_R2.fastq.gz",
        uread1="results/biosample/{ID}_unpaired_R1.fastq.gz",
        uread2="results/biosample/{ID}_unpaired_R2.fastq.gz",
        refdb=expand("results/kmc_ref_dbs/{species}.{ext}", ext=["kmc_pre", "kmc_suf"], species=lookup(query="BioSample == '{ID}'", within = reads, cols="Species")),
    output:
        filt1=temp("results/no_contam_reads/{ID}_stage1.fastq"),
        list=temp("results/input_{ID}.txt")
    params:
        refdb = lambda wildcards, input: os.path.splitext(str(input.refdb[0]))[0],
        contamMatchLimitCount = config["contam_match_limit_count"],
    conda:
        "../envs/kmc.yaml"
    shell:
        """
        if [ ! -d "results/no_contam_reads" ]; then
            mkdir -p results/no_contam_reads/
        fi

        # create file list
        echo {input.pread1} {input.pread2} {input.uread1} {input.uread2} | tr ' ' '\n' > {output.list}

        # filter reads for contamination
        kmc_tools -t{threads} filter {params.refdb} @{output.list} -ci{params.contamMatchLimitCount} -cx1000000 {output.filt1}
        """

rule kmc_combine_dbs:
    input:
        expand("results/kmc_db_{ID}.{suffix}", ID=lookup(query="Species == {species}", within=reads, cols="BioSample"), suffix=["kmc_pre", "kmc_suf"]),
    output:
        expand("results/kmc_combine_dbs/{{species}}.{suffix}", suffix=["kmc_pre", "kmc_suf"]),
        temp("results/kmc_combine_dbs/{species}.complex")
    conda:
        "../envs/kmc.yaml"
    shell:
        """
        mkdir -p results/kmc_combine_dbs

        {{
            echo "INPUT:"
            printf '%s\\n' {input.dbs} | grep '\\.kmc_pre$' | sed 's/\\.kmc_pre$//' | awk '{{print "set" NR " = " $0 " -ci1"}}'
            echo "OUTPUT:"
            printf "results/kmc_combine_dbs/{wildcards.species} = "
            printf '%s\\n' {input.dbs} | grep '\\.kmc_pre$' | sed 's/\\.kmc_pre$//' | awk '{{printf "%sset%d", (NR>1?" + ":""), NR}} END{{print ""}}'
            echo "OUTPUT_PARAMS:"
            echo "-cs10000000000"
        }} > {output.complex}

        kmc_tools -t{threads} complex {output.complex}
        """

rule dump_combined_kmers:
    input:
        expand("results/kmc_combine_dbs/{{species}}.{suffix}", suffix=["kmc_pre", "kmc_suf"]),
    output:
        "results/dump_combined_kmers/{species}.txt"
    conda:
        "../envs/kmc.yaml" 
    shell:
        """
        # dump all k-mers to text file
        kmc_tools transform results/kmc_combine_dbs/{wildcards.species} dump {output}
        """

rule prejoin:
    input:
        comb=expand("results/dump_combined_kmers/{species}.txt", species = lookup(query="BioSample == {ID}", within=reads, cols="Species")),
        counts="results/kmc/{ID}.txt"
    output:
        "results/prejoin/{ID}.txt"
    shell:
        "join -t $'\t' -a1 -a2 -e '0' -o auto {input.comb} {input.sample} | cut -f 3 > {output}"

rule paste:
    input:
        kmer_list="results/dump_combined_kmers/{species}.txt",
        kmer_dumps=expand("results/prejoin/{ID}.txt", ID=lookup(query="Species == {species}", within=reads, cols="BioSample"))
    output:
       "results/paste/{species}.txt"
    shell:
        "paste <(cut -f 1 {input.kmer_list}) {input.kmer_dumps} > {output}"
