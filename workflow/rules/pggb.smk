rule pggb:
    input:
        gz="{chr}.fa.gz",
        tbi="{chr}.fa.gz.tbi",
        triangle="{chr}_triangle.txt"
    output:
        "{chr}_pggb_results/{chr}.fa.final.gfa"
    params:
        outdir = "{chr}_pggb_results"
    shell:
        """
        echo Getting max divergence
        MAXDIV=$(sed 1,1d  {input.triangle} | tr '\t' '\n' | grep "#" -v  | grep "Max" -v | LC_ALL=C sort -g -k 1nr | uniq | head -n 1)
        RESULT=$(echo "$MAXDIV - 2.0" | bc -l)

        echo $MAXDIV
        echo $RESULT
        
        pggb -i {input.fa} \       # input file in FASTA format
             -o {params.outdir} \      # output directory
             -t {threads} \          # number of threads
             -p $RESULT \          # minimum average nucleotide identity for segments
             -s 5k \          # segment length for scaffolding the graph
        """
