rule pggb:
    input:
        gz="results/separate_chrom/{chr}.fa.gz",
        tbi="results/separate_chrom/{chr}.fa.gz.tbi",
        triangle="results/mash/{chr}_triangle.txt"
    output:
        "results/pggb/{chr}/{chr}.fa.final.gfa"
    params:
        outdir = "results/pggb/{chr}"
    shell:
        """
        echo Getting max divergence
        MAXDIV=$(sed 1,1d  {input.triangle} | tr '\t' '\n' | grep "#" -v  | grep "Max" -v | LC_ALL=C sort -g -k 1nr | uniq | head -n 1)
        RESULT=$(echo "$MAXDIV - 2.0" | bc -l)

        echo $MAXDIV
        echo $RESULT
        
        pggb -i {input.gz} \
             -o {params.outdir} \
             -t {threads} \
             -p $RESULT \
             -s 5k
        """
