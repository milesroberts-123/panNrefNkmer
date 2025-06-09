rule kmer_distances:
    input:
        "cbf_table.txt"
    output:
        "kmer_distances.txt"
    conda:
        "../envs/cbf.yaml"
    shell:
        "scripts/kmer_distances.py --input {input} --output {output}"
