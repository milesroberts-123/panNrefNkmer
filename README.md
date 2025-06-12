# vg-kmer-ref-snakes

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)


A snakemake workflow to analyze short reads with a either a pangenome reference, linear reference, or no reference

# To-do

- [x] add mitochondrial and chloroplast genomes as contaminants

- [x] figure out vg surject

- [x] add kmer counts

- [x] add genotype calling for linear references

- [x] sort bam files

- [x] add MK test in degenotate?

- [x] create linear references for pangenome snp calling depending on path: Use `vg paths -F` to get fasta sequences for paths, then use that as the linear reference for genotype calling

- [x] counting bloom filter

- [x] k-mer distances

- [x] gene expression for different references

- [x] add chromosomes.tsv: two column file with assembly prefix in one column and chromosome names in the other column, useful to split workflow by chromosome 

- [x] split genotype calling by chromosome?

- [x] add bcftools concat to concat genotype calls split by chromosome

- [x] find read duplicates with samtools

- [x] multiqc

- [x] tajimas d

- [x] pi

- [x] watterson's theta

- [x] add profiles

- [x] for counting bloom filter, just read input file one line at a time

- [x] include outgroup in bcftools merge

- [ ] bcftools stats for ingroup and outgroup, variant and invariant sites

- [ ] samtools stats for read alignment, bam files

- [ ] benchmarks for various steps

- [ ] split fastp into two steps to avoid storing extra temporary files

- [ ] check if salmon quantmerge is fixed to allow custom column merge and gene level quantification

- [ ] add localrules

- [ ] add schema

- [ ] split vcfs into variant and invariant sites

- [ ] variant filtering

- [ ] dn/ds between references?

- [ ] remove k-mers from cds sequences?

- [ ] add graph aligner as alternative to vg giraffe?

- [ ] add vg map as alternative to vg giraffe?
