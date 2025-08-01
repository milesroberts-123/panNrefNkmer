# Notes

* Screening for more contaminants requires more memory because you have to load the whole database of contaminating k-mers into memory

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

- [x] samtools stats for read alignment, bam files

- [x] bcftools stats for ingroup and outgroup, variant and invariant sites

- [x] split vcfs into variant and invariant sites

- [x] variant filtering

- [x] add fastqc after removing contamination step

- [x] add localrules

- [x] add schema

- [x] before bwa mem, check if paired end files need to have names rewritten so that read names match

- [x] estimate DFE with fastDFE

- [x] add rules to use pggb

- [ ] estimate substitution rate in a manner that's robust to repeat sequences

- [ ] add target rules to only process dna-seq or only rna-seq

- [ ] add dN/dS calculation by comparing 1:1 orthologs with arabidopsis lyrata reference: download outgroup, get protein sequences, apply orthofinder, get single copies, use pal2nal to get codon alignments, calculate dN/dS

- [ ] add single-cell analysis with salmon-alevin

- [ ] rWCVP

- [ ] figure out batching, or maybe job grouping?

- [ ] increase priority for jobs that eliminate temporary files

- [ ] use fastk k-mer map to determine regions covered by k-mers

- [ ] benchmarks for various steps

- [ ] use smudgeplot to calculate lower and upper coverage bounds

- [ ] split fastp into two steps to avoid storing extra temporary files

- [ ] check if salmon quantmerge is fixed to allow custom column merge and gene level quantification

- [ ] remove k-mers from cds sequences?

- [ ] add graph aligner as alternative to vg giraffe?

- [ ] add vg map as alternative to vg giraffe?
