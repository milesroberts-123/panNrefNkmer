# vg-kmer-ref-snakes

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

- [ ] gene expression for different references

- [ ] add chromosomes.tsv: two column file with assembly prefix in one column and chromosome names in the other column, useful to split workflow by chromosome 

- [ ] dn/ds between references?

- [ ] tajima's D by gene

- [ ] synonymous diversity by gene

- [ ] nonsynonymous diversity by gene

- [ ] split genotype calling by chromosome?

- [ ] add bcftools concat to concat genotype calls split by chromosome

- [ ] find read duplicates with samtools

- [ ] remove k-mers from cds sequences?

- [ ] add graph aligner as alternative to vg giraffe?

- [ ] add vg map as alternative to vg giraffe?
