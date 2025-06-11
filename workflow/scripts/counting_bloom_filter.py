# counting bloom filter function from: https://github.com/williarj/kmers2024/blob/main/diversity_metrics/measure_diversity.py
import click
import pandas as pd
import sys
import numpy as np
import mmh3

#def load_kmers(filename):
#    print("Loading k-mer counts...")
#    df = pd.read_csv(filename,sep="\t",header=None,names=["kmer","counts"])
#
#    return df

def counting_bloom_filter(filename,num_hash,array_size):
    print("Calculating counting bloom filter...")

    # initialize final array
    final_array = np.zeros(array_size,dtype=np.uint64)
    
    lines_read = 0
    with open(filename) as f:
        for line in f:
            lines_read += 1

            if lines_read % 1000000 == 0:
                print(f"Processed {lines_read} lines...")
                   
            split_line = line.split()
            kmer = split_line[0]
            count = int(split_line[1])
            
            for k in range(0, num_hash):
                index = mmh3.hash(kmer,k,signed=False)%array_size
                final_array[index] += count

    # iterate over k-mers and calculate hashes for each
#    for i in range(0,total):
#        kmer = kmers[i]
#        count = counts[i]
#        for k in range(0, num_hash):
#            index = mmh3.hash(kmer,k,signed=False)%array_size
#            final_array[index] += count

    return final_array

#def bray_curtis(df1, df2):

#def jaccard(df1, df2):

#def cosine(df1, df2):

# define click options
@click.command(context_settings={'show_default': True})
@click.option("-i", "--input", default=None, help="Path to k-mer count file", multiple=False)
@click.option("-n","--num-hash", default=None, help="Number of hash functions", type = click.INT)
@click.option("-s","--array-size", default=None, help="Number of elements in array", type = click.INT)
@click.option("-o", "--output", default=None, help="Path to output file")

def main(input, output, num_hash, array_size):
    # load k-mers
    #df = load_kmers(input)

    # calculating counting bloom filter
    cbf = counting_bloom_filter(input,num_hash,array_size)

    # write counting bloom filter to disk
    print("Writing CBF to disk...")
    np.savetxt(output, cbf, delimiter=',', fmt='%d')
    print("Done! :D")

if __name__ == '__main__':
    main()
