# counting bloom filter function from: https://github.com/williarj/kmers2024/blob/main/diversity_metrics/measure_diversity.py
print("Importing packages...")
import click
import pandas as pd
import sys
import numpy as np
import mmh3

# Counting bloom filter function
def counting_bloom_filter(input,num_hash,array_size):
    print("Calculating counting bloom filter...")

    # initialize final array
    # enforce that only large integers can be in array
    final_array = np.zeros(array_size,dtype=np.uint64)
    
    # track number of lines read
    lines_read = 0
    
    # loop over every line of k-mer count file
    with open(input) as f:
        for line in f:

            # output progress every few lines        
            lines_read += 1

            if lines_read % 1000000 == 0:
                print(f"Processed {lines_read} k-mers...")
            
            # get k-mer and count           
            split_line = line.split()
            kmer = str(split_line[0])
            count = int(split_line[1])
            
            # loop over hash functions
            # insert count at the index determined by hash
            for k in range(0, num_hash):
                index = mmh3.hash(kmer,k,signed=False)%array_size
                final_array[index] += count

    return final_array

# define click options
@click.command(context_settings={'show_default': True})
@click.option("-i", "--input", default=None, help="Path to k-mer count file", multiple=False)
@click.option("-n","--num-hash", default=None, help="Number of hash functions", type = click.INT)
@click.option("-s","--array-size", default=None, help="Number of elements in array", type = click.INT)
@click.option("-o", "--output", default=None, help="Path to output file")

def main(input, output, num_hash, array_size):
    # calculating counting bloom filter
    cbf = counting_bloom_filter(input,num_hash,array_size)

    # write counting bloom filter to disk
    print("Writing CBF to disk...")
    np.savetxt(output, cbf, delimiter=',', fmt='%d')
    print("Done! :D")

if __name__ == '__main__':
    main()
