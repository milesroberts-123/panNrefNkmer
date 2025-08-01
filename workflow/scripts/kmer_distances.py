# counting bloom filter function from: https://github.com/williarj/kmers2024/blob/main/diversity_metrics/measure_diversity.py
print("Importing packages...")
import click
import pandas as pd
import sys
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
import itertools
import math

def load_matrix(filename):
    print("Loading k-mer count matrix...")
    df = pd.read_table(filename,header=None,sep=" ")
    return df

def bray_curtis(x, y):
    return 1 - (2*sum(np.minimum(x, y)))/sum(np.add(x, y))

#def jaccard(df1, df2):

#def cosine(df1, df2):

# define click options
@click.command(context_settings={'show_default': True})
@click.option("-i", "--input", default=None, help="Path to k-mer count matrix", multiple=False)
@click.option("-o", "--output", default=None, help="Path to output file")

def main(input, output):
    # load matrix
    df = load_matrix(input)
    print(df)
    # loop over pairs of columns and calcualte distances
    n = len(df.columns)
    num_pairs = math.comb(n, 2)

    bc_total = 0
    cos_total = 0
    
    print(f"Number of columns: {n}")
    print(f"Number of column pairs: {num_pairs}")
    
    print("Looping over pairs...")
    pairs_done = 0
    for col1, col2 in itertools.combinations(df.columns, 2):
        x = df[col1].values
        y = df[col2].values

        #print(x)
        #print(y)
        bc_total += bray_curtis(x,y)

        #cos_total += 1 - cosine_similarity([x],[y])
        
        pairs_done += 1
        
        if pairs_done % 1000 == 0:
            print(f"Processed {pairs_done} k-mers...")

    # calculate average across all pairs
    final_result = [str(bc_total/num_pairs)]

    # write to file
    print(f"Write result to {output}...")
    with open(output, 'w') as file:
        file.write(','.join(final_result))

    print("Done! :D")

if __name__ == '__main__':
    main()
