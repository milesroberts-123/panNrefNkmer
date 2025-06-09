# counting bloom filter function from: https://github.com/williarj/kmers2024/blob/main/diversity_metrics/measure_diversity.py
import click
import pandas as pd
import sys
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
import itertools

def load_matrix(filename):
    print("Loading k-mer count matrix...")
    df = pd.read_csv(filename,sep="\t",header=None)
    return df

def bray_curtis(df1, df2):
    return 1 - 2*np.sum(df1, df2)/sum(np.minimum(df1, df2))

#def jaccard(df1, df2):

#def cosine(df1, df2):

# define click options
@click.command(context_settings={'show_default': True})
@click.option("-i", "--input", default=None, help="Path to k-mer count matrix", multiple=False)
@click.option("-o", "--output", default=None, help="Path to output file")

def main(input, output):
    # load matrix
    df = load_matrix(input)

    # loop over pairs of columns and calcualte distances
    n = len(df.columns)
    num_pairs = math.comb(n, 2)

    bc_total = 0
    cos_total = 0

    for col1, col2 in itertools.combinations(df.columns, 2):
        x = df[col1].values
        y = df[col2].values

        bc_total += bray_curtis(x,y)

        cos_total += 1 - cosine_similarity(x,y)

    # calculate average across all pairsimport math
    final_result = [bc_total/num_pairs, cos_total/num_pairs]

    # write to file
    with open(output, 'w') as file:
        file.write(','.join(final_result))

    print("Done! :D")

if __name__ == '__main__':
    main()
