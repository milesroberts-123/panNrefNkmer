"""Mean pairwise k-mer distance metrics (Bray-Curtis, cosine) over a count matrix."""

print("Importing modules...")
import click
import pandas as pd
import numpy as np
import itertools
import math


def load_matrix(filename):
    """Load a whitespace-delimited k-mer count matrix (rows = k-mers, cols = samples)."""
    return pd.read_table(filename, header=None, sep=" ")


def bray_curtis(x, y):
    """Bray-Curtis dissimilarity between two count vectors."""
    x_norm = x / x.sum()
    y_norm = y / y.sum()
    return 1 - (2 * sum(np.minimum(x_norm, y_norm))) / sum(np.add(x_norm, y_norm))


def cosine_distance(x, y):
    """Cosine distance between two count vectors."""
    dot = np.dot(x, y)
    norm_x = np.linalg.norm(x)
    norm_y = np.linalg.norm(y)
    if norm_x == 0 or norm_y == 0:
        return 1.0
    return 1 - dot / (norm_x * norm_y)


@click.command(context_settings={'show_default': True})
@click.option("-i", "--input", required=True, help="Path to k-mer count matrix")
@click.option("-o", "--output", required=True, help="Path to output file")
def main(input, output):
    """Average Bray-Curtis and cosine distances over all column pairs."""
    print(f"Loading k-mer count matrix: {input}...")
    df = load_matrix(input)

    print("Counting number of columns...")
    n = len(df.columns)

    print("Calculating number of pairwise comparisions...")
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

        bc_total += bray_curtis(x, y)
        cos_total += cosine_distance(x, y)

        pairs_done += 1

        if pairs_done % 1000 == 0:
            print(f"Processed {pairs_done} pairs...")

    final_result = [str(bc_total / num_pairs), str(cos_total / num_pairs)]

    print(f"Write result to {output}...")
    with open(output, 'w') as file:
        file.write(','.join(final_result))

    print("Done! :D")


if __name__ == '__main__':
    main()
