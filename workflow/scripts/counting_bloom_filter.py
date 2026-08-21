"""Counting bloom filter over a k-mer count file."""

import click
import numpy as np
import mmh3


def counting_bloom_filter(input_path, num_hash, array_size):
    """Sum k-mer counts into a counting bloom filter array."""
    print("Calculating counting bloom filter...")

    # enforce that only large integers can be in array
    final_array = np.zeros(array_size, dtype=np.uint64)

    lines_read = 0

    with open(input_path) as f:
        for line in f:
            lines_read += 1

            if lines_read % 1000000 == 0:
                print(f"Processed {lines_read} k-mers...")

            kmer, count = line.split()
            count = int(count)

            # insert count at the index determined by each hash
            for k in range(num_hash):
                index = mmh3.hash(kmer, k, signed=False) % array_size
                final_array[index] += count

    return final_array


@click.command(context_settings={'show_default': True})
@click.option("-i", "--input", required=True, help="Path to k-mer count file")
@click.option("-n", "--num-hash", required=True, help="Number of hash functions", type=click.INT)
@click.option("-s", "--array-size", required=True, help="Number of elements in array", type=click.INT)
@click.option("-o", "--output", required=True, help="Path to output file")
def main(input, output, num_hash, array_size):
    """Write a counting bloom filter to disk."""
    cbf = counting_bloom_filter(input, num_hash, array_size)

    print("Writing CBF to disk...")
    np.savetxt(output, cbf, delimiter=',', fmt='%d')
    print("Done! :D")


if __name__ == '__main__':
    main()
