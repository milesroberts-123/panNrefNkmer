import fastdfe as fd
import click

# define click options
@click.command(context_settings={'show_default': True})
@click.option("-i", "--input", default=None, help="Path to k-mer count file", multiple=False)
@click.option("-n","--num-hash", default=None, help="Number of hash functions", type = click.INT)
@click.option("-s","--array-size", default=None, help="Number of elements in array", type = click.INT)
@click.option("-o", "--output", default=None, help="Path to output file")

def main(input, output, num_hash, array_size):

    # instantiate parser
    p = fd.Parser(
        n=8,
        vcf=url + "resources/genome/betula/biallelic.polarized.subset.50000.vcf.gz?raw=true",
        fasta=url + "resources/genome/betula/genome.subset.1000.fasta.gz?raw=true",
        gff=url + "resources/genome/betula/genome.gff.gz?raw=true",
        annotations=[
            fd.DegeneracyAnnotation(),
            fd.MaximumLikelihoodAncestralAnnotation(
                outgroups=["ERR2103730"],
                n_ingroups=10
            )
        ],
        stratifications=[fd.DegeneracyStratification()]
    )

    # parse SFS
    spectra: fd.Spectra = p.parse()

    # visualize SFS
    spectra.plot();

# create inference object
inf = fd.BaseInference(
    sfs_neut=spectra[['neutral.*']].merge_groups(1),
    sfs_sel=spectra[['selected.*']].merge_groups(1),
    do_bootstrap=True,
)

# run inference
inf.run();

inf.plot_discretized();

inf.plot_sfs_comparison();

if __name__ == '__main__':
    main()
