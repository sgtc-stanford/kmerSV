# main script for kmerSV

import argparse
from .kmerSV_plot import kmersv_plot_function
from .kmerSV_pangenome import kmersv_pangenome_function


def main():
    parser = argparse.ArgumentParser(description='''
        KmerSV: A tool for kmer-based visualization and annotation of structural variants (SVs) and pangenomes.
        
        This tool provides two modes:
        1. plot: Generates a kmer-based SV plot using a given reference genome.
        2. pangenome: Performs a kmer-based SV plot for a pangenome reference.

        Choose a mode and type '-h' after the mode for mode-specific help.
    ''')
    subparsers = parser.add_subparsers(dest="mode")

    # kmersv plot subparser 
    plot_parser = subparsers.add_parser('plot', help='Generate a kmer-based SV plot.')
    plot_parser.add_argument('-r', '--reference', help='reference sequence', required=True)
    plot_parser.add_argument('-i', '--input', help='target sequence', required=True)
    plot_parser.add_argument('-o', '--output', help='output SV plot', required=True)
    plot_parser.add_argument('-a', '--annotation', help='Genomic Annotation File', required=False, default='data/annotation/annotation.bed')
    plot_parser.add_argument('-c', '--chr', help='chromosome region')
    plot_parser.add_argument('-s', '--start', help='start position')
    plot_parser.add_argument('-e', '--end', help='end position')
    plot_parser.add_argument('-k', '--kmer', help='kmer legnth', required=True, default=31)
    plot_parser.add_argument('-f', '--text', help='output text file with detailed information', required=True)

    # kmersv pangenome subparser
    pangenome_parser = subparsers.add_parser('pangenome', help='Generate kmer-based SV plot for pangenome reference.')
    pangenome_parser.add_argument('-i', '--input', help='input fasta file with all sample paths', required=True)
    pangenome_parser.add_argument('-o', '--output', help='output SV plot using pangenome reference', required=True)
    pangenome_parser.add_argument('-a', '--annotation', help='Genomic Annotation File', default='data/annotation/annotation.bed')
    pangenome_parser.add_argument('-c', '--chr', help='chromosome region', required=True)
    pangenome_parser.add_argument('-s', '--start', help='start position', required=True)
    pangenome_parser.add_argument('-e', '--end', help='end position', required=True)
    pangenome_parser.add_argument('-k', '--kmer', help='kmer length', default=31)

    args = parser.parse_args()

    if args.mode == "plot":
        kmersv_plot_function(args.reference, args.input, args.output, args.annotation, args.chr, args.start, args.end, args.kmer, args.text)
    elif args.mode == "pangenome":
        kmersv_pangenome_function(args.input, args.output, args.annotation, args.chr, args.start, args.end, args.kmer)

if __name__ == "__main__":
    main()