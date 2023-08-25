import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patches as patches
import re
import argparse
from collections import Counter
import pandas as pd
import SV_pangenome as SV

if __name__ == '__main__':
    # take arguments from command line
    parser = argparse.ArgumentParser(description='Plot kmer SVs')
    parser.add_argument('-i', '--input', help='input fasta file with all sample paths', required=True)
    parser.add_argument('-o', '--output', help='output SV plot using pangenome reference', required=True)
    parser.add_argument('-a', '--annotation', help='Genomic Annotation File', required=False, default='data/annotation/MANE_exon_intron.st.bed')
    parser.add_argument('-c', '--chr', help='chromosome region')
    parser.add_argument('-s', '--start', help='start position')
    parser.add_argument('-e', '--end', help='end position')
    parser.add_argument('-k', '--kmer', help='kmer legnth', required=True, default=31)
    args = parser.parse_args()

    k = int(args.kmer)
    path = args.input
    # read the first sequence in fasta file as reference genome
    ref_genome = SV.get_genome(path)

    # number of samples in pangenome
    num_samples = len(ref_genome.keys())
    # choose a colormap
    cmap = cm.get_cmap('rainbow', num_samples)

    # extract the chromosome name and the start and end position from the reference genome
    ref_name = list(ref_genome.keys())[0]
    ref_chr = ref_name.split(':')[0]
    ref_start = int(ref_name.split(':')[1].split('-')[0])
    ref_end = int(ref_name.split(':')[1].split('-')[1])
    ref_strand = '+'

    ref_kmer_list, ref_kmer_pos, ref_kmer_list_unique, \
    ref_kmer_pos_unique, ref_kmer_list_nonunique, ref_kmer_strand = \
        SV.get_kmer_interval(ref_genome[ref_name], k, 1, ref_chr, ref_start, ref_end, ref_strand)

    # plot the kmer positions in the reference genome and the query genome (with unique kmers in green and non-unique kmers in red)
    # without using the scientific notation on y-axis
    fig, ax = plt.subplots(figsize=(5, 5), dpi=1000)
    # get the max of the list 
    max_val = max(ref_kmer_pos_unique)
    min_val = min(ref_kmer_pos_unique)
    plt.plot(np.arange(min_val, max_val, 1), np.arange(0, max_val-min_val, 1), markersize=1.5, marker='o', color='black')
    # extract everything until second #
    label_ref = ref_name.split('#')[0] + ' chr19 \n'
    label_var = ['', '', '', '', '']
    color_dict = ['red', 'green', 'orange', 'blue', 'purple', 'black']
    variant_idx = {}
    count = [0, 0, 0, 0, 0, 0]
    start_var = [0, 0, 0, 0, 0, 0]

    ref_kmer_pos_unique_first = ref_kmer_pos_unique[0]
    overlap_max  = 1000000000000
    idx = 0
    # for the rest of sequences in the file, extract the query genome
    for i in range(0, num_samples):
        # extract the query genome of interest from a fasta file
        query_chr_name = list(ref_genome.keys())[i]
        query_chr = query_chr_name.split(':')[0]
        query_start = int(query_chr_name.split(':')[1].split('-')[0])
        query_end = int(query_chr_name.split(':')[1].split('-')[1])
        query_strand = '+'
        query_genome = ref_genome
        # use the uniqe kmers in the reference to find the positions of the kmers in the query genome
        kmer_list_unique_new, ref_kmer_pos_unique_new, query_kmer_pos_unique, kmer_strand_unique, query_count_unique, \
        kmer_pos_unique_nonrepeat, ref_kmer_pos_unique_nonrepeat, query_count_unique_nonrepeat, kmer_strand_unique_nonrepeat = \
        SV.get_kmer_pos(ref_kmer_list_unique, ref_kmer_pos_unique, query_genome, query_chr_name, query_start)

        ref_kmer_pos_unique_new = np.array(ref_kmer_pos_unique_new)
        query_kmer_pos_unique = np.array(query_kmer_pos_unique)
        offset = 0
        if query_kmer_pos_unique[0] > query_start:
            offset = query_kmer_pos_unique[0] - query_start

        if len(query_kmer_pos_unique) == 0:
            continue
        # offset the query genome position so that query_kmer_pos_unique[0] = ref_kmer_pos_unique_new[0]
        
        query_kmer_pos_unique = query_kmer_pos_unique - query_kmer_pos_unique[0] + ref_kmer_pos_unique_new[0]
        # find the index where query_kmer_pos_unique is different from ref_kmer_pos_unique_new
        idx_diff = np.where(ref_kmer_pos_unique_new != query_kmer_pos_unique)[0]
        # find the index where query_kmer_pos_unique is the same as ref_kmer_pos_unique_new
        idx_same = np.where(ref_kmer_pos_unique_new == query_kmer_pos_unique)[0]
        # offset both ref_kmer_pos_unique_new and query_kmer_pos_unique so the first element is 0

        query_kmer_pos_unique = query_kmer_pos_unique - np.min(query_kmer_pos_unique)
        # ref_kmer_pos_unique_new = ref_kmer_pos_unique_new - np.min(ref_kmer_pos_unique_new)
        offset = np.min(ref_kmer_pos_unique_new)

        # grid off
        plt.grid(False)
        plt.rcParams['axes.formatter.useoffset'] = False
        color = 'black'
        if len(idx_diff) != 0:
            if query_kmer_pos_unique.max() not in variant_idx.keys():
                variant_idx[query_kmer_pos_unique.max()] = idx
                idx += 1
            idx_var = variant_idx[query_kmer_pos_unique.max()]
            color = color_dict[idx_var]
            count[idx_var] += 1
            start_var[idx_var] = ref_kmer_pos_unique_new[idx_diff[0]]
            if color == 'green':
                # remove the first three points from idx_diff
                idx_diff = idx_diff[6:]
            plt.plot(ref_kmer_pos_unique_new[idx_diff], query_kmer_pos_unique[idx_diff], 'o', markersize=1.5, color=color)
            plt.plot([ref_kmer_pos_unique_new[idx_same[-1]], ref_kmer_pos_unique_new[idx_diff[0]]], \
                [query_kmer_pos_unique[idx_same[-1]], query_kmer_pos_unique[idx_diff[0]]], \
                    color=color, linestyle='dashed', linewidth=1.5)
            overlap_max = min(ref_kmer_pos_unique_new[idx_same[-1]], overlap_max)
            label_sex = 'paternal' if query_chr_name.split('#')[1] == '1' else 'maternal'
            label_var[idx_var] += query_chr_name.split('#')[0]
            if query_chr_name.split('#')[0] != 'CHM13':
                label_var[idx_var] += ' ' + label_sex
            elif query_chr_name.split('#')[0] == 'CHM13':
                label_var[idx_var] += '  chr19'
            label_var[idx_var] += '\n'
        else:
            count[-1] += 1
            if query_chr_name.split('#')[0] != 'CHM13':
                continue
                # label_ref += ' ' + label_sex
            elif query_chr_name.split('#')[0] == 'CHM13':
                label_sex = 'paternal' if query_chr_name.split('#')[1] == '1' else 'maternal'
                label_ref += query_chr_name.split('#')[0]
                label_ref += ' chr19'
            label_ref += '\n'
        plt.xlabel('Pangenome Reference', weight='bold', fontsize=13)
        plt.ylabel('Target Sequence', weight='bold', fontsize=13)
        # change the font size of the tick labels
        plt.tick_params(axis='both', which='major', labelsize=13)
        # change the font size of lables and title
        plt.rcParams.update({'font.size': 10})
        plt.locator_params(axis='x', nbins=5)
    # for each label, remove the last \n
    label_ref = label_ref[:-1]
    for i in range(len(label_var)):
        label_var[i] = label_var[i][:-1]
    plt.plot([], [], color='black', linewidth=1.5, label = label_ref)
    plt.plot([], [], color='red', linewidth=1.5, label = label_var[0])
    plt.plot([], [], color='green', linewidth=1.5, label = label_var[1])
    plt.plot([], [], color='orange', linewidth=1.5, label = label_var[2])
    plt.plot([], [], color='blue', linewidth=1.5, label = label_var[3])
    plt.plot([], [], color='purple', linewidth=1.5, label = label_var[4])
    # plt.legend(loc='best', prop={'size': 6})

    # add the count in counts to the plot 
    variant_idx[ref_kmer_pos_unique_new[-1] - offset] = 0
    for i in range(len(count)):
        plt.text(ref_kmer_pos_unique_new[-1] + 50, list(variant_idx.keys())[i], 'N='+str(count[i]),\
                ha='center', va='center', fontsize=9, color=color_dict[i], weight='bold')

    # Add gene annotation
    # Load BED file
    df = pd.read_csv("/mnt/ix1/Projects/M077_210115_kmer_assembly_analysis/P05_kmer_SVplot/MANE_exon_intron.st.bed", sep='\t')
    # User input for chromosome name and region
    chr_name = args.chr
    start_pos = int(args.start)	
    end_pos = int(args.end)
    plot_start = start_pos
    plot_end = end_pos
    # Filter dataframe based on user input
    filtered_df = df[(df['#chr'] == chr_name) & (df['start'] <= start_pos) & (df['end'] >= end_pos)]
    print(filtered_df)

    # if filtered_df is empty, find the closet two genes using the start and end position

    filtered_df_before = df[(df['#chr'] == chr_name) & (df['start'] <= start_pos)]
    gene_before = filtered_df_before.iloc[-1]['gene_name']
    gene_before_dist = start_pos - filtered_df_before.iloc[-1]['start']
    filtered_df_after = df[(df['#chr'] == chr_name) & (df['end'] >= end_pos)]
    gene_after = filtered_df_after.iloc[0]['gene_name']
    gene_after_dist = filtered_df_after.iloc[0]['end'] - end_pos
    # Prepare data for plotting
    plot_data = [(row['start'], row['end'], row['gene_name'], row['type'], row['exon_number']) for index, row in filtered_df.iterrows()]
    min_y = 0
    max_y = 1000


    if filtered_df.empty:
        # add horizontal arrows and texts to show the two closest genes , also add the distance in the texts
        plt.arrow(plot_start+30, min_y-40, -50, 0, head_width=10, head_length=6, fc='k', ec='k')
        plt.arrow(plot_end-30, min_y-40, 50, 0, head_width=10, head_length=6, fc='k', ec='k')
        # show the gene_before_dist and gene_after_dist in scientific notation
        gene_before_dist = "{:.2e}".format(gene_before_dist)
        gene_after_dist = "{:.2e}".format(gene_after_dist)

        plt.text(plot_start+40, min_y - 20, str(gene_before_dist)+'bp from '+gene_before, ha='center', va='center', fontsize=5, color='red', weight='bold')
        plt.text(plot_end-25, min_y - 20, str(gene_after_dist)+'bp to ' +gene_after, ha='center', va='center', fontsize=5, color='red', weight='bold')
    else:
        # Add regions to the plot
        for i, (start, end, gene_name, type, exon_number) in enumerate(plot_data):
            start = max(plot_start, start)
            end = min(plot_end, end)
            ax.add_patch(patches.Rectangle((start, min_y-60), end - start + 120, 12, edgecolor='blue', facecolor='blue', alpha=0.5))
            # add vertical bars on the x-axis to indicate the start and end of the region
            # plt.axvline(x=start, color='gray', linestyle='--', linewidth=0.75)
            # plt.axvline(x=end, color='gray', linestyle='--', linewidth=0.75)
            ax.text(start + 40 + (end - start) / 2, min_y - 25, f'{gene_name}:EX{exon_number}', ha='center', va='center', fontsize=9, color='blue', weight='bold')

    plt.ylim(min_y - 60, max_y+50)
    plt.xlim(start_pos-30, end_pos+180)
    plt.legend(loc='best', prop={'size': 4})

    plt.savefig(args.output, dpi=1000)



