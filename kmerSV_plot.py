import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.cm as cm
import pandas as pd
import re
import argparse
import SV_utilities as SV

def kmersv_plot_function(reference, input, output, annotation, chr, start, end, kmer, text):
    """
    Generates a kmer-based structural variant (SV) plot using a given reference genome.

    Parameters:
    - reference (str): Path to the reference FASTA file.
    - input_file (str): Path to the input FASTA file containing all sample paths.
    - output_file (str): Path to save the output SV plot.
    - annotation (str, optional): Path to the Genomic Annotation File. Default is 'data/annotation/MANE_exon_intron.st.bed'.
    - chr_region (str, optional): Chromosome region for analysis.
    - start_pos (int, optional): Start position in the chromosome for analysis.
    - end_pos (int, optional): End position in the chromosome for analysis.
    - kmer_length (int, optional): Length of the kmer. Default is 31.
    - text (str, optional): Path to save the output text file with detailed information. Default is 'data/output/SV.txt'.

    Returns:
    None. The function saves the plot to the specified output_file.
    """

    k = int(kmer)
    ref_path = reference
    query_path = input

    # extract the reference genome from a fasta file
    ref_genome = SV.get_genome(ref_path)
    # extract the chromosome name and the start and end position from the reference genome
    ref_name = list(ref_genome.keys())[0]
    ref_chr = ref_name.split(':')[0]
    ref_start = int(ref_name.split(':')[1].split('-')[0])
    ref_end = int(ref_name.split(':')[1].split('-')[1])
    ref_strand = '+'

    # extract the query genome of interest from a fasta file
    query_genome = SV.get_genome(query_path)
    # extract the contig name and the start and end position from the query genome
    query_chr_name = list(query_genome.keys())[0]
    query_chr = query_chr_name.split(':')[0]
    query_start = int(query_chr_name.split(':')[1].split('-')[0])
    query_end = int(query_chr_name.split(':')[1].split('-')[1])
    query_strand = '+'

    ref_len = ref_end - ref_start
    l = max(1, int(ref_len/1000))



    ref_kmer_list, ref_kmer_pos, ref_kmer_list_unique, \
    ref_kmer_pos_unique, ref_kmer_list_nonunique, ref_kmer_strand = \
        SV.get_kmer_interval(ref_genome[ref_name], k, l, ref_chr, ref_start, ref_end, ref_strand)

    # use the uniqe kmers in the reference to find the positions of the kmers in the query genome
    kmer_list_unique_new, ref_kmer_pos_unique_new, query_kmer_pos_unique, kmer_strand_unique, query_count_unique, \
    kmer_pos_unique_nonrepeat, ref_kmer_pos_unique_nonrepeat, query_count_unique_nonrepeat, kmer_strand_unique_nonrepeat = \
    SV.get_kmer_pos(ref_kmer_list_unique, ref_kmer_pos_unique, query_genome, query_chr_name, query_start)



    # find all the ranges of index that have negative strand
    # and reverse the order of the query kmer position in these ranges
    # this is to make sure the query kmer position is in the same order as the reference kmer position
    # so that the kmer position can be plotted in the same order
    major_strand, minor_strand = SV.get_minor_strand(kmer_strand_unique_nonrepeat)
    if major_strand != '+':
        # change the '+' to '-' and '-' to '+' in the kmer_strand_unique list"
        kmer_strand_unique_nonrepeat = [kmer_strand.replace('+', 'x').replace('-', '+').replace('x', '-') for kmer_strand in kmer_strand_unique]
        kmer_strand_unique = [kmer_strand.replace('+', 'x').replace('-', '+').replace('x', '-') for kmer_strand in kmer_strand_unique]
    neg_strand_ranges = []
    for i in range(len(kmer_strand_unique)):
        if kmer_strand_unique[i] == minor_strand:
            if i == 0:
                neg_strand_ranges.append([i])
            elif kmer_strand_unique[i-1] == major_strand:
                neg_strand_ranges.append([i])
            else:
                neg_strand_ranges[-1].append(i)
    query_kmer_pos_unique_inv = query_kmer_pos_unique.copy()
    query_kmer_pos_inversion = []
    ref_kmer_pos_inversion = []
    for neg_strand_range in neg_strand_ranges:
        query_kmer_pos_unique_inv[neg_strand_range[0]:neg_strand_range[-1]+1] = query_kmer_pos_unique[neg_strand_range[0]:neg_strand_range[-1]+1][::-1]
        query_kmer_pos_inversion.extend(query_kmer_pos_unique[neg_strand_range[0]:neg_strand_range[-1]+1][::-1])
        ref_kmer_pos_inversion.extend(ref_kmer_pos_unique_new[neg_strand_range[0]:neg_strand_range[-1]+1])


    # detect the SVs in the query genome
    variants = SV.detect_structural_variants(ref_kmer_pos_unique_nonrepeat, kmer_pos_unique_nonrepeat, \
                                        kmer_strand_unique_nonrepeat, query_count_unique_nonrepeat, l)



    # remove SV that are too small
    v_max = max([v[3] for v in variants])
    # v_max = max(v_max, max([v[6] for v in variants]))
    threshold = 45
    variants = [v for v in variants if v[3] > threshold]


    # for variants that are insertions
    # check if the insertion is a duplication
    # if it is a duplication, then change the type to tandem_duplication
    # and change the length to the length of the duplication
    # specifcally check if the inserted sequence has lots of non-unique kmers in the nearby region
    # if it does, then it is a duplication
    # if it does not, then it is not a duplication
    inserted_kmer_query_pos = []
    inserted_kmer_ref_pos = []
    for i in range(len(variants)):
        if variants[i][0] == "insertion":
            # find the start position of insertion on the query genome
            # find the nearest
            insertion_start = query_kmer_pos_unique[np.abs(np.array(ref_kmer_pos_unique_new)-variants[i][1]).argmin()] - query_start
            insertion_end = query_kmer_pos_unique[np.abs(np.array(ref_kmer_pos_unique_new)-variants[i][2]).argmin()] - query_start
            insertion_start_nearby = max(0, insertion_start - 500)
            insertion_end_nearby = min(query_end, insertion_end + 500)
            # get the inserted sequence
            inserted_seq = query_genome[query_chr_name][insertion_start:insertion_end]
            # get the kmers in the inserted sequence
            inserted_kmer_list, inserted_kmer_pos, inserted_kmer_list_unique, \
            inserted_kmer_pos_unique, inserted_kmer_list_nonunique, inserted_kmer_strand = \
                SV.get_kmer_interval(inserted_seq, k, 1, query_chr, insertion_start+query_start, query_end, query_strand)
            # get the kmers in the nearby region
            nearby_kmer_list, nearby_kmer_pos, nearby_kmer_list_unique, \
            nearby_kmer_pos_unique, nearby_kmer_list_nonunique, nearby_kmer_strand = \
                SV.get_kmer_interval(query_genome[query_chr_name][insertion_start_nearby:insertion_end_nearby], k, 1, query_chr, insertion_start_nearby+query_start, query_end, query_strand)
            
            # get the number of non-unique kmers in the inserted sequence
            # save the coordinates of the non-unique kmers in the inserted sequence
            inserted_kmer_list_dup = [kmer for kmer in inserted_kmer_list if nearby_kmer_list.count(kmer) > 1]
            # if the number of non-unique kmers is more than 50% of the total number of kmers in the inserted sequence
            # then it is a duplication
            dup_flag = 1
            if len(inserted_kmer_list_dup) * 1.0 > len(inserted_kmer_list) * 0.5:
                variants[i] = ("duplication", variants[i][1], variants[i][2], variants[i][3])
            # if the number of non-unique kmers is more than 20% of the total number of kmers in the inserted sequence
            # then it is a complex SV
            elif len(inserted_kmer_list_dup) * 1.0 > len(inserted_kmer_list) * 0.1:
                variants[i] = ("complex:DUP_INS", variants[i][1], variants[i][2], variants[i][3])
            else:
                dup_flag = 0
            if dup_flag == 1:
                # find the position of the unique kmer in the inserted sequence
                # find the position of them in the reference genome
                for j in range(len(inserted_kmer_list)):
                    if inserted_kmer_list[j] in inserted_kmer_list_dup:
                        matches = re.finditer(inserted_kmer_list[j], ref_genome[ref_name])
                        num_match = 0
                        for match in matches:
                            if num_match < 1 and match.start() + ref_start >= variants[i][1] and match.start() + ref_start <= variants[i][2]:
                                inserted_kmer_query_pos.append(inserted_kmer_pos[j])
                                inserted_kmer_ref_pos.append(match.start() + ref_start)
                                num_match += 1
    print(variants)
    # write variants to the output file (tab-delimited)
    with open(text, 'w') as f:
        # write a header
        f.write('\t'.join(['SV type', 'ref_start', 'ref_end', 'length']) + '\n')
        for variant in variants:
            f.write('\t'.join([str(v) for v in variant]) + '\n')

    query_genome_rev = SV.reverse_complement(query_genome[query_chr_name])
    ref_pos_nonunique = []
    query_pos_nonunique = []
    ref_pos_nonunique_dup = []
    query_pos_nonunique_dup = []
    for kmer_dup in ref_kmer_list_nonunique:
        # find the indexes of all kmer_dup in the ref_kmer_list
        # then use the index to find the corresponding value in the ref_kmer_pos
        # ref_kmer_list is a list
        ref_pos_dup = []
        for i in range(len(ref_kmer_list)):
            if ref_kmer_list[i] == kmer_dup:
                ref_pos_dup.append(ref_kmer_pos[i])

        rev_flag = 1
        num_match = 0
        query_kmer_pos_dup = []
        for match in re.finditer(kmer_dup, query_genome[query_chr_name]):
            query_kmer_pos_dup.append(match.start()+query_start)
            rev_flag = 0
            num_match += 1
        if rev_flag == 1:
            for match in re.finditer(kmer_dup, query_genome_rev):
                query_kmer_pos_dup.append(match.start()+query_start)
                num_match += 1
        if num_match != 0: 
            query_kmer_pos_dup_used = []
            for i in range(len(ref_pos_dup)):
                # find the index of nearest value in the ref_kmer_pos_unique to the ref_pos_first
                # then use the index to find the corresponding value in the query_kmer_pos_unique
                query_pos_nearest = query_kmer_pos_unique[np.abs(np.array(ref_kmer_pos_unique_new)-ref_pos_dup[i]).argmin()]
                # find the nearest value in the query_kmer_pos_dup to the query_pos_nearest
                query_pos_dup_nearest = query_kmer_pos_dup[np.abs(np.array(query_kmer_pos_dup)-query_pos_nearest).argmin()]
                # add this point to the non-unique kmer plot
                ref_pos_nonunique.append(ref_pos_dup[i])
                query_pos_nonunique.append(query_pos_dup_nearest)
                # remove the value from the query_kmer_pos_dup
                query_kmer_pos_dup_used.append(query_pos_dup_nearest)
                # remove the redundant point with the same query_pose_nonqunique 
                # from the ref_pos_unique and query_pos_unique
                idx = SV.find_index(query_kmer_pos_unique, query_pos_dup_nearest)
                if idx != -1:
                    ref_kmer_pos_unique_new.pop(idx)
                    query_kmer_pos_unique.pop(idx)
                    kmer_list_unique_new.pop(idx)
            for k in list(set(query_kmer_pos_dup_used)):
                query_kmer_pos_dup.remove(k)
            # then add the rest of the points to the non-unique kmer plot
            for i in range(len(query_kmer_pos_dup)):
                # find the index of nearest value in the query_pos_nonunique to the query_kmer_pos_dup[i]
                # then use the index to find the corresponding value in the ref_pos_nonunique
                ref_pos_nearest = ref_pos_nonunique[np.abs(np.array(query_pos_nonunique)-query_kmer_pos_dup[i]).argmin()]
                # find the nearest value in the ref_pos_dup to the ref_pos_nearest
                ref_pos_dup_nearest = ref_pos_dup[np.abs(np.array(ref_pos_dup)-ref_pos_nearest).argmin()]
                # add it as duplicate to the nonunique kmer plot
                ref_pos_nonunique_dup.append(ref_pos_dup_nearest)
                query_pos_nonunique_dup.append(query_kmer_pos_dup[i])

    # Load BED file
    path_annotate = annotation
    df = pd.read_csv(path_annotate, sep='\t')

    # User input for chromosome name and region
    chr_name = chr
    start_pos = int(start)
    end_pos = int(end)


    # Filter dataframe based on user input
    filtered_df = df[(df['#chr'] == chr_name) & (df['start'] >= start_pos) & (df['end'] <= end_pos)]

    # if filtered_df is empty, find the closet two genes using the start and end position

    filtered_df_before = df[(df['#chr'] == chr_name) & (df['start'] <= start_pos)]
    gene_before = filtered_df_before.iloc[-1]['gene_name']
    gene_before_dist = start_pos - filtered_df_before.iloc[-1]['start']
    filtered_df_after = df[(df['#chr'] == chr_name) & (df['end'] >= end_pos)]
    gene_after = filtered_df_after.iloc[0]['gene_name']
    gene_after_dist = filtered_df_after.iloc[0]['end'] - end_pos

    # Prepare data for plotting
    plot_data = [(row['start'], row['end'], row['gene_name'], row['type'], row['exon_number']) for index, row in filtered_df.iterrows()]


    # plot the kmer positions in the reference genome and the query genome (with unique kmers in green and non-unique kmers in red)
    # without using the scientific notation on y-axis
    fig, ax = plt.subplots(figsize=(5, 5), dpi = 1000)
    # grid off
    plt.grid(False)
    plt.rcParams['axes.formatter.useoffset'] = False
    plt.plot(ref_kmer_pos_unique_new, query_kmer_pos_unique_inv, 'o', markersize=2.5)
    # plt.plot(query_pos_nonunique, ref_pos_nonunique, 'o', markersize=1, color='red')
    plt.plot(inserted_kmer_ref_pos, inserted_kmer_query_pos, 'o', markersize=2.5, color='orange')
    # plt.plot(ref_kmer_pos_inversion, query_kmer_pos_inversion, 'o', markersize=2.5, color='orange')
    plt.xlabel('Reference Sequence', weight='bold', fontsize=12)
    plt.ylabel('Target Sequence', weight='bold', fontsize = 12)
    # change the font size of the tick labels
    plt.tick_params(axis='both', which='major', labelsize=12)
    # change the font size of lables and title
    plt.rcParams.update({'font.size': 12})
    # plt.savefig(args.plot_path + "/" + query_chr_name + "_" +\
    #             "ref:" + ref_name + ".png", dpi=1200)
    # make the ticks on x-axis sparse
    plt.locator_params(axis='x', nbins=5)

    # find the min in query_kmer_pos_unique_inv
    min_y = min(query_kmer_pos_unique_inv)*0.9996
    max_y = max(query_kmer_pos_unique_inv)*1.0001

    if filtered_df.empty:
        # add horizontal arrows and texts to show the two closest genes , also add the distance in the texts
        plt.arrow(start_pos, min_y+400, -1000, 0, head_width=100, head_length=100, fc='k', ec='k')
        plt.arrow(end_pos, min_y+400, 1000, 0, head_width=100, head_length=100, fc='k', ec='k')
        # show the gene_before_dist and gene_after_dist in scientific notation
        gene_before_dist = "{:.2e}".format(gene_before_dist)
        gene_after_dist = "{:.2e}".format(gene_after_dist)

        plt.text(start_pos+500, min_y + 700, str(gene_before_dist)+'bp from '+gene_before, ha='center', va='center', fontsize=5, color='red', weight='bold')
        plt.text(end_pos-500, min_y + 700, str(gene_after_dist)+'bp to ' +gene_after, ha='center', va='center', fontsize=5, color='red', weight='bold')

    else:
        # Add regions to the plot
        offset = -9
        for i, (start, end, gene_name, type, exon_number) in enumerate(plot_data):
            ax.add_patch(patches.Rectangle((start, min_y+100), end - start, 40, edgecolor='blue', facecolor='blue', alpha=0.5))
            # add vertical bars on the x-axis to indicate the start and end of the region
            plt.axvline(x=start, color='gray', linestyle='--', linewidth=0.75)
            plt.axvline(x=end, color='gray', linestyle='--', linewidth=0.75)
            if type == 'exon':
                ax.text(start + (end - start) / 2 + offset, min_y + 720, f'{gene_name}:EX{exon_number}', ha='center', va='center', fontsize=7, color='blue', weight='bold')
                offset = 9
            else:
                ax.text(start + (end - start) / 2 + 4, min_y + 20, f'{gene_name}:INT', ha='center', va='center', fontsize=7, color='blue', weight='bold')
        # add horizontal arrows and texts to show the two closest genes , also add the distance in the texts
        plt.arrow(start_pos, min_y+960, -1000, 0, head_width=50, head_length=20, fc='k', ec='k', linewidth=2)
        plt.arrow(end_pos, min_y+960, 1000, 0, head_width=50, head_length=30, fc='k', ec='k', linewidth=2)
        # show the gene_before_dist and gene_after_dist in scientific notation
        gene_before_dist = "{:.2e}".format(gene_before_dist)
        gene_after_dist = "{:.2e}".format(gene_after_dist)
        plt.text(start_pos+1680, min_y + 2000, str(gene_before_dist)+'bp from '+gene_before, ha='center', va='center', fontsize=7, color='red', weight='bold')
        plt.text(end_pos-1530, min_y + 2000, str(gene_after_dist)+'bp to ' +gene_after, ha='center', va='center', fontsize=7, color='red', weight='bold')
    plt.ylim(min_y, max_y)
    plt.tight_layout()
    plt.savefig(output, dpi=1200)
