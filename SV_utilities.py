import re
from collections import Counter

def reverse_complement(dna):
    complement = str.maketrans('ATCG', 'TAGC')
    return dna.translate(complement)[::-1]

# a function to retreive the reference genome froma fasta file
# input: fasta file name
# output: a dictionary with key as chromosome name and value as the sequence
def get_genome(fasta_file):
    ref_genome = {}
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                chr_name = line.strip().split()[0][1:]
                ref_genome[chr_name] = ''
            else:
                ref_genome[chr_name] += line.strip()
    return ref_genome   


# a function to extract the kmer from a sequence with interval length l 
# the output are two lists: one for the unique kmer and one for the non-unique kmer     
# input: a sequence, kmer length, interval length, chromosome name, start position, end position
# output: a list of unique kmers at every l away, a list of kmer positions, and its strand
#         a list of non-unique kmers at every l away, a list of kmer positions, and its strand
def get_kmer_interval(seq, k, l, chr, start, end, strand):
    kmer_list = []
    kmer_pos = []
    kmer_list_unique = []
    kmer_pos_unique = []
    kmer_strand = []
    kmer_list_nonunique = []
    for i in range(0, len(seq)-k+1, l):
        kmer = seq[i:i+k]
        kmer_list.append(kmer)
        kmer_pos.append(start+i)
        if kmer not in kmer_list_unique and kmer not in kmer_list_nonunique:
            kmer_list_unique.append(kmer)
            kmer_pos_unique.append(start+i)
        elif kmer in kmer_list_unique and kmer not in kmer_list_nonunique:
            kmer_list_nonunique.append(kmer)
            kmer_pos_unique.pop(kmer_list_unique.index(kmer))
            kmer_list_unique.pop(kmer_list_unique.index(kmer))
    return kmer_list, kmer_pos, kmer_list_unique, kmer_pos_unique,\
          kmer_list_nonunique, strand

# For each kmer in the kmer list, find the position of the kmer in the query genome. 
# Search for both the query genome and the reverse complement of the query genome.
# if there are multiple matches, record all the positions
# if the kmer is found in the query genome, record the position of the kmer in the query genome and its strand
# if the kmer is found in the reverse complement of the query genome, record the position of the kmer in the the query genome and its strand
# if the kmer is not found in the query genome or its reverse complement, record the position of the kmer as -1 and its strand as 'NA'
# input: a list of kmers, a query genome, chromosome name, start position, end position
# output: a list of positions of the kmers in the query genome, a list of strands of the kmers in the query genome
def get_kmer_pos(kmer_list, ref_kmer_pos, query_genome, chr_name, start):
    kmer_pos = []
    kmer_pos_nonrepeat = []
    ref_kmer_pos_nonrepeat = []
    query_count_nonrepeat = []
    kmer_strand_nonrepeat = []
    kmer_strand = []
    kmer_list_new = []
    query_count = []
    ref_kmer_pos_new = []
    query_genome_rev = reverse_complement(query_genome[chr_name])
    for i in range(len(kmer_list)):
        rev_flag = 1
        num_match = 0
        # if it finds more than one match, skip it for now
        if sum(1 for _ in re.finditer(kmer_list[i], query_genome[chr_name])) > 1:
            continue
        for match in re.finditer(kmer_list[i], query_genome[chr_name]):
            if (match.start()+start) not in kmer_pos:
                kmer_pos.append(match.start()+start)
                if num_match == 0: 
                    kmer_pos_nonrepeat.append(match.start()+start)
                    ref_kmer_pos_nonrepeat.append(ref_kmer_pos[i])
                    kmer_strand_nonrepeat.append('+')
                kmer_strand.append('+')
                kmer_list_new.append(kmer_list[i])
                ref_kmer_pos_new.append(ref_kmer_pos[i])
                rev_flag = 0
                num_match += 1
        if rev_flag == 1:
            if sum(1 for _ in re.finditer(kmer_list[i], query_genome_rev)) > 1:
                continue
            for match in re.finditer(kmer_list[i], query_genome_rev):
                if (match.start()+start) not in kmer_pos:
                    kmer_pos.append(match.start()+start)
                    if num_match == 0: 
                        kmer_pos_nonrepeat.append(match.start()+start)
                        ref_kmer_pos_nonrepeat.append(ref_kmer_pos[i])
                        kmer_strand_nonrepeat.append('-')
                    kmer_strand.append('-')
                    kmer_list_new.append(kmer_list[i])
                    ref_kmer_pos_new.append(ref_kmer_pos[i])
                    num_match += 1
        if num_match != 0: 
            query_count.extend([num_match]*num_match)
            query_count_nonrepeat.append(num_match)
    return kmer_list_new, ref_kmer_pos_new, kmer_pos, kmer_strand, query_count, \
        kmer_pos_nonrepeat, ref_kmer_pos_nonrepeat, query_count_nonrepeat, kmer_strand_nonrepeat


# find whether the positive or negative strand is dominant in the query genome
# if the positive strand is dominant, the query genome is in the same orientation as the reference genome
# if the negative strand is dominant, the query genome is in the reverse complement orientation as the reference genome
# input: a list of strands of the kmers in the query genome
# output: the dominant strand
def get_minor_strand(kmer_strand):
    strand_count = Counter(kmer_strand)
    if strand_count['+'] > strand_count['-']:
        return '+', '-'
    else:
        return '-', '+'
        

# a function to detect the SVs in the query genome using the kmer positions, the strand of the query genome, and the reference genome
# input: a list of kmer positions in the reference genome, its position in the query genome, and its strand in the query genome, 
# and number of occurances in the query genome
# output: a list of SVs in the query genome
def detect_structural_variants(ref_positions, contig_positions, strand_info, kmer_counts, l):
    variants = []
    ref_idx = 0
    contig_idx = 0

    inv_begin_flag = 0
    inv_finish_flag = 0
    inversion_length = 0
    inversion_start_ref = -1
    inversion_start_ctg = -1
    inversion_end_ref =  -1
    inversion_end_ctg =  -1

    dup_start_ctg = -1
    dup_end_ctg =  -1
    dup_start_ref = -1
    dup_length = 0
    dup_begin_flag = 0
    dup_finish_flag = 0
    dup_count = []

    while ref_idx < len(ref_positions) - 1 and contig_idx < len(contig_positions) - 1:
        ref_diff = ref_positions[ref_idx + 1] - ref_positions[ref_idx]
        contig_diff = contig_positions[contig_idx + 1] - contig_positions[contig_idx]

        # check if there are any tandem or duplicated regions
        if kmer_counts[ref_idx] == 1 and kmer_counts[ref_idx+1] > 1:
            dup_begin_flag = 1
            dup_start_ref = ref_positions[ref_idx]
        elif kmer_counts[ref_idx] >1 and kmer_counts[ref_idx+1] == 1:
            dup_finish_flag = 1
            dup_end_ref = ref_positions[ref_idx]

        # check if the inversion is finished or not
        if strand_info[contig_idx] == '+' and strand_info[contig_idx + 1] == '-':
            inv_begin_flag = 1
            inversion_start_ref = ref_positions[ref_idx + 1]
            inversion_start_ctg = contig_positions[contig_idx + 1]
        elif strand_info[contig_idx] == '-' and strand_info[contig_idx + 1] == '+':
            inv_finish_flag = 1
            inversion_end_ref = ref_positions[ref_idx]
            inversion_end_ctg = contig_positions[contig_idx]
        
        # update the dup_count
        if dup_begin_flag == 1 and kmer_counts[ref_idx]>1:
            dup_count.append(kmer_counts[ref_idx])
             
        if ref_diff == contig_diff:
            if inv_begin_flag == 1:
                inversion_length = inversion_length + ref_diff
            if dup_begin_flag == 1:
                dup_length = dup_length + contig_diff
            ref_idx += 1
            contig_idx += 1
        elif ref_diff < contig_diff and dup_begin_flag == 0:
            insertion_length = contig_diff - ref_diff
            variants.append(("insertion", ref_positions[ref_idx], ref_positions[ref_idx+1], insertion_length))
            ref_idx += 1
            contig_idx += 1
            
        elif ref_diff > contig_diff and dup_begin_flag == 0:
            deletion_length = ref_diff - contig_diff
            variants.append(("deletion", ref_positions[ref_idx], ref_positions[ref_idx+1], deletion_length))
            ref_idx += 1
            contig_idx += 1
        else:
            ref_idx += 1
            contig_idx += 1

        if inv_finish_flag:
                variants.append(("inversion", inversion_start_ref, inversion_end_ref, inversion_length))
                inversion_length = 0
                inv_begin_flag = 0
                inv_finish_flag = 0
        if dup_finish_flag:
                variants.append(("tandem_duplication", dup_start_ref, dup_end_ref, dup_end_ref - dup_start_ref, Counter(dup_count).most_common(1)[0][0] - 1))
                dup_length = 0
                dup_begin_flag = 0
                dup_finish_flag = 0
                dup_count = []

    return variants

def find_index(lst, value):
    try:
        index = lst.index(value)
    except ValueError:
        index = -1  
    return index