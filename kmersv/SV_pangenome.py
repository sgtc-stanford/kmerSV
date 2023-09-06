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