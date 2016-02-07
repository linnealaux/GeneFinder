# -*- coding: utf-8 -*-
"""
Mini Project 1

@author: Linnea Laux

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide
        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    >>> get_complement('T')
    'A'
    >>> get_complement('G')
    'C'
    """
    if nucleotide == 'A':
        return 'T'
    if nucleotide == 'C':
        return 'G'
    if nucleotide == 'T':
        return 'A'
    if nucleotide == 'G':
        return 'C'


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence
        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    reverse = ''
    for nucleotide in dna[::-1]:
        a = get_complement(nucleotide)
        reverse = reverse + a
    return reverse

        

    


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.
        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    >>> rest_of_ORF("ATGAGA")
    'ATGAGA'
    """
    index = 0
    stop1 = 'TAG'
    stop2 = 'TAA'
    stop3 = 'TGA'
    while index < len(dna):
        codon = dna[index:index+3]
        if codon ==  stop1 or codon==stop2 or codon==stop3:
            return dna[0:index]
        else:
            index = index +3
    return dna

                

def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    index = 0
    ORF_list = []
    while index < len(dna):
        codon = dna[index:index+3]
        if codon == 'ATG':
            a = rest_of_ORF(dna[index:])
            ORF_list.append(a)
            index = index + len(a)
        else:
            index = index + 3
    return ORF_list



def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    dna1 = dna[1:]
    dna2 = dna[2:]
    one = find_all_ORFs_oneframe(dna)
    two = find_all_ORFs_oneframe(dna1)
    three = find_all_ORFs_oneframe(dna2)
    all_ORFs = one + two + three
    return all_ORFs


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    
    a = find_all_ORFs(dna)
    backwards = get_reverse_complement(dna)
    b = find_all_ORFs(backwards)
    all_orf_list = a+b
    return all_orf_list


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    ORFs = find_all_ORFs_both_strands(dna)
    lengths = []
    for item in ORFs:
        lengths.append(len(item))
    longest = max(lengths)
    for item in ORFs:
        if len(item) == longest:
            return item
        else:
            continue
    return False







def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    more_orfs = []
    orf_lengths = []
    for i in range(1,num_trials):
        nda = shuffle_string(dna)
        more_orfs.append(longest_ORF(nda))
    for q in range(1,len(more_orfs)):
        orf_lengths.append(len(more_orfs[q]))
    return max(orf_lengths)


    


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    index = 0
    AA_string = ''
    while index < len(dna):
        codon = dna[index:index+3]
        if len(codon) == 3:
            aa = aa_table[codon]
            AA_string = AA_string + aa
            index = index + 3
        else:
            break
    return AA_string


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    threshold = longest_ORF_noncoding(dna, 150)
    all_orf_list = find_all_ORFs_both_strands(dna)
    above_threshold = ''
    final_aa_list = []
    for j in all_orf_list:
        if len(j)>threshold:
            above_threshold = above_threshold + j
            final_aa_list.append(coding_strand_to_AA(above_threshold))
    return final_aa_list

# if __name__ == "__main__":
#     import doctest
#     doctest.testmod()
from load import load_seq
dna = load_seq("./data/X73525.fa")
print gene_finder(dna)
