# -*- coding: utf-8 -*-
"""
Genefinder assignment
Software Design Spring 2017
First REAL PYTHON CODE WAAT WAAAT
@author: Alisha Pegan Jan 26th, 2017


"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffle the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
        did not add addition docstring for doctest
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
        """

    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'T':
        return 'A'
    elif nucleotide == 'C':
        return 'G'
    elif nucleotide == 'G':
        return 'C'
    else:
        return (' ')


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence
        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
        did not add addition docstring for doctest
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'

    """
    # declare variable outside of loop or else it will be lost
    # seq will be blank string until further declared
    seq = ''

    for c in dna:
    # seq is old one plus new nucleotide
        seq += get_complement(c)
    # if return is one indent in, then return will pop out of the loop
    return seq[::-1]
    # dnaseq = list(range(len(dna)))


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
        did not add addition docstring for doctest
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    orfseq = ''
    # go through the dna sequence in threes. one three is a codeon
    # if codeon is a stop code, return the sequence
    # if codeon is not a stop code, add on to the end of the sequence
    for i in range(0, len(dna), 3):
        codeon1 = dna[i:i+3]
        if codeon1 == 'TAG' or codeon1 =='TAA' or codeon1 == 'TGA':
            return orfseq
        else:
            orfseq = orfseq + codeon1
    return orfseq


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.
        did not add addition docstring for doctest, because unsure what the
        correct answer should be

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']

    """

    i = 0
    frameseq = ''
    framelist = []

# go through the dna in threes, if one of the codeon is ATG
# send it to the rest_of_orf function, and assign the return value
# to dnaseq. then, increase the index to the end of the dnaseq
# then put the result in a list called framelist
# if the codeon is not ATG, then increase the index to the next three
    while i < len(dna):
        codeon1 = dna[i:i+3]
        if codeon1 == 'ATG':
            frameseq = dna[i:]
            dnaseq = rest_of_ORF(frameseq)
            i = i + len(dnaseq)
            framelist.append(dnaseq)
        else:
            i = i + 3
    return framelist


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
        did not add addition docstring for doctest

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    dnalist = []
    # use find_all_orfs function for all 3 potential frames
    frame1 = find_all_ORFs_oneframe(dna[0:])
    frame2 = find_all_ORFs_oneframe(dna[1:])
    frame3 = find_all_ORFs_oneframe(dna[2:])

    # add the results to the dnalist
    dnalist.extend(frame1)
    dnalist.extend(frame2)
    dnalist.extend(frame3)
    return dnalist


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
        did not add addition docstring for doctest
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    strandlist = []
    strand1 = find_all_ORFs(dna)
    reversestrand = get_reverse_complement(dna)
    strand2 = find_all_ORFs(reversestrand)

    strandlist.extend(strand1)
    strandlist.extend(strand2)
    return strandlist


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
        did not add addition docstring for doctest
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    longest = max(find_all_ORFs_both_strands(dna), key=len)
    return longest


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence
        no docstring to test for doctest because every dna strand is different
        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF
        with open will open any files
        for strings can use "for s in string"
        for numbers can use range
        do not use longest[i] because it will think you are trying to find
        something that does not exist, rather, use .append
        """
    with open(dna, 'r') as f:
        read_data = f.read()
        # print(read_data)
    local_longest = []
    for i in range(num_trials):
        # print(i)
        trial = shuffle_string(read_data)
        local_longest.append(longest_ORF(trial))
    # print(local_longest)
#    max(find_all_ORFs_both_strands(dna), key=len)
    maxlength = max([len(s) for s in local_longest])
    # print(maxlength)
    return maxlength


#longest_ORF_noncoding('data/X73525.fa', 10)


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment
        did not add addition docstring for doctest

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    aa_sequence = ''
    i = 0
    while i+1 < len(dna):
        codeon = dna[i:i+3]
        if len(codeon) == 3:
            # print(codeon)
            amino_acid = aa_table[codeon]
            aa_sequence += amino_acid
            i += 3
        else:
            break
    return aa_sequence


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    threshold = longest_ORF_noncoding(dna, 1500)
#    print(threshold)
    if len(dna) >= threshold:
        ORF = longest_ORF(dna)
        seq = coding_strand_to_AA(ORF)
        return seq
    else:
        print('Protein not found')



if __name__ == "__main__":
    import doctest
    # doctest.run_docstring_examples(gene_finder, globals(), verbose= True)
    doctest.testmod()
    gene_finder('data/X73525.fa')
