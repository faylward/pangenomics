# -*- coding: utf-8 -*-
__author__ = 'Jessica_Bryant'
__email__ = 'jessawbryant@gmail.com'


"""
This module holds dictionaries that contain codon and amino acid tables, and functions that
 caculate GC, N-ARSC, C-ARSC and Nc.

Original citations for calculated metrics:

   Baudouin-Cornu P, Surdin-Kerjan Y, Marliere P, Thomas D. 2001. Molecular evolution of protein atomic
    composition. Science 293 297–300.

   Wright F. 1990. The 'effective number of codons' used in a gene. Gene 87 23-29.

Updated: 5/29/2017 
"""

import itertools as itertools
from operator import add as add
import random as random
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd
import copy as cp


codon_dictionary = {

    # lists the synonymous family type (SF*), encoded amino acid (aa), ranking of nitrogen use of each codon relative
    # to other the synonymous codons encoding the same aa (GC_rank), number of nitrogen atoms on each encoded aa side
    # chain(N),  number of sulfur atoms on each encoded aa side chain (S), molecular weight of the encoded aa (MW).
    # refs:
    # Molecular Weight, N and S counts come from page 30 of 'Understanding Bioinformatics' by Zvelbil and Baum

    #D
    'GAT': {'SF': 'SF2', 'aa': 'D', 'GC_rank': 1, 'N': 0, 'S': 0, 'C': 2, 'MW': 133.1032},
    'GAC': {'SF': 'SF2', 'aa': 'D', 'GC_rank': 0, 'N': 0, 'S': 0, 'C': 2, 'MW': 133.1032},

    #E
    'GAA': {'SF': 'SF2', 'aa': 'E', 'GC_rank': 1, 'N': 0, 'S': 0,'C': 3, 'MW': 147.1299},
    'GAG': {'SF': 'SF2', 'aa': 'E', 'GC_rank': 0, 'N': 0, 'S': 0,'C': 3, 'MW': 147.1299},

    #S
    'TCT': {'SF': 'SF6', 'aa': 'S', 'GC_rank': 1, 'N': 0, 'S': 0, 'C': 1, 'MW': 105.0930},
    'TCC': {'SF': 'SF6', 'aa': 'S', 'GC_rank': 0, 'N': 0, 'S': 0, 'C': 1, 'MW': 105.0930},
    'TCA': {'SF': 'SF6', 'aa': 'S', 'GC_rank': 1, 'N': 0, 'S': 0, 'C': 1, 'MW': 105.0930},
    'TCG': {'SF': 'SF6', 'aa': 'S', 'GC_rank': 0, 'N': 0, 'S': 0, 'C': 1, 'MW': 105.0930},
    'AGT': {'SF': 'SF6', 'aa': 'S', 'GC_rank': 1, 'N': 0, 'S': 0, 'C': 1, 'MW': 105.0930},
    'AGC': {'SF': 'SF6', 'aa': 'S', 'GC_rank': 0, 'N': 0, 'S': 0, 'C': 1, 'MW': 105.0930},

    #T
    'ACT': {'SF': 'SF4', 'aa': 'T', 'GC_rank': 1, 'N': 0, 'S': 0, 'C': 2, 'MW': 119.1197},
    'ACC': {'SF': 'SF4', 'aa': 'T', 'GC_rank': 0, 'N': 0, 'S': 0, 'C': 2, 'MW': 119.1197},
    'ACA': {'SF': 'SF4', 'aa': 'T', 'GC_rank': 1, 'N': 0, 'S': 0, 'C': 2, 'MW': 119.1197},
    'ACG': {'SF': 'SF4', 'aa': 'T', 'GC_rank': 0, 'N': 0, 'S': 0, 'C': 2, 'MW': 119.1197},

    #Y
    'TAT': {'SF': 'SF2', 'aa': 'Y', 'GC_rank': 1, 'N': 0, 'S': 0, 'C': 7, 'MW': 181.1894},
    'TAC': {'SF': 'SF2', 'aa': 'Y', 'GC_rank': 0, 'N': 0, 'S': 0, 'C': 7, 'MW': 181.1894},

    #A
    'GCT': {'SF': 'SF4', 'aa':'A', 'GC_rank': 1, 'N': 0, 'S': 0, 'C': 1, 'MW': 89.0935},
    'GCC': {'SF': 'SF4', 'aa':'A', 'GC_rank': 0, 'N': 0, 'S': 0, 'C': 1, 'MW': 89.0935},
    'GCA': {'SF': 'SF4', 'aa':'A', 'GC_rank': 1, 'N': 0, 'S': 0, 'C': 1, 'MW': 89.0935},
    'GCG': {'SF': 'SF4', 'aa':'A', 'GC_rank': 0, 'N': 0, 'S': 0, 'C': 1, 'MW': 89.0935},

    #V
    'GTT': {'SF': 'SF4', 'aa': 'V', 'GC_rank': 0, 'N': 0, 'S': 0, 'C': 3, 'MW': 117.1469},
    'GTC': {'SF': 'SF4', 'aa': 'V', 'GC_rank': 1, 'N': 0, 'S': 0, 'C': 3, 'MW': 117.1469},
    'GTA': {'SF': 'SF4', 'aa': 'V', 'GC_rank': 0, 'N': 0, 'S': 0, 'C': 3, 'MW': 117.1469},
    'GTG': {'SF': 'SF4', 'aa': 'V', 'GC_rank': 1, 'N': 0, 'S': 0, 'C': 3, 'MW': 117.1469},

    #L
    'TTA': {'SF': 'SF6', 'aa': 'L', 'GC_rank': 1, 'N': 0, 'S': 0, 'C': 4, 'MW': 131.1736},
    'TTG': {'SF': 'SF6', 'aa': 'L', 'GC_rank': .5, 'N': 0, 'S': 0, 'C': 4, 'MW': 131.1736},
    'CTT': {'SF': 'SF6', 'aa': 'L', 'GC_rank': .5, 'N': 0, 'S': 0, 'C': 4, 'MW': 131.1736},
    'CTC': {'SF': 'SF6', 'aa': 'L', 'GC_rank': 0, 'N': 0, 'S': 0, 'C': 4, 'MW': 131.1736},
    'CTA': {'SF': 'SF6', 'aa': 'L', 'GC_rank': .5, 'N': 0, 'S': 0, 'C': 4, 'MW': 131.1736},
    'CTG': {'SF': 'SF6', 'aa': 'L', 'GC_rank': 0, 'N': 0, 'S': 0, 'C': 4, 'MW': 131.1736},

    #I
    'ATT': {'SF': 'SF3', 'aa': 'I', 'GC_rank': 1, 'N': 0, 'S': 0, 'C': 4, 'MW': 131.1736},
    'ATC': {'SF': 'SF3', 'aa': 'I', 'GC_rank': 0, 'N': 0, 'S': 0, 'C': 4, 'MW': 131.1736},
    'ATA': {'SF': 'SF3', 'aa': 'I', 'GC_rank': 1, 'N': 0, 'S': 0, 'C': 4, 'MW': 131.1736},

    #P
    'CCT': {'SF': 'SF4', 'aa': 'P', 'GC_rank': 1, 'N': 0, 'S': 0, 'C': 3, 'MW': 115.1310},
    'CCC': {'SF': 'SF4', 'aa': 'P', 'GC_rank': 0, 'N': 0, 'S': 0, 'C': 3, 'MW': 115.1310},
    'CCA': {'SF': 'SF4', 'aa': 'P', 'GC_rank': 1, 'N': 0, 'S': 0, 'C': 3, 'MW': 115.1310},
    'CCG': {'SF': 'SF4', 'aa': 'P', 'GC_rank': 0, 'N': 0, 'S': 0, 'C': 3, 'MW': 115.1310},

    #F
    'TTT': {'SF': 'SF2', 'aa': 'F', 'GC_rank': 1, 'N': 0, 'S': 0, 'C': 7, 'MW': 165.1900},
    'TTC': {'SF': 'SF2', 'aa': 'F', 'GC_rank': 0, 'N': 0, 'S': 0, 'C': 7, 'MW': 165.1900},

    #G
    'GGT': {'SF': 'SF4', 'aa': 'G', 'GC_rank': 1, 'N': 0, 'S': 0, 'C': 0, 'MW': 75.0669},
    'GGC': {'SF': 'SF4', 'aa': 'G', 'GC_rank': 0, 'N': 0, 'S': 0, 'C': 0, 'MW': 75.0669},
    'GGA': {'SF': 'SF4', 'aa': 'G', 'GC_rank': 1, 'N': 0, 'S': 0, 'C': 0, 'MW': 75.0669},
    'GGG': {'SF': 'SF4', 'aa': 'G', 'GC_rank': 0, 'N': 0, 'S': 0, 'C': 0, 'MW': 75.0669},

    #C
    'TGT': {'SF': 'SF2', 'aa': 'C', 'GC_rank': 1, 'N': 0, 'S': 1, 'C': 1, 'MW': 121.1590},
    'TGC': {'SF': 'SF2', 'aa': 'C', 'GC_rank': 0, 'N': 0, 'S': 1, 'C': 1, 'MW': 121.1590},

    #M
    'ATG': {'SF': 'SF1', 'aa': 'M', 'GC_rank': None, 'N': 0, 'S': 1, 'C': 3, 'MW': 149.2124},

    #K
    'AAA': {'SF': 'SF2', 'aa': 'K', 'GC_rank': 1, 'N': 1, 'S': 0, 'C': 4, 'MW': 146.1882},
    'AAG': {'SF': 'SF2', 'aa': 'K', 'GC_rank': 0, 'N': 1, 'S': 0, 'C': 4, 'MW': 146.1882},

    #W
    'TGG': {'SF': 'SF1', 'aa': 'W', 'GC_rank': None, 'N': 1, 'S': 0, 'C': 9, 'MW': 204.2262},

    #N
    'AAT': {'SF': 'SF2', 'aa': 'N', 'GC_rank': 1, 'N': 1, 'S': 0, 'C': 2, 'MW': 132.1184},
    'AAC': {'SF': 'SF2', 'aa': 'N', 'GC_rank': 0, 'N': 1, 'S': 0, 'C': 2, 'MW': 132.1184},

    #Q
    'CAA': {'SF': 'SF2', 'aa': 'Q', 'GC_rank': 1, 'N': 1, 'S': 0, 'C': 3, 'MW': 146.1451},
    'CAG': {'SF': 'SF2', 'aa': 'Q', 'GC_rank': 0, 'N': 1, 'S': 0, 'C': 3, 'MW': 146.1451},

    #H
    'CAT': {'SF': 'SF2', 'aa': 'H', 'GC_rank': 1, 'N': 2, 'S': 0, 'C': 4, 'MW': 155.1552},
    'CAC': {'SF': 'SF2', 'aa': 'H', 'GC_rank': 0, 'N': 2, 'S': 0, 'C': 4, 'MW': 155.1552},

    #R
    'CGT': {'SF': 'SF6', 'aa': 'R', 'GC_rank': .5, 'N': 3, 'S': 0, 'C': 4, 'MW': 174.2017},
    'CGC': {'SF': 'SF6', 'aa': 'R', 'GC_rank': 0, 'N': 3, 'S': 0, 'C': 4, 'MW': 174.2017},
    'CGA': {'SF': 'SF6', 'aa': 'R', 'GC_rank': .5, 'N': 3, 'S': 0, 'C': 4, 'MW': 174.2017},
    'CGG': {'SF': 'SF6', 'aa': 'R', 'GC_rank': 0, 'N': 3, 'S': 0, 'C': 4, 'MW': 174.2017},
    'AGA': {'SF': 'SF6', 'aa': 'R', 'GC_rank': 1, 'N': 3, 'S': 0, 'C': 4, 'MW': 174.2017},
    'AGG': {'SF': 'SF6', 'aa': 'R', 'GC_rank': .5, 'N': 3, 'S': 0, 'C': 4, 'MW': 174.2017},

    # stop codon table
    }


def calculate_SCU(gene, errorfile_handle, codon_dictionary = codon_dictionary):
    """
    This script takes a nucleotide sequence as input and returns measure of synonymous codon usage (SCU) Nc and
    average GC_rank (codons ranked by # nitrogen atoms relative to other synonymous codons) of the sequence.

    code_11 is a dictionary containing general information about each codon.
    """

    code_11 = {
        'SF_types': {
            'SF1': ['M', 'W'],
            'SF2': ['F', 'Y', 'H', 'Q', 'N', 'K', 'D', 'E', 'C'],
            'SF3': ['I'],
            'SF4': ['A', 'P', 'T', 'V', 'G'],
            'SF6': ['L', 'S', 'R']
        },  # complete sub-dictionary indexed by 'SF types'

        'aas': {
            'M': {'codons': ['ATG'], 'SF': 'SF1'},
            'W': {'codons': ['TGG'], 'SF': 'SF1'},
            'F': {'codons': ['TTT', 'TTC'], 'SF': 'SF2'},
            'Y': {'codons': ['TAT', 'TAC'], 'SF': 'SF2'},
            'H': {'codons': ['CAT', 'CAC'], 'SF': 'SF2'},
            'Q': {'codons': ['CAA', 'CAG'], 'SF': 'SF2'},
            'N': {'codons': ['AAT', 'AAC'], 'SF': 'SF2'},
            'K': {'codons': ['AAA', 'AAG'], 'SF': 'SF2'},
            'D': {'codons': ['GAT', 'GAC'], 'SF': 'SF2'},
            'E': {'codons': ['GAA', 'GAG'], 'SF': 'SF2'},
            'C': {'codons': ['TGT', 'TGC'], 'SF': 'SF2'},
            'I': {'codons': ['ATT', 'ATC', 'ATA'], 'SF': 'SF3'},
            'A': {'codons': ['GCT', 'GCC', 'GCA', 'GCG'], 'SF': 'SF4'},
            'P': {'codons': ['CCT', 'CCC', 'CCA', 'CCG'], 'SF': 'SF4'},
            'T': {'codons': ['ACT', 'ACC', 'ACA', 'ACG'], 'SF': 'SF4'},
            'V': {'codons': ['GTT', 'GTC', 'GTA', 'GTG'], 'SF': 'SF4'},
            'G': {'codons': ['GGT', 'GGC', 'GGA', 'GGG'], 'SF': 'SF4'},

            'L': {'codons': ['TTA', 'TTG', 'CTT','CTC','CTA','CTG'], 'SF': 'SF6'},
            'S': {'codons': ['TCT', 'TCC', 'TCA','TCG','AGT','AGC'], 'SF': 'SF6'},
            'R': {'codons': ['CGT', 'CGC', 'CGA', 'CGG','AGA','AGG'], 'SF': 'SF6'}
            }  # close aa information dictionary

        }  # close entire codon information dictionary

    # This dictionary tallies codon and aa counts for the gene
    codon_usage_dictionary = {'codons': {
        # Met, M
        'ATG': {'count': 0, 'frequency': 0},

        # Trp, W
        'TGG': {'count': 0, 'frequency': 0},

        # Phe, F
        'TTT': {'count': 0, 'frequency': 0},
        'TTC': {'count': 0, 'frequency': 0},

        # Tyr, Y
        'TAT': {'count': 0, 'frequency': 0},
        'TAC': {'count': 0, 'frequency': 0},

        # His, H
        'CAT': {'count': 0, 'frequency': 0},
        'CAC': {'count': 0, 'frequency': 0},

        # Gln,
        'CAA': {'count': 0, 'frequency': 0},
        'CAG': {'count': 0, 'frequency': 0},

        # Asn, N
        'AAT': {'count': 0, 'frequency': 0},
        'AAC': {'count': 0, 'frequency': 0},

        # Lys, K
        'AAA': {'count': 0, 'frequency': 0},
        'AAG': {'count': 0, 'frequency': 0},

        # Asp, D
        'GAT': {'count': 0, 'frequency': 0},
        'GAC': {'count': 0, 'frequency': 0},

        # Glu, E
        'GAA': {'count': 0, 'frequency': 0},
        'GAG': {'count': 0, 'frequency': 0},

        # Cys, C
        'TGT': {'count': 0, 'frequency': 0},
        'TGC': {'count': 0, 'frequency': 0},

       # Ile, I
        'ATT': {'count': 0, 'frequency': 0},
        'ATC': {'count': 0, 'frequency': 0},
        'ATA': {'count': 0, 'frequency': 0},

        # Ala, A
        'GCT': {'count': 0, 'frequency': 0},
        'GCC': {'count': 0, 'frequency': 0},
        'GCA': {'count': 0, 'frequency': 0},
        'GCG': {'count': 0, 'frequency': 0},

        # Pro, P
        'CCT': {'count': 0, 'frequency': 0},
        'CCC': {'count': 0, 'frequency': 0},
        'CCA': {'count': 0, 'frequency': 0},
        'CCG': {'count': 0, 'frequency': 0},

        # Thr, T
        'ACT': {'count': 0, 'frequency': 0},
        'ACC': {'count': 0, 'frequency': 0},
        'ACA': {'count': 0, 'frequency': 0},
        'ACG': {'count': 0, 'frequency': 0},

        # Val, V
        'GTT': {'count': 0, 'frequency': 0},
        'GTC': {'count': 0, 'frequency': 0},
        'GTA': {'count': 0, 'frequency': 0},
        'GTG': {'count': 0, 'frequency': 0},

        # Gly, G
        'GGT': {'count': 0, 'frequency': 0},
        'GGC': {'count': 0, 'frequency': 0},
        'GGA': {'count': 0, 'frequency': 0},
        'GGG': {'count': 0, 'frequency': 0},

        # Leu, L
        'TTA': {'count': 0, 'frequency': 0},
        'TTG': {'count': 0, 'frequency': 0},
        'CTT': {'count': 0, 'frequency': 0},
        'CTC': {'count': 0, 'frequency': 0},
        'CTA': {'count': 0, 'frequency': 0},
        'CTG': {'count': 0, 'frequency': 0},

        # Ser, S
        'TCT': {'count': 0, 'frequency': 0},
        'TCC': {'count': 0, 'frequency': 0},
        'TCA': {'count': 0, 'frequency': 0},
        'TCG': {'count': 0, 'frequency': 0},
        'AGT': {'count': 0, 'frequency': 0},
        'AGC': {'count': 0, 'frequency': 0},

        # Arg, R
        'CGT': {'count': 0, 'frequency': 0},
        'CGC': {'count': 0, 'frequency': 0},
        'CGA': {'count': 0, 'frequency': 0},
        'CGG': {'count': 0, 'frequency': 0},
        'AGA': {'count': 0, 'frequency': 0},
        'AGG': {'count': 0, 'frequency': 0},

        },
        'amino_acids': {
        'M': {'count': 0, 'frequency': 0, 'Ne': 0},
        'W': {'count': 0, 'frequency': 0, 'Ne': 0},
        'F': {'count': 0, 'frequency': 0, 'Ne': 0},
        'Y': {'count': 0, 'frequency': 0, 'Ne': 0},
        'H': {'count': 0, 'frequency': 0, 'Ne': 0},
        'Q': {'count': 0, 'frequency': 0, 'Ne': 0},
        'N': {'count': 0, 'frequency': 0, 'Ne': 0},
        'K': {'count': 0, 'frequency': 0, 'Ne': 0},
        'D': {'count': 0, 'frequency': 0, 'Ne': 0},
        'E': {'count': 0, 'frequency': 0, 'Ne': 0},
        'C': {'count': 0, 'frequency': 0, 'Ne': 0},
        'I': {'count': 0, 'frequency': 0, 'Ne': 0},
        'A': {'count': 0, 'frequency': 0, 'Ne': 0},
        'P': {'count': 0, 'frequency': 0, 'Ne': 0},
        'T': {'count': 0, 'frequency': 0, 'Ne': 0},
        'V': {'count': 0, 'frequency': 0, 'Ne': 0},
        'G': {'count': 0, 'frequency': 0, 'Ne': 0},
        'L': {'count': 0, 'frequency': 0, 'Ne': 0},
        'S': {'count': 0, 'frequency': 0, 'Ne': 0},
        'R': {'count': 0, 'frequency': 0, 'Ne': 0},
        },

        'SF_type': {
        'absent': [],
        'SF1': [],
        'SF2': [],
        'SF3': [],
        'SF4': [],
        'SF6': []
        },  # complete sub-dictionary indexed by 'SF types'

        'gene_codon_length': 0,
        'number_of_positions_with_GC_variability': 0,
        'sum_GC_rank': 0
        }

    # gene parsing begins here.

    # remove start and stop codon from current DNA sequence
    current_gene_sequence = gene.seq[3:(len(gene.seq)-3)]

    # Parse through codons
    for i in xrange(0, len(current_gene_sequence), 3):
        # count each codon
        current_codon = current_gene_sequence[i:(i+3)]

        #identify possible codons that will cause problems
        if current_codon in ['TAA', 'TAG', 'TGA']:
            error1 = gene.id.strip() + '\t' + 'internal stop codon' + '\t' + current_codon + '\n'
            errorfile_handle.writelines(error1)
            return ['Nan', 'Nan']

        if 'N' in current_codon:
            error2 = gene.id.strip() + '\t' + 'N present' + '\t' + current_codon + '\n'
            errorfile_handle.writelines(error2)
            return ['Nan', 'Nan']

        if current_codon not in codon_usage_dictionary['codons']:
            print current_codon, 'hu', codon_usage_dictionary['codons'][current_codon]
            error3 = gene.id.strip() + '\t' + 'other issue' + '\t' + current_codon + '\n'
            errorfile_handle.writelines(error3)
            return ['Nan', 'Nan']

        codon_usage_dictionary['codons'][current_codon]['count'] += 1
        codon_usage_dictionary['amino_acids'][codon_dictionary[current_codon]['aa']]['count'] += 1
        codon_usage_dictionary['gene_codon_length'] += 1

        # count the number of codons with GC variability and their rank
        if codon_dictionary[current_codon]['GC_rank'] != None:
            codon_usage_dictionary['number_of_positions_with_GC_variability'] += 1
            codon_usage_dictionary['sum_GC_rank'] += codon_dictionary[current_codon]['GC_rank']

    # calculate codon frequencies and update codon_usage_dictionary
    for k in codon_usage_dictionary['codons'].keys():
        amino_acid_total_usage = float(codon_usage_dictionary['amino_acids'][codon_dictionary[k]['aa']]['count'])

        if amino_acid_total_usage > 1:
            codon_usage_dictionary['codons'][k]['frequency'] = float(codon_usage_dictionary['codons'][k]['count']) \
                                                              / amino_acid_total_usage

        #if amino acid is rarely used (aka one or less times in protein)
        if amino_acid_total_usage < 2:
           codon_usage_dictionary['codons'][k]['frequency'] = 'Rare'
           #codon_usage_dictionary['amino_acids'][codon_dictionary[k]['aa']]['Ne'] == 'absent'

    # Calculate Ne for each aa and record in codon_usage_dictionary['amino_acids']['aa']['Ne']
    for aa in code_11['aas'].keys():
        if codon_usage_dictionary['amino_acids'][aa]['count'] < 2:
            codon_usage_dictionary['amino_acids'][aa]['Ne'] = 'absent'
            codon_usage_dictionary['SF_type']['absent'].append(aa)


        elif codon_usage_dictionary['amino_acids'][aa]['count'] > 1:
            squared_fequencies = []

            for each_codon in code_11['aas'][aa]['codons']:
                p = float(codon_usage_dictionary['codons'][each_codon]['frequency'])
                squared_fequencies.append(p*p)

            n = codon_usage_dictionary['amino_acids'][aa]['count']
            F = (((n * (sum(squared_fequencies)))-1) * 1.0/(n-1.0))


            # amino acid is rarely used if the numerator or denominator of equation for F == 0
            #
            if (((n * (sum(squared_fequencies)))-1) * 1.0/(n-1.0)) == 0:
                codon_usage_dictionary['amino_acids'][aa]['Ne'] = 'absent'
                codon_usage_dictionary['SF_type']['absent'].append(aa)

            if round(F, 10) != 0:
                Ne = 1 / F
                codon_usage_dictionary['amino_acids'][aa]['Ne'] = Ne
                codon_usage_dictionary['SF_type'][code_11['aas'][aa]['SF']].append(aa)


        #print aa, n, F, Ne, squared_fequencies, codon_usage_dictionary['SF_type']

    # calculate GC rank sums
    GC_variability_summary = (float(codon_usage_dictionary['sum_GC_rank'])\
                              / float(codon_usage_dictionary['number_of_positions_with_GC_variability'])\
                             if codon_usage_dictionary['number_of_positions_with_GC_variability'] != 0 else 'Nan')

    # parse through SF values to get final Nc
    Nc = 2.0
    #print gene.id

    for sf_type in ['SF2', 'SF4', 'SF6']:
        aas_for_sf_type = codon_usage_dictionary['SF_type'][sf_type]
        sf_float = float(len(aas_for_sf_type))

        # if one of the groups other than SF3 isn't represented, then dont return an Nc value
        if sf_float == 0:
            return(['Nan', str(GC_variability_summary)])

        #This is actually F list
        F_list = [1/codon_usage_dictionary['amino_acids'][x]['Ne'] for x in aas_for_sf_type]
        av_F_for_sf_type = sum(F_list) / sf_float
###
        if sf_type == 'SF2':
            Nc += 9/av_F_for_sf_type
        if sf_type == 'SF4':
            Nc += 5/av_F_for_sf_type
        if sf_type == 'SF6':
            Nc += 3/av_F_for_sf_type

        codon_usage_dictionary[''.join(sf_type + '_Ne')] = av_F_for_sf_type
              
###
    #SF3 is a special case because only I is SF3
    if len(codon_usage_dictionary['SF_type']['SF3']) == 1:
        Nc += codon_usage_dictionary['amino_acids']['I']['Ne']

    if len(codon_usage_dictionary['SF_type']['SF3']) == 0:
        Nc += 1/((codon_usage_dictionary['SF2_Ne'] + codon_usage_dictionary['SF4_Ne']) / float(2))
                
    if Nc > 61:
       Nc = 61
    return [str(round(Nc, 4)), str(GC_variability_summary)] #, codon_usage_dictionary, code_11


def codon_grouper(n, iterable, fillvalue=None):
    """grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"""
    args = [iter(iterable)] * n
    return itertools.izip_longest(fillvalue=fillvalue, *args)


def screen_out_ambiguous_codons(sequence):
    #removes codons with Ns
    return reduce(add, [reduce(add, x) for x in codon_grouper(3, sequence) if len(set(x) - set(['A', 'T', 'C', 'G'])) < 1])


def random_permutation(iterable):
    """Random selection from itertools.permutations(iterable, r)"""
    pool = tuple(iterable)
    r = len(iterable)
    # check for a stop codon. if there is one, pick a new sequence.
    random_sequence = reduce(add, random.sample(pool, r))
    for codon in codon_grouper(3, random_sequence):
        codon = reduce(add, codon)
        if codon in ['TAA', 'TAG', 'TGA']:
            return random_permutation(iterable)
    return random_sequence

def return_coding_sequence(fna_handle):
    """creates a generator that cycles through CDS features in gbff files"""
    x=0
    g = SeqIO.parse(fna_handle,'genbank')
    while True:
        genome = g.next()
        for gene in genome.features:
            if gene.type != "CDS":
                continue
            if 'pseudo' in gene.qualifiers:
                continue
            if 'exception' in gene.qualifiers:
                continue
            is_entire_gene = gene.location.__repr__()
            if 'BeforePosition' in is_entire_gene:
                continue
            if 'AfterPosition' in is_entire_gene:
                continue
            sequence = gene.extract(genome.seq)
            if gene.qualifiers['codon_start'][0] != '1':
                print gene.qualifiers['codon_start'][0]
                first_codon = int(gene.qualifiers['codon_start'][0]) - 1
                sequence = sequence[first_codon:]
                print gene.qualifiers['locus_tag'][0]
            gene_seq = SeqRecord(sequence, id = gene.qualifiers['locus_tag'][0] + ' ' + gene.qualifiers['product'][0])
            yield gene_seq


def scan_for_stop_codons(DNA_sequence):
    for codon in codon_grouper(3, DNA_sequence[0:-3]):
        codon = reduce(add, codon)
        if reduce(add, codon) in ['TAA', 'TAG', 'TGA']:
            return 'True'
    return False

def ARSC_and_MW_from_amino_acids(protein_sequence):
    '''
    This functions takes an amino acid sequence coded in single letters and returns N/C ARSC and Molecular Weight
    N and S counts come from page 30 of 'Understanding Bioinformatics' by Zvelbil and Baum
    molecular weights from http://www.webqc.org/aminoacids.php
    '''

    # create dictionary for amino acids, MW = MW of entire aa, 'N', 'S' counts only include side chains
    aa_dictionary = {
        'K': {'N': 1, 'S': 0, 'MW': 146.1882, 'C': 4},
        'R': {'N': 3, 'S': 0, 'MW': 174.2017, 'C': 4},
        'H': {'N': 2, 'S': 0, 'MW': 155.1552, 'C': 4},
        'D': {'N': 0, 'S': 0, 'MW': 133.1032, 'C': 2},
        'E': {'N': 0, 'S': 0, 'MW': 147.1299, 'C': 3},
        'N': {'N': 1, 'S': 0, 'MW': 132.1184, 'C': 2},
        'Q': {'N': 1, 'S': 0, 'MW': 146.1451, 'C': 3},
        'S': {'N': 0, 'S': 0, 'MW': 105.0930, 'C': 1},
        'T': {'N': 0, 'S': 0, 'MW': 119.1197, 'C': 2},
        'Y': {'N': 0, 'S': 0, 'MW': 181.1894, 'C': 7},
        'A': {'N': 0, 'S': 0, 'MW': 89.0935,  'C': 1},
        'V': {'N': 0, 'S': 0, 'MW': 117.1469, 'C': 3},
        'L': {'N': 0, 'S': 0, 'MW': 131.1736, 'C': 4},
        'I': {'N': 0, 'S': 0, 'MW': 131.1736, 'C': 4},
        'P': {'N': 0, 'S': 0, 'MW': 115.1310, 'C': 3},
        'F': {'N': 0, 'S': 0, 'MW': 165.1900, 'C': 7},
        'M': {'N': 0, 'S': 1, 'MW': 149.2124, 'C': 3},
        'W': {'N': 1, 'S': 0, 'MW': 204.2262, 'C': 9},
        'G': {'N': 0, 'S': 0, 'MW': 75.0669,  'C': 0},
        'C': {'N': 0, 'S': 1, 'MW': 121.1590, 'C': 1},
        'U': {'N': 0, 'S': 0, 'MW': 168.07,   'C': 1},
        'J': {'N': 0, 'S': 0, 'MW': 131.1736, 'C': 4},
        'B': {'N': 0.5, 'S': 0, 'MW': 132.6108, 'C': 2},
        'Z': {'N': 0.5, 'S': 0, 'MW': 146.6375, 'C': 3}
    }

    # remove whitespaces and '*' (termination) that prodigal adds to end of aa sequences
    protein_sequence = protein_sequence.strip().strip('*')

    # remove ambigous aas from the string
    protein_sequence_no_Xs = ''.join([x for x in protein_sequence if x != 'X' if x != '-'])

    # caculate ARSC and Molecular Weight
    N_ARSC = sum(map(lambda x: aa_dictionary[x]['N'], protein_sequence_no_Xs)) / float(len(protein_sequence_no_Xs))
    C_ARSC = sum(map(lambda x: aa_dictionary[x]['C'], protein_sequence_no_Xs)) / float(len(protein_sequence_no_Xs))
    av_molecular_weight = sum(map(lambda x: aa_dictionary[x]['MW'], protein_sequence_no_Xs))/float(len(protein_sequence_no_Xs))

    return[str(round(N_ARSC, 4)), str(round(av_molecular_weight, 4)), str(round(C_ARSC, 4))]

def ARSC_MW_from_nucleotides(sequence):
    '''
    This functions takes a protein sequence coded in nucleotides and returns N-ARSC.
    Internal stop codons will cause this script problems
    '''

    # remove whitespaces and stop codon at end of sequences
    sequence = sequence.strip()[0:-3]

    # caculate NARSC
    aa_sequence_length = len(sequence) / 3.0
    total_nitrogen_atoms = sum(map(lambda x: codon_dictionary[reduce(add, x)]['N'], codon_grouper(3, sequence)))
    av_ARSC = total_nitrogen_atoms / aa_sequence_length

    total_carbon_atoms = sum(map(lambda x: codon_dictionary[reduce(add, x)]['C'], codon_grouper(3, sequence)))
    av_C_ARSC = total_carbon_atoms / aa_sequence_length

    total_molecular_weight = sum(map(lambda x: codon_dictionary[reduce(add, x)]['MW'], codon_grouper(3, sequence)))
    av_molecular_weight = total_molecular_weight / aa_sequence_length

    return str(round(av_ARSC, 4)), str(round(av_molecular_weight, 4)), str(round(av_C_ARSC, 4)), str(round(av_ARSC/av_C_ARSC, 4))




import os, re, subprocess, shlex
import numpy as np
from Bio import SeqUtils
from Bio.Seq import Seq

genome_dir = "/home/frankaylward/Desktop/marker_gene_benchmarking/Marinimicrobia/all_genomes/genomes/"
protein_dir = "/home/frankaylward/Desktop/marker_gene_benchmarking/Marinimicrobia/all_genomes/proteins/"

gc_file = open("GC_stats.txt", "w")
gc_file.write("Genome\tGC\tSize\n")

arsc_file = open("ARSC_stats.txt", "w")
arsc_file.write("Genome\tNARSC\tCARSC\n")

for i in os.listdir(genome_dir):
	if i.endswith(".fna") or i.endswith(".fa"):
		print i

		full_record = Seq("")
		input_file = os.path.join(genome_dir, i)
		for record in SeqIO.parse(input_file, "fasta"):
			seq = record.seq.strip().strip('*')
			full_record = full_record + seq

		# get final gc
		final_gc = SeqUtils.GC(full_record)
		final_length = len(full_record)
		print i, final_gc, final_length
		gc_file.write(i +"\t"+ str(final_gc) +"\t"+ str(final_length) +"\n")


for i in os.listdir(protein_dir):
	if i.endswith(".faa"):
		protein_file = os.path.join(protein_dir, i)
		proteins = open(protein_file, "r")

		full_record = Seq("")
		for record in SeqIO.parse(proteins, "fasta"):
			seq = record.seq.strip().strip('*')
			full_record = full_record + seq

		#print(full_record.seq)
		#full_record = re.sub("J", "I", full_record)

		full_arsc = ARSC_and_MW_from_amino_acids(full_record)
		full_narsc = full_arsc[0]
		full_carsc = full_arsc[2]
		print i, full_narsc, full_carsc
		arsc_file.write(i +"\t"+ str(full_narsc) +"\t"+ str(full_carsc) +"\n")

	
	
	
	
	
	
	
	
