#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Get CRISPR position from BLAST result
date: 2022-09-29
'''
__author__ = 'Yan Li'
__email__= 'yan.li2@outlook.com'

import sys
import pandas as pd

MAX_SPACER_LENGTH = 150 # allow to skip one or two DR

def get_crispr_alleles(blast_dr_result_file):
    """
    Read the DR Blast result, write CRISPR alleles and direct repeats to a GFF file.
    Args:
        blast_dr_result_file (str): File name of a sorted BLAST result
    """

    blast_out = pd.read_csv(blast_dr_result_file, sep='\t', names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'slen', 'qseq'])
    
    # Decide strand direction
    contig_direction = {}
    print("# Decide strand direction:")
    for contig in pd.unique(blast_out['qseqid']):
        strand_direction_votes = []
        blast_out_per_contig = blast_out[blast_out['qseqid']==contig]
        strand_direction_vote_values = (blast_out_per_contig['send'] - blast_out_per_contig['sstart']).tolist()
        strand_direction_votes = ['Forward' if v>0 else 'Backward' for v in strand_direction_vote_values]
        strand_direction = max(strand_direction_votes, key=strand_direction_votes.count)
        print("The direction of contig {0} is {1}".format(contig, strand_direction))
        contig_direction[contig] = strand_direction

    # Find CRISPR allele start and end


    # If overlap, use the higher bitscore


    # Decide wether it's a new variant (pident == 100)


    # Out put Gff file

    return

'''

def select_true_hit(blast_out):
    """
    Select the hit that 100% identical to the DR database

    Args:
        blast_out (dataframe): Blast result that read with pandas

    Returns:
        (dataframe) Select Blast result
    """

    blast_out.sort_values(by=['qseqid', 'qstart'], inplace=True, ignore_index=True)

    rows_to_return = []
    
    for index, row in blast_out.iterrows():
    # Assume all the query records are always 'qstart' < 'qend'



    return(selected_blast_out)
'''


if __name__ == '__main__':
    get_crispr_alleles(sys.argv[1])