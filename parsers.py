#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Parse the inpot arguments
date: 2022-10-03
'''
__author__ = 'Yan Li'
__email__= 'yan.li2@outlook.com'

import sys, os
import argparse

CURRENT_PATH = os.path.dirname(__file__)
DB_PATH = CURRENT_PATH+'/db/'

# Parse the user input arguments
def get_arguments():

    parser = argparse.ArgumentParser(description="Salmonella Whole genome CRISPR Typing pipeline: Detect CRISPR1 and CRISPR2 alleles of Salmonella genomes")

    # Input & output options
    io_group = parser.add_argument_group("Input/Output")
    io_group.add_argument("-i", "--input", required=True, metavar="fasta",
                            help="Input FASTA file (required)")
    io_group.add_argument("-o", "--output", required=True, metavar="filename",
                            help="Output file base names (required)")
    io_group.add_argument('-O', '--outdir', required=False, metavar='path',
                              default='.',
                              help='Output directory')
    
    # Other options
    other_group = parser.add_argument_group('Other')
    other_group.add_argument('-d', '--database_path', metavar='db_path',
                             default=DB_PATH,
                             help='BLAST databases of Salmonella CRISPR DR and spacers, default={}'.format(DB_PATH))
    other_group.add_argument('--algorithm', choices=['blast', 'minced'],
                            default = 'blast',
                            help='The algorithm used for searching CRISPR alleles. They may generate different results. Default="blast". "minced" is an option when the "blast" does not work.')
    
    # If no arguments were used, print the entire help
    if len(sys.argv) == 1:
        parser.print_help(file=sys.stderr)
        sys.exit(1)
    
    # Go
    args = parser.parse_args()

    if args.outdir == None:
        args.outdir = '.'
    
    return(args)