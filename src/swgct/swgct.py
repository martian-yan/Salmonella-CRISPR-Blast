#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
The MAIN module of Salmonella Whole genome CRISPR Typing pipeline
date: 2022-10-03
'''
__author__ = 'Yan Li'
__email__ = 'yan.li2@outlook.com'
__version__ = '1.0.0'

from .parsers import get_arguments
from .get_crispr_alleles import blast_crispr_dr, sort_crispr_alleles, annotate_crispr
from .misc import random_string

import os
import shutil

def main():

    args = get_arguments()
    # Processing the raw reads
    # TODO: add mince
    
    ## Input
    fasta_file = args.input
    dr_db = args.database_path.rstrip('/') + '/DR_Salmonella'
    spacer_db = args.database_path.rstrip('/') + '/spacers_Salmonella'

    ## Output
    out_dir = args.outdir.rstrip('/')
    out_name = args.prefix
    # creat tmp dir
    tmpdir = "{outdir}/{tmp_dir}".format(outdir=out_dir,
                                         tmp_dir="tmp_" + random_string())
    os.makedirs(tmpdir)

    # Run
    # if args.algorithm == 'blast'
    blast_crispr_dr(fasta_file, dr_db, out_dir, out_name, tmpdir)
    # if args.algorithm == 'minced'
    sort_crispr_alleles(fasta_file, out_dir, out_name, tmpdir)
    annotate_crispr(fasta_file, spacer_db, out_dir, out_name, tmpdir)

    # remove tmp dir
    if not args.keep:
        shutil.rmtree(tmpdir)

    return


if __name__ == "__main__":
    main()