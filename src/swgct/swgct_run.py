#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
The MAIN module of Salmonella Whole genome CRISPR Typing pipeline
date: 2022-10-03
'''
__author__ = 'Yan Li'
__email__ = 'yan.li2@outlook.com'

import os
import shutil

from swgct.get_crispr_alleles import blast_crispr_dr, sort_crispr_alleles, annotate_crispr
from swgct.misc import random_string


def run(args):

        
    if args.outdir == None:
        args.outdir = '.'.join(args.input.split('.')[:-1])
    if args.prefix == None:
        args.prefix = os.path.basename('.'.join(args.input.split('.')[:-1]))

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
    # Processing the raw reads
    # TODO: add mince
    # if args.algorithm == 'blast'
    blast_crispr_dr(fasta_file, dr_db, out_dir, out_name, tmpdir)
    # if args.algorithm == 'minced'
    sort_crispr_alleles(fasta_file, out_dir, out_name, tmpdir)
    annotate_crispr(fasta_file, spacer_db, out_dir, out_name, tmpdir)

    # remove tmp dir
    if not args.keep:
        shutil.rmtree(tmpdir)

    return

