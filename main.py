#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
The MAIN module of Salmonella Whole genome CRISPR Typing pipeline
date: 2022-10-03
'''
__author__ = 'Yan Li'
__email__ = 'yan.li2@outlook.com'

from parsers import get_arguments
from get_crispr_alleles import blast_crispr_dr, sort_crispr_alleles


def main():

    args = get_arguments()
    # Processing the raw reads
    # TODO: add mince
    # if args.algorithm == 'blast'
    blast_crispr_dr(args)
    # if args.algorithm == 'minced'
    sort_crispr_alleles(args)

    return


if __name__ == "__main__":
    main()