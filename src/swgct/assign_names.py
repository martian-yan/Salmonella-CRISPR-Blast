#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Assign CT names according to the spacer array
date: 2022-10-20
'''
__author__ = 'Yan Li'
__email__= 'yan.li2@outlook.com'

import pandas as pd
import sys

def main(summary_file, ct_file):

    summary = pd.read_csv(summary_file, sep='\t')
    ct = pd.read_csv(ct_file, sep='\t')

    for index, row in summary.iterrows():
        for i, row_ct in ct.iterrows():
            if (row['CRISPR1'] == row_ct['CRISPR1']) & (row['CRISPR2'] == row_ct['CRISPR2']):
                summary.at[index, 'CT'] = row_ct['name']


    summary.to_csv(summary_file[:-3]+'named.tsv', sep='\t', index=False)
    return

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
