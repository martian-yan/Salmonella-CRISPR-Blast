#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Summarize the result from SWGCT, do:
1. summarize new spacers and new spacer variants
2. make a table of CRISPR Types
3. make a table of presence/absence of spacers, which could make a NJ tree

Note: The result folders are expected to be under the same directory, such as:
all_results
|-- result_1/
    |-- result_1.CRISPR.fasta
    |-- result_1.CRISPR.gff
    |-- result_1.CRISPR.tsv
    |-- ...
|-- result_2/
    |-- ...
|-- result_3/
    |-- ...
|-- ...

ct_summary is expected to run at the dir `all_results`

date: 2022-10-19
'''
__author__ = 'Yan Li'
__email__= 'yan.li2@outlook.com'

import pandas as pd
import numpy as np
import os
import sys
import subprocess
from Bio import SeqIO

def summary(args):

    work_dir = args.dir.rstrip('/')
    # Check results
    all_files = []
    for root, dirs, files in os.walk(work_dir):
        for dirname in files:
            all_files.append(os.path.join(root, dirname))
 
    swgct_output = [f for f in all_files if f.endswith('.CRISPR.tsv')]
    swgct_output.sort()

    if not swgct_output:
        exit("CANNOT find any swgct results in current working directory")

    print("Found {n} results from \"swgct run\":\n".format(n=len(swgct_output)))
    for n in swgct_output[0:5]:
        print(n)
    if len(swgct_output) > 5:
        print("...")
    print("Start summarize...\n")

    # Step 1: summarize new spacers and new spacer variants
    # 
    find_new_spacers("new_spacer", all_files, work_dir)
    find_new_spacers("new_spacer_var", all_files, work_dir)

    # Step 2: make a table of CRISPR
    # Find all files
    crispr_summary_list = []
    crispr_summary = pd.DataFrame(columns=['CRISPR1', 'CRISPR2', 'unsorted'])
    for tsv in swgct_output:
        filename = os.path.basename(tsv)
        filename = filename.split(".")[0]
        with open(tsv, "r") as f:
            read_lines = f.read().splitlines()
            crispr_type_dict = {}
            unsorted = []
            for line in read_lines:
                allele = line.split('\t')[0]
                array = line.split('\t')[1]
                if allele == "CRISPR1" or allele == "CRISPR2":
                    crispr_type_dict[allele] = array
                if allele == "CRISPR":
                    unsorted.append(array)

            if unsorted:
                crispr_type_dict['unsorted'] = unsorted

            crispr_type = pd.Series(crispr_type_dict, name=filename, dtype=str)
            crispr_summary_list.append(crispr_type)
    
    crispr_summary_add = pd.DataFrame(crispr_summary_list)
    crispr_summary = pd.concat([crispr_summary, crispr_summary_add]) # Fix the bug that the column "unsorted" may not present by using a empty df with column names 

    print("\nCreated CRISPR summary table:")
    print(crispr_summary)
    crispr_summary.to_csv(work_dir+'/crispr_summary.unsorted.tsv', sep='\t')

    unsorted_before = len([x for x in list(crispr_summary['unsorted']) if not pd.isna(x)])

    # Put unsorted crispr to CRISPR1 or CRISPR2 if it matches any of the sorted CRISPR array
    # Potential BUG: If any weird assembly have reversed CRISPR1 and CRISPR2 positions, all the unsorted results cannot be sorted correctly.
    # TODO: A cross comparing of CRISPR1 and CRISPR2 - see if any CRISPR2 are in CRISPR1 column
    for index, row in crispr_summary.iterrows():
        if (not pd.isna(row['unsorted'])) and (len(row['unsorted']) == 2):
            # If more than 2 CRISPR alleles, they can be fragments, thus cannot be sorted
            return_value = list(row['unsorted']) # crispr_summary.at[index, 'unsorted'] somehow turned into tuple by pandas
            for crispr in row['unsorted']:
                if crispr in list(crispr_summary['CRISPR1']):
                    crispr_summary.at[index, 'CRISPR1'] = crispr
                    return_value.remove(crispr) 
                elif crispr in list(crispr_summary['CRISPR2']):
                    crispr_summary.at[index, 'CRISPR2'] = crispr
                    return_value.remove(crispr)
            if len(return_value) == 0:
                crispr_summary.at[index, 'unsorted'] = np.nan
            else:
                crispr_summary.at[index, 'unsorted'] = tuple(return_value)

    unsorted_after = len([x for x in list(crispr_summary['unsorted']) if not pd.isna(x)])

    print("\nThe CRISPR alleles of {} strains were sorted according to its array".format(unsorted_before-unsorted_after))
    print("Remaining unsorted CRISPRs:", len([x for x in list(crispr_summary['unsorted']) if not pd.isna(x)]))
    crispr_summary.to_csv(work_dir+"/CRISPR_summaries.tsv", sep='\t')
    # find new spacers
    # TODO: assign new spacer names to the CRISPR results

    # Step 3 Make a table of uniq_CRISPR_Types
    crispr_summary_sorted = crispr_summary[['CRISPR1', 'CRISPR2']]
    uniq_ct = crispr_summary_sorted.value_counts().to_frame().reset_index() # df.value_counts generate a Series which could turn into df
    uniq_ct.to_csv(work_dir+"/uniq_CRISPR_Types.tsv", sep='\t')

    # Step 4 make a table of presence and absence of each spacer
    crispr1_matrix = get_spacer_matrix(uniq_ct['CRISPR1'])
    crispr2_matrix = get_spacer_matrix(uniq_ct['CRISPR2'])
    crispr_matrix = pd.concat([crispr1_matrix, crispr2_matrix], axis=1)
    crispr_matrix.to_csv(work_dir+'/crispr_matrix.tsv', sep='\t')

    return

def find_new_spacers(mode, all_files, work_dir):
    """
    Find new spacers or spacer variants

    Args:
        mode (str): one of "new_spacer", "new_spacer_var"
        all_files (list): all possible files
        work_dir (str): current working directory
    """    

    if mode == "new_spacer":
        msg = "new spacers"
    else:
        msg = "new spacer variants"

    new_spacer_results = [f for f in all_files if f.endswith('{}.fasta'.format(mode))]

    if not new_spacer_results:
        print("\nNo {}".format(msg))
    else:
        cmd = "cat {0} | seqkit rmdup -s | seqkit rename | seqkit sort > {1}/{2}s.fasta".format(' '.join(new_spacer_results), work_dir, mode)
        print(cmd)
        return_value = subprocess.run(cmd, shell=True)
        
        if return_value.returncode != 0:
            print("\nNo {}".format(msg))
            os.remove("{0}/{1}s.fasta".format(work_dir, mode))
        else:
            cmd = "grep '>' {}s.fasta | wc".format(mode)
            return_value = subprocess.run(cmd, shell=True, capture_output=True)
            spacer_num = return_value.stdout.decode().lstrip().split()[0]
            print("\nFound {0} {1}".format(spacer_num, msg))
    
    return

def get_spacer_matrix(col):

    col = col.apply(lambda x: x.split('-'))

    spacers = [spacer for sublist in col for spacer in sublist] # flatten a list of lists
    spacers = list(set(spacers))
    spacers.sort()


    df_to_return = pd.DataFrame(0, columns=spacers, index=col.index)
    for index, spacer_array in col.items():
        spacer_array = list(set(spacer_array)) # remove duplicated spacers, especially labeled "new_spacers"
        df_to_return.loc[index, spacer_array] = 1

    return df_to_return