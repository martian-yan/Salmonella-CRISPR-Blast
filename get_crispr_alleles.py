#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Get CRISPR position from BLAST result
date: 2022-09-29
'''
__author__ = 'Yan Li'
__email__ = 'yan.li2@outlook.com'

import logging
import os
import shutil

import pandas as pd

from misc import cmd, current_time, random_string, write_log

MIN_REPEATS = 2  # define CRISPR allele with at least 2 DR in a row
MAX_SPACER_LENGTH = 150  # allow to skip one or two DR
MAX_OVERLAPPING = 1  # dont allow overlaps of DRs


def blast_crispr_dr(ARGS):
    """
    Read the DR Blast result, write CRISPR alleles and direct repeats to a GFF file.
    Output files:
        log_file = out_file+'.log',
        sorted_blast_out = out_file+'.blast_DR.sorted.tsv',
        gff_file = out_file+'.DR_blast.gff'
    """
    ## Input
    fasta_file=ARGS.input,
    dr_db=ARGS.database_path.rstrip('/') + '/DR_Salmonella'

    ## Output
    outdir = ARGS.outdir.rstrip('/')
    out_name = ARGS.output
    out_file = "{outdir}/{out_name}".format(outdir=outdir, out_name=out_name)

    # If Outdir not exist, create one
    outdir_exist = os.path.exists(outdir)
    if not outdir_exist:
        os.mkdir(outdir)

    ## Log
    log_file = out_file + '.log'
    # Start log
    screen_log = "{time}: START looking for CRISPR in {fasta}".format(time=current_time(), fasta=ARGS.input)
    write_log(screen_log, log_file)

    ## Step 1, BLAST against the DR database
    sorted_blast_out = out_file + '.blast_DR.sorted.tsv'
    blast_cmd = "blastn -db {dr_db} -query {fasta_file} -out {blast_out} -outfmt '6 std slen qseq' -perc_identity 90 -culling_limit 5; awk '$4==$13' {blast_out} | sort -k 7 -n > {sorted_blast_out}".format(
        fasta_file=fasta_file,
        dr_db=dr_db,
        blast_out=out_file + '.blast_DR.tsv',
        sorted_blast_out=sorted_blast_out,
    )
    cmd(blast_cmd, "Blast DR sequences", log_file)

    ## Step 2, get 2 CRISPR allele positions
    # Print log
    screen_log = "Looking for CRISPR alleles"
    write_log(screen_log, log_file)

    blast_out = pd.read_csv(sorted_blast_out,
                            sep='\t',
                            names=[
                                'qseqid', 'sseqid', 'pident', 'length',
                                'mismatch', 'gapopen', 'qstart', 'qend',
                                'sstart', 'send', 'evalue', 'bitscore', 'slen',
                                'qseq'
                            ])

    # Remove overlaps - Find perfect match or use the higher bitscore
    blast_out = rm_overlaps(blast_out)

    # Find CRISPR allele start and end
    # Crispr_alleles is a list of DataFrame then output as GFF
    crispr_alleles = find_crispr_alleles(blast_out)

    # Out put Minced style Gff file
    if crispr_alleles == []:
        screen_log = "\nCould not find any CRISPR allele"
        write_log(screen_log, log_file, level=logging.CRITICAL)
        quit()
    else:
        screen_log = "\nFound {n} CRISPR allele(s)".format(
            n=len(crispr_alleles))
        write_log(screen_log, log_file)

    gff_records = []
    for crispr_allele in crispr_alleles:
        gff_records += blast2gff(crispr_allele)

    gff_file = out_file + '.DR_blast.gff'
    with open(gff_file, 'w') as fo:
        fo.write('# GFF file created from blast result {}\n'.format(
            sorted_blast_out))
        fo.write('# BLAST cmd: {}\n'.format(blast_cmd))
        fo.write('# The float in column 6 is identity%\n')
        fo.write('#\n')
        for record in gff_records:
            fo.write('{}\n'.format('\t'.join(record)))

    return


def sort_crispr_alleles(ARGS):
    """
    1. Decide CRISPR1 and CRISPR2 sites
    2. Output GFF and a simple table with CRISPR arrays
    """    

    # input
    outdir = ARGS.outdir.rstrip('/')
    out_name = ARGS.output
    out_file = "{outdir}/{out_name}".format(outdir=outdir, out_name=out_name)
    dr_gff = out_file + '.DR_blast.gff'

    # output
    log_file = out_file + '.log'
    # creat tmp dir
    tmpdir = "{outdir}/{tmp_dir}".format(outdir=outdir, tmp_dir="tmp_"+random_string())

    crispr_alleles = pd.read_csv(dr_gff, sep='\t', header=None, comment='#', names=['qseqid', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes'])
    crispr_alleles = crispr_alleles[crispr_alleles['feature'] == 'CRISPR']
    crispr_sorted = False

    strand = "Unclear"
    sort_ascending = True
    if all(strand == '+' for strand in crispr_alleles['strand']):
        strand = "Forward"
    elif all(strand == '-' for strand in crispr_alleles['strand']):
        strand = "Backward"
        sort_ascending = False
    
    crispr_alleles.sort_values(by=['qseqid', 'qstart'],
                          inplace=True,
                          ascending= [True, sort_ascending],
                          ignore_index=True)

    # TODO: Possible upgrade - train a ML model to dicide whether it's CRISPR1 or 2
    # Step 1 decide CRISPR1 and 2 through the order if they on the SAME contig
    if (len(crispr_alleles) == 2) and (strand != "Unclear") and all(qseqid == crispr_alleles.loc[0, 'qseqid'] for qseqid in crispr_alleles['qseqid']):
        crispr_alleles.loc[0, 'feature'] = 'CRISPR1'
        crispr_alleles.loc[1, 'feature'] = 'CRISPR2'
        crispr_sorted = True
        
    
    screen_log = "CRISPR allele sorted according to their position:\n\
                CRISPR1 is at: {0}..{1}\n\
                CRISPR2 is at: {2}..{3}".format(crispr_alleles.loc[0, 'start'], crispr_alleles[0, 'end'],
                                               crispr_alleles.loc[1, 'start'], crispr_alleles[1, 'end'])
    write_log(screen_log, log_file)

    # Step 2 confirm/make the decision based on BLAST of adjacent genes


    # Step 3 annotate each of CRISPR1 and 2

    # Output Fasta, Gff and Spacer Arrays
    

    # remove tmp dir
    shutil.rmtree(tmpdir)

    return


def rm_overlaps(blast_out):
    """
    Remove overlapped blast hits
    Use constant variable MAX_OVERLAPPING

    Args:
        blast_out (Dataframe): columns: ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'slen', 'qseq']

    Returns:
        blast_out: blast_out that removed overlaps
    """

    blast_out.sort_values(by=['qseqid', 'qstart'],
                          inplace=True,
                          ignore_index=True)

    rows_to_return = []

    for index, row in blast_out.iterrows():
        # Assume all the query records are forward strand, in other words, always 'qstart' < 'qend'
        if rows_to_return == []:
            rows_to_return.append(row)
            continue

        last_row = rows_to_return[-1]
        # if the 2 blast record not overlapped
        if (row['qseqid'] != last_row['qseqid']) or (
                last_row['qend'] - row['qstart'] < MAX_OVERLAPPING):
            rows_to_return.append(row)
            continue

        # else if they are overlapped,
        # 1. If one of them have 100%, keep it
        # 2. If more than 1 matches are 100% identity, or none is
        #     a. keep "DR" caused it's the most common DR sequence
        #     b. keep the higher bitscore
        ## The name of DR actually doesn't matter, just to make sure find a proper match for identify CRISPR allele.
        if (last_row['pident'] == 100) and (row['pident'] != 100):
            continue
        elif (last_row['pident'] != 100) and (row['pident'] == 100):
            rows_to_return[-1] = row
            continue
        elif last_row['sseqid'] == "DR":
            continue
        elif row['pident'] == "DR":
            rows_to_return[-1] = row
            continue
        elif row['bitscore'] > last_row['bitscore']:
            rows_to_return[-1] = row
            continue
        # should cover all possibilities
        else:
            continue
            
    blast_out = pd.DataFrame.from_records(rows_to_return)

    return blast_out


def find_crispr_alleles(blast_out):
    
    # Crispr_alleles is a list of DataFrame then output as GFF
    crispr_alleles = []

    # Rows return to one CRISPR allele
    crispr_allele_rows = []
    for index, row in blast_out.iterrows():
        if len(crispr_allele_rows) == 0:
            crispr_allele_rows.append(row)
            continue

        last_row = crispr_allele_rows[-1]
        # if current row is close to the previous one
        # Loop could finish within this condition without saving the crispr_allele_rows to returning list
        if (row['qseqid'] == last_row['qseqid']
            ) and (row['qstart'] - last_row['qend']) <= MAX_SPACER_LENGTH:
            crispr_allele_rows.append(row)
            continue
        # else, if current row is far away to the previous one
        # If found less than 2 DR, drop the rows, start again with current row
        elif len(crispr_allele_rows) < MIN_REPEATS:
            crispr_allele_rows = [row]
            continue
        # Else, put the allele to the returning list, and start again with current row
        else:
            crispr_allele = pd.DataFrame.from_records(crispr_allele_rows)
            crispr_alleles.append(crispr_allele)
            crispr_allele_rows = [row]
            continue

    # When loop finishes, if crispr_allele_rows have more than 2 DR
    if len(crispr_allele_rows) >= MIN_REPEATS:
        crispr_allele = pd.DataFrame.from_records(crispr_allele_rows)
        crispr_alleles.append(crispr_allele)

    return crispr_alleles


def blast2gff(crispr_allele):
    # decide strand
    crispr_allele['strand'] = crispr_allele.apply(
        lambda row: '+' if row['qstart'] < row['qend'] else '-', axis=1)
    strand = max(list(crispr_allele['strand']),
                 key=list(crispr_allele['strand']).count)

    # first line - the region of the allele
    gff_records_to_return = []
    crispr = [
        str(crispr_allele.loc[0, 'qseqid']), 'blast', 'CRISPR',
        str(crispr_allele.loc[0, 'qstart']),
        str(crispr_allele.loc[len(crispr_allele) - 1, 'qend']),
        str(len(crispr_allele)), strand, '.', "name=CRISPR"
    ]
    gff_records_to_return.append(crispr)

    # following lines - each DR
    for index, row in crispr_allele.iterrows():
        gff_record = [
            str(row['qseqid']), 'blast', 'DR',
            str(row['qstart']),
            str(row['qend']), '{:.2f}'.format(row['pident']), row['strand'],
            '.',
            "name={id};length={len}".format(id=row['sseqid'],
                                            len=row['qend']-row['qstart']+1)
        ]
        gff_records_to_return.append(gff_record)

    return gff_records_to_return
