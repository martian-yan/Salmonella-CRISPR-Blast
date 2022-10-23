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
import numpy as np

from Bio import SeqIO
from Bio.Seq import Seq

from swgct.misc import cmd, current_time, write_log

MIN_REPEATS = 2  # define CRISPR allele with at least 2 DR in a row
MAX_SPACER_LENGTH = 150  # allow to skip one or two DR
MIN_NEW_SPACER_LENGTH = 15 # find new spacers
MAX_OVERLAPPING = 5  # dont allow overlaps of DRs
FLANKING_LENGTH = 1000
CRISPR_SORTED = False

## Main modules

def blast_crispr_dr(fasta_file, dr_db, out_dir, out_name, tmpdir):
    """
    Read the DR Blast result, write CRISPR alleles and direct repeats to a GFF file.

    Args:
        fasta_file (str): input fasta file
        dr_db (str): BLAST database of direct repeats
        out_dir (str): output file folder
        out_name (str): output file name
        tmpdir (str): temporary directory

    Output files:
        log_file = out_file+'.log',
        sorted_blast_out = out_file+'.blast_DR.sorted.tsv',
        gff_file = out_file+'.DR_blast.gff'
    """

    ## Output
        ## log
    tmp_file = "{tmpdir}/{out_name}".format(tmpdir=tmpdir, out_name=out_name)
    out_file = "{outdir}/{out_name}".format(outdir=out_dir, out_name=out_name)
    log_file = out_file + '.log'

    # If Outdir not exist, create one
    outdir_exist = os.path.exists(out_dir)
    if not outdir_exist:
        os.mkdir(out_dir)

    # Start log
    screen_log = "{time}: START looking for CRISPR in {fasta}".format(
        time=current_time(), fasta=fasta_file)
    write_log(screen_log, log_file)

    ## Step 1, BLAST against the DR database
    blast_out_file=tmp_file + '.blast_DR.tsv'
    blast_cmd = "blastn -db {dr_db} -query {fasta_file} -out {blast_out} -outfmt '6 std slen qseq' -task blastn -evalue 0.001 -perc_identity 90 -max_target_seqs 10000 -culling_limit 5; awk '$4==$13' {blast_out} > {blast_out}.tmp && mv {blast_out}.tmp {blast_out}".format(
        fasta_file=fasta_file,
        dr_db=dr_db,
        blast_out=blast_out_file,
    )
    cmd(blast_cmd, "Blast DR sequences", log_file)

    ## Step 2, get 2 CRISPR allele positions
    # Print log
    screen_log = "Looking for CRISPR alleles"
    write_log(screen_log, log_file)

    blast_out = pd.read_csv(blast_out_file,
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
        gff_records += repeats_to_gff(crispr_allele)

    gff_file = tmp_file + '.DR_blast.gff'
    with open(gff_file, 'w') as fo:
        fo.write('# GFF file created from blast result {}\n'.format(
            blast_out_file))
        fo.write('# BLAST cmd: {}\n'.format(blast_cmd))
        fo.write('# The float in column 6 is identity%\n')
        fo.write('#\n')
        for record in gff_records:
            fo.write('{}\n'.format('\t'.join(record)))

    return


def sort_crispr_alleles(fasta_file, out_dir, out_name, tmpdir):
    """
    1. Decide CRISPR1 and CRISPR2 sites
    2. Output GFF and a simple table with CRISPR arrays
        
    Args:
        fasta_file (str): input fasta file
        out_dir (str): output file folder
        out_name (str): output file name
        tmpdir (str): temporary directory
    """

    # input
    out_file = "{outdir}/{out_name}".format(outdir=out_dir, out_name=out_name)
    tmp_file = "{tmpdir}/{out_name}".format(tmpdir=tmpdir, out_name=out_name)
    dr_gff = tmp_file + '.DR_blast.gff'
    log_file = out_file + '.log'

    crispr_alleles = read_gff(dr_gff)
    crispr_alleles = crispr_alleles[crispr_alleles['feature'] == 'CRISPR']

    strand = "Unclear"
    sort_ascending = True
    if all(strand == '+' for strand in crispr_alleles['strand']):
        strand = "Forward"
    elif all(strand == '-' for strand in crispr_alleles['strand']):
        strand = "Backward"
        sort_ascending = False

    crispr_alleles.sort_values(by=['qseqid', 'start'],
                               inplace=True,
                               ascending=[True, sort_ascending],
                               ignore_index=True)

    # TODO: Possible upgrade - train a ML model to dicide whether it's CRISPR1 or 2
    # Step 1 decide CRISPR1 and 2 through the order if they on the SAME contig
    if (len(crispr_alleles) == 2) and (strand != "Unclear") and \
        all(qseqid == crispr_alleles.loc[0, 'qseqid'] for qseqid in crispr_alleles['qseqid']):
        crispr_alleles.loc[0, 'feature'] = 'CRISPR1'
        crispr_alleles.loc[0, 'attributes'] = 'CRISPR1'
        crispr_alleles.loc[1, 'feature'] = 'CRISPR2'
        crispr_alleles.loc[1, 'attributes'] = 'CRISPR2'
        CRISPR_SORTED = True

        print(crispr_alleles)
        screen_log = "CRISPR alleles sorted according to their position:\n\
                    CRISPR1 is at: {0}..{1}\n\
                    CRISPR2 is at: {2}..{3}".format(crispr_alleles.loc[0, 'start'],
                                                    crispr_alleles.loc[0, 'end'],
                                                    crispr_alleles.loc[1, 'start'],
                                                    crispr_alleles.loc[1, 'end'])
    else:
        screen_log = "CRISPR alleles CAN NOT be sorted\n"
 
    write_log(screen_log, log_file)

    # Step 2 confirm/make the decision based on BLAST of adjacent genes

    # Output Fasta, Gff and Spacer Arrays
    crispr_alleles_flanking = crispr_alleles.copy()
    crispr_alleles_flanking['start'] = crispr_alleles_flanking['start'] - FLANKING_LENGTH
    crispr_alleles_flanking['end'] = crispr_alleles_flanking['end'] + FLANKING_LENGTH
    crispr_alleles_flanking_gff = tmp_file+".crisprs.flanking.gff"

    df_to_gff(crispr_alleles_flanking, crispr_alleles_flanking_gff)

    crispr_gff = tmp_file+".crispr.gff"
    df_to_gff(crispr_alleles, crispr_gff)

    # getfasta
    crispr_alleles_fasta = out_file + ".CRISPR.fasta"
    get_fasta = "bedtools getfasta -fi {fasta_file} -fo {crispr_alleles_fasta} -bed {crispr_alleles_gff} -nameOnly ; sed -i 's/>/>{out_name}_/g' {crispr_alleles_fasta}".format(
        fasta_file=fasta_file,
        crispr_alleles_fasta=crispr_alleles_fasta,
        crispr_alleles_gff=crispr_alleles_flanking_gff,
        out_name=out_name)
    cmd(get_fasta, "Get CRISPR fasta with {}bp length flanking sequences".format(FLANKING_LENGTH), log_file)

    return


def annotate_crispr(fasta_file, spacer_db, out_dir, out_name, tmpdir):
    """
    Annotate the spacers

    Args:
        fasta_file (str): input fasta file
        spacer_db (str): BLAST database
        out_dir (str): output file folder
        out_name (str): output file name
        log_file (str): log file
        tmpdir (str): temporary directory
    """
    tmp_file = "{tmpdir}/{out_name}".format(tmpdir=tmpdir, out_name=out_name)
    out_file = "{outdir}/{out_name}".format(outdir=out_dir, out_name=out_name)
    log_file = out_file + '.log'

    ## Step 1, BLAST against the spacer database
    blast_out_file = tmp_file + '.blast_spacer.tsv'

    blast_cmd = "blastn -db {spacer_db} -query {fasta_file} -out {blast_out} -outfmt '6 std slen qseq' -task blastn -evalue 0.001 -perc_identity 90 -max_target_seqs 10000 -culling_limit 5; awk '$4==$13' {blast_out} > {blast_out}.tmp && mv {blast_out}.tmp {blast_out}".format(
        fasta_file=fasta_file,
        spacer_db=spacer_db,
        blast_out=blast_out_file
    )
    cmd(blast_cmd, "Blast spacers", log_file)

    blast_out = pd.read_csv(blast_out_file,
                        sep='\t',
                        names=[
                            'qseqid', 'sseqid', 'pident', 'length',
                            'mismatch', 'gapopen', 'qstart', 'qend',
                            'sstart', 'send', 'evalue', 'bitscore', 'slen',
                            'qseq'
                        ])
    # change datatype of qseqid to avoid bug caused by comparing numeric contig names
    blast_out = blast_out.astype({'qseqid': 'str'})

    ## Sort spacers by each CRISPR allele
    crispr_gff = tmp_file+".crispr.gff"
    crispr_alleles = read_gff(crispr_gff)
    new_spacers = []
    new_spacer_var = []
    gff_to_write = []
    spacer_arrays = []
    for index, crispr_allele in crispr_alleles.iterrows():
        new_spacers_itter = []
        spacers = blast_out[(blast_out['qseqid']==crispr_allele['qseqid']) & (blast_out['qstart']>=crispr_allele['start']) & (blast_out['qend']<=crispr_allele['end'])]
        spacers = rm_overlaps(spacers) 

        # Check new spacer var: identity < 100
        for index, row in spacers.iterrows():
            if row['pident'] < 100:
                row['sseqid'] = row['sseqid'] + '_var'
                spacers.at[index, 'sseqid'] = row['sseqid']
                # save new var to olsutput
                new_spacer_var.append(row)
        
        # Check gaps which could be new spacers
        for i in range(1, len(spacers)):
            if spacers.loc[i, 'qstart'] - spacers.loc[i-1, 'qend'] > MIN_NEW_SPACER_LENGTH:
                # save new spacers
                new_spacer = get_new_spacer(fasta_file, spacers.loc[i, 'qseqid'], crispr_allele['strand'], spacers.loc[i-1, 'qend'], spacers.loc[i, 'qstart'])
                new_spacers_itter.append(new_spacer)
                new_spacers.append(new_spacer)
        
        if new_spacers_itter:
            new_spacers_df = pd.DataFrame.from_records(new_spacers_itter)
            # Put new spacers in the output GFF file
            spacers = pd.concat([spacers, new_spacers_df], ignore_index=True)
            spacers.sort_values(by=['qseqid', 'qstart'],
                          inplace=True,
                          ignore_index=True)

        # prepare gff files        
        gff_records = spacer_to_gff(crispr_allele, spacers)
        gff_to_write.extend(gff_records)
        
        # prepare output spacer arrays
        spacer_array = spacers['sseqid'].tolist()
        spacer_array = [x for x in spacer_array if "DR" not in x]
        if crispr_allele['strand'] == '-':
            spacer_array.reverse()
        spacer_array = "-".join(spacer_array)
        spacer_arrays.append("{0}\t{1}\n".format(crispr_allele['feature'], spacer_array))


    ## Write files
    out_gff = out_file + '.CRISPR.gff'
    with open(out_gff, 'w') as fo:
        fo.write('# GFF file created from blast result {}\n'.format(
            blast_out_file))
        fo.write('# BLAST cmd: {}\n'.format(blast_cmd))
        fo.write('# The float in column 6 is identity%\n')
        fo.write('#\n')
        for record in gff_to_write:
            fo.write('{}\n'.format('\t'.join(record)))

    out_summary = out_file + '.CRISPR.tsv'
    with open(out_summary, 'w') as fo:
        fo.writelines(spacer_arrays)

    if new_spacers:
        # output new spacers
        new_spacers_df = pd.DataFrame.from_records(new_spacers)
        new_spacer_gff = blast_to_gff(new_spacers_df)
        new_spacer_gff_file = out_file + '.new_spacer.gff'
        df_to_gff(new_spacer_gff, new_spacer_gff_file)
        new_spacer_fasta_file = out_file + '.new_spacer.fasta'
        screen_log = "\nFound {n} new spacers".format(n=len(new_spacers))
        write_log(screen_log, log_file)

        get_fasta = "bedtools getfasta -fi {fasta_file} -fo {new_spacer_fasta} -bed {new_spacer_gff} -nameOnly".format(
            fasta_file=fasta_file,
            new_spacer_fasta=new_spacer_fasta_file,
            new_spacer_gff=new_spacer_gff_file,)
        cmd(get_fasta, "Get fasta of new spacers", log_file)  

    if new_spacer_var:
        # out put new_spacer_var
        new_spacer_var_df = pd.DataFrame.from_records(new_spacer_var)
        new_spacer_var_gff = blast_to_gff(new_spacer_var_df)
        new_spacer_var_gff_file = out_file + '.new_spacer_var.gff'
        df_to_gff(new_spacer_var_gff, new_spacer_var_gff_file)
        new_spacer_var_fasta_file = out_file + '.new_spacer_var.fasta'
        screen_log = "\nFound {n} new spacers var".format(n=len(new_spacer_var))
        write_log(screen_log, log_file)

        get_fasta = "bedtools getfasta -fi {fasta_file} -fo {new_var_fasta} -bed {new_var_gff} -nameOnly".format(
            fasta_file=fasta_file,
            new_var_fasta=new_spacer_var_fasta_file,
            new_var_gff=new_spacer_var_gff_file,)
        cmd(get_fasta, "Get fasta of new spacers vars", log_file)

    return


## Supp modules

def rm_overlaps(blast_out):
    """
    Remove overlapped blast hits
    Use constant variable MAX_OVERLAPPING

    Args:
        blast_out (Dataframe): columns: ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'slen', 'qseq']

    Returns:
        blast_out: blast_out that removed overlaps
    """

    blast_out = blast_out.sort_values(by=['qseqid', 'qstart'],
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
        elif row['sseqid'] == "DR":
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

def repeats_to_gff(crispr_allele):
    """Turn DR blast results into gff output

    Args:
        crispr_allele (dataframe): a dataframe of DR blast results

    Returns:
        list: a list of gff records
    """    
    # decide strand
    crispr_allele['strand'] = crispr_allele.apply(
        lambda row: '+' if row['sstart'] < row['send'] else '-', axis=1)
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
            '.', "name={id};seq={seq}length={len}".format(id=row['sseqid'],
                                                 len=row['qend']-row['qstart']+1,
                                                 seq=row['qseq'])
        ]
        gff_records_to_return.append(gff_record)

    return gff_records_to_return

def spacer_to_gff(crispr_allele, spacers):
    """make gff records from the blast results

    Args:
        crispr_allele (pd.series): a row of gff file shows whole CRISPR allele
        spacers (DataFrame): a dataframe of spacers BLAST result

    Returns:
        list: a list of gff records
    """    

    gff_records_to_return = []
    crispr_allele_str = [str(int(x)) if type(x)!=str else x for x in crispr_allele ]
    gff_records_to_return.append(crispr_allele_str)

    # following lines - each spacer
    crispr = crispr_allele['feature']
    for index, row in spacers.iterrows():
        if "DR" in row['sseqid']:
            feature = 'DR'
        else:
            feature = 'spacer'

        if row['sstart'] < row['send']:
            strand = '+'
        else:
            strand = '-'

        gff_record = [
            str(row['qseqid']), 'blast', feature,
            str(int(row['qstart'])),
            str(int(row['qend'])), '{:.2f}'.format(row['pident']), strand,
            '.', "name={id};parent={crispr};seq={seq};length={len}".format(id=row['sseqid'],
                                                 crispr=crispr,
                                                 len=row['qend']-row['qstart']+1,
                                                 seq=row['qseq'])
        ]
        gff_records_to_return.append(gff_record)

    return gff_records_to_return


def blast_to_gff(blast_out):

    gff_records_to_return = []

    for index, row in blast_out.iterrows():
        
        seq = row['qseq']

        if row['sstart'] < row['send']:
            strand = '+'

        else:
            strand = '-'
            seq = Seq(seq)
            seq = seq.reverse_complement()
            seq = str(seq).upper()

        gff_record = {
            'qseqid': str(row['qseqid']),
            'source': 'blast',
            'feature': row['sseqid'],
            'start': row['qstart'],
            'end': row['qend'], 
            'score': '{:.2f}'.format(row['pident']), 
            'strand': strand,
            'frame': '.', 
            'attribute': "name={id};seq={seq};length={len}".format(id=row['sseqid'],
                                                    len=row['qend']-row['qstart']+1,
                                                    seq=seq)
        }
        gff_records_to_return.append(gff_record)

    gff_records_to_return = pd.DataFrame.from_records(gff_records_to_return)

    return gff_records_to_return

def df_to_gff(df, gff_file):
    """Output a df as a gff file

    Args:
        df (DataFrame): dataframe with gff style column names 
        gff_file (str): path to gff file
    """
    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)     
    df.to_csv(gff_file, sep='\t', index=False, header=False, )
    return

def read_gff(gff_file):
    """Read Gff files as a dataframe

    Args:
        gff_file (str): path to gff file
    """    

    df = pd.read_csv(gff_file,
                    sep='\t',
                    header=None,
                    comment='#',
                    names=[
                        'qseqid', 'source', 'feature', 'start', 'end',
                        'score', 'strand', 'frame', 'attributes'
                    ],
                    dtype={
                        'qseqid': 'str',
                        'source': 'str',
                        'feature': 'str',
                        'start': np.float64,
                        'end': np.float64,
                        'score': 'str',
                        'strand': 'str',
                        'frame': 'str',
                        'attributes': 'str'
                    })

    return df

def get_new_spacer(fasta, contig, strand, start, end):

    length = end-start+1
    sstart = 1
    send = length

    fastas = list(SeqIO.parse(fasta, 'fasta'))
    seq = [seq_record for seq_record in fastas if seq_record.id == contig]
    seq = seq[0].seq[start:end-1]
    if strand == '-':
        seq = seq.reverse_complement()
        sstart = length
        send = 1
    seq = str(seq).upper()


    series_to_return={
        'qseqid': contig, 'sseqid': "new_spacer", 'pident': 100, 'length': length, 'mismatch': 0, 'gapopen': 0, 'qstart': start+1, 'qend': end-1, 'sstart': sstart, 'send': send, 'evalue': 0, 'bitscore': 100, 'slen': length, 'qseq': seq
    }
    
    return series_to_return