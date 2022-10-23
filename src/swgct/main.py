#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Parse the inpot arguments
date: 2022-10-03
'''
__author__ = 'Yan Li'
__email__ = 'yan.li2@outlook.com'

import sys, os
import argparse
from swgct import __version__
from swgct.swgct_run import run
from swgct.swgct_summary import summary

CURRENT_PATH = os.path.dirname(__file__)
DB_PATH = CURRENT_PATH + '/db/'


# Parse the user input arguments
def get_arguments():

    parser = argparse.ArgumentParser(
        description=
        "Salmonella Whole genome CRISPR Typing pipeline: Detect CRISPR1 and CRISPR2 alleles of Salmonella genomes",
        epilog = "Use \"swgct [command] -h\" for more information about a command." 
    )
    parser.add_argument('-v',
                        '--version',
                        action='version',
                        version='%(prog)s ' + __version__)

    # Add sub commands
    subparsers = parser.add_subparsers(title='Available commands',
                                       required=True)

    ## 1. parsers for sub command "swgct run"

    parser_run = subparsers.add_parser(
        'run',
        help="Run Salmonella Whole genome CRISPR Typing on each genome.")

    # Input & output options
    io_group = parser_run.add_argument_group("Input/Output")
    io_group.add_argument("input",
                          metavar="fasta",
                          help="Input FASTA file (required)")
    io_group.add_argument(
        '-o',
        '--outdir',
        required=False,
        metavar='outdir',
        help=
        'Output directory. By default it\'s a folder with the same name as input'
    )
    io_group.add_argument(
        "-p",
        "--prefix",
        required=False,
        metavar="prefixed",
        help=
        "Prefixed output file base names. By default it\'s the same with '--outdir'"
    )

    # Other options
    other_group = parser_run.add_argument_group('Other')
    other_group.add_argument(
        '-d',
        '--database_path',
        metavar='db_path',
        default=DB_PATH,
        help='BLAST databases of Salmonella CRISPR DR and spacers, default={}'.
        format(DB_PATH))
    other_group.add_argument(
        '--algorithm',
        choices=['blast', 'minced'],
        default='blast',
        help=
        'The algorithm used for searching CRISPR alleles. They may generate different results. Default="blast". "minced" is an option when the "blast" does not work.'
    )
    other_group.add_argument(
        '-k',
        '--keep',
        action='store_true',
        help='Keep the tmp files, which will be cleaned by default')

    parser_run.set_defaults(func=run)

    # 2. Subcommand "swgct summary"

    parser_summary = subparsers.add_parser(
        'summary',
        help="Summarize the result from SWGCT",
        description=
        "Summarize the result from SWGCT, expected to be run under a directory containing all the results from 'swgct run'."
    )

    parser_summary.add_argument(
        '-d',
        '--dir',
        default='.',
        metavar='dir',
        help='Directory to work with, default is current directory.'
    )

    parser_summary.set_defaults(func=summary)

    # If no arguments were used, print the entire help
    if (len(sys.argv) == 1):
        parser.print_help(file=sys.stderr)
        sys.exit(1)

    # Go
    args = parser.parse_args()

    return args


def main():

    args = get_arguments()
    args.func(args)

    return


if __name__ == '__main__':
    main()
