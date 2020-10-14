#!/usr/bin/env python

import argparse
import os
from Bio import SeqIO
import sys
import gzip

def main(args):
    '''Main program'''
    
    # open handles
    if args.i == 'stdin':
        fh_in = sys.stdin
    elif args.i.endswith('.gz'):
        fh_in = gzip.open(args.i, 'rt')
    else:
        fh_in = open(args.i, 'r')
    
    if args.o == 'stdout':
        fh_out = sys.stdout
    else:
        fh_out = open(args.o, 'w')
    
    # parse
    for record in SeqIO.parse(fh_in, 'fasta'):
        if len(record) >= args.min:      
            c = SeqIO.write(record, fh_out, 'fasta')
    
    fh_in.close()
    fh_out.close()
    

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(prog='fasta_filterlength.py',
                                 formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=25, width=140),
                                 description='''Filter fasta for length''', 
                                 usage='%(prog)s [options]')
    
    # add the arguments
    parser.add_argument('--i', help='input fasta [stdin]', default='stdin')
    parser.add_argument('--min', help='min length [0]', default=0, type=int)
    parser.add_argument('--o', help='output file [stdout]', default='stdout')
    
    args = parser.parse_args()
    #args = parser.parse_args('--i contigs.fna.gz --min 2000 --o test.out'.split())
    
    main(args)

