#!/usr/bin/env python

import vamb
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--i', help='input file')
parser.add_argument('--o', help='output file')
args = parser.parse_args()

input_f = args.i
output_f = args.o
#i_filename = os.path.split(input)[1]
#o_filename = os.path.split(output)[1]

with open(input_f, "r") as file:
    depths = vamb.vambtools.load_jgi(file)

vamb.vambtools.write_npz(output_f, depths)
