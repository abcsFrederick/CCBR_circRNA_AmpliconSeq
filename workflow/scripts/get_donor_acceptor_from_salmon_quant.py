import HTSeq
import sys
import argparse
import os

def read_bed_file(filename):
	bedfile=open(filename,'r')
	primers=dict()
	for f in bedfile.readlines():
		f=f.strip().split("\t")
		primer=f[3]
		primers[primer]=dict()
		primers[primer]["chrom"]=f[0]
		primers[primer]["start"]=int(f[1])
		primers[primer]["end"]=int(f[2])
		primers[primer]["strand"]=f[5]
	return primers

parser = argparse.ArgumentParser(description='Create fasta file from divergent primers bed file')
parser.add_argument('--bed', dest='primerbed', type=str, required=True,
                    help='Divergent primers bed file with min. 4 columns (chr,start,end,name) name is expected to have _AS/_SS suffix')
parser.add_argument('--reffa', dest='reffa', type=str, required=True,
                    help='reference fasta')
parser.add_argument('--outfa', dest='outfa', type=str, required=True,
                    help='output fasta')
parser.add_argument('--scanlength', dest='scanlen', type=int, required=False, default=200,
                    help='scan length...even positive number(default 200)')
parser.add_argument('--flankmax', dest='flankmax', type=int, required=False, default=30,
                    help='flankmax (default 30)')
args = parser.parse_args()
# sequences = dict( (s.name, s) for s in HTSeq.FastaReader(args.reffa) )
sequences = dict( (s[1], s[0]) for s in HTSeq.FastaReader(args.reffa, raw_iterator=True) )
# sequences = dict( (s[1], s[0]) for s in HTSeq.FastaReader(args.reffa) )
primers = read_bed_file(args.primerbed)