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

def complement(c):
	if c=="A":
		return "T"
	elif c=="C":
		return "G"
	elif c=="G":
		return "C"
	elif c=="T":
		return "A"
	elif c=="N":
		return "N"
	else:
		print("Unknown char %s. expecting A,C,G,T or N"%(c))
		exit()

def revcom(seq):
	rc=seq[::-1].upper()
	rc="".join(list(map(lambda x:complement(x),rc)))
	return rc

def convertnt(c):
	if c=="A":
		return "1"
	elif c=="C":
		return "4"
	elif c=="G":
		return "3"
	elif c=="T":
		return "2"
	elif c=="N":
		return "0"
	else:
		return "-1"

parser = argparse.ArgumentParser(description='Create fasta file from divergent primers bed file')
parser.add_argument('--bed', dest='primerbed', type=str, required=True,
                    help='Divergent primers bed file with min. 4 columns (chr,start,end,name), name is expected to have _AS/_SS suffix')
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
# print(primers)
# for primer in primers:
# 	print(primers[primer]["chrom"])
outfastafile = open( args.outfa, "w" )
offset = "N"*300
for primer in primers:
	ASorSS=primer.split("##")[-1]
	seq=sequences[primers[primer]["chrom"]]
	for i in range(primers[primer]["start"]-args.scanlen,primers[primer]["start"],1):
		for j in range(primers[primer]["end"],primers[primer]["end"]+args.scanlen,1):
			k=int((j-i)/2)
			if k >= args.flankmax:
				bsjflanka=seq[ j - args.flankmax : j ]
				bsjflankb=seq[ i : i + args.flankmax ]
			else:
				bsjflanka=seq[ j-k : j ]
				bsjflankb=seq[ i : i+k ]
			if ASorSS == "SS":
				five_da=bsjflanka[-2:]+"-"+seq[j:j+2]
				three_da=seq[i-2:i]+"-"+bsjflankb[:2]
				da=seq[j:j+2]+"-"+seq[i-2:i]
			elif ASorSS == "AS":
				five_da=revcom(seq[i-2:i])+"-"+revcom(bsjflankb[:2])
				three_da=revcom(bsjflanka[-2:])+"-"+revcom(seq[j:j+2])
				da=revcom(seq[i-2:i])+"-"+revcom(seq[j:j+2])
			else:
				print("Primer %s does not have AS or SS suffix"%(primer))
				exit()
			bsjflank=bsjflanka+bsjflankb
			bsjflank=bsjflank.upper()
			bsjflank_nt="".join(list(map(lambda x:convertnt(x),bsjflank)))				
			sname = "##".join([primer,primers[primer]["chrom"],str(i),str(j),five_da,three_da,da,bsjflank_nt])
			# print(sname)
			# print(bsjflank)
			myseq = HTSeq.Sequence( bytes(offset+bsjflank+offset, 'utf-8'), sname)
			myseq.write_to_fasta_file(outfastafile)
outfastafile.close()

