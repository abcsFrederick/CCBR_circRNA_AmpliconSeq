import sys
primerstsv=open(sys.argv[1]).readlines()
primerstsv.pop(0)
primers=dict()
for i in primerstsv:
	i=i.strip().split("\t")
	circRNAnamesplit=i[0].split("_")
	circRNAname="_".join(circRNAnamesplit[:-1])
	if not circRNAname in primers:
		primers[circRNAname]=dict()
		primers[circRNAname]["coordinates"]=list()
	ForR=circRNAnamesplit[-1][0]
	primers[circRNAname]["chrom"]=i[1]
	strand=i[4]
	primers[circRNAname]["coordinates"].append(int(i[2]))
	primers[circRNAname]["coordinates"].append(int(i[3]))
	primers[circRNAname]["coordinates"].sort()
	if ForR == "F":
		primers[circRNAname]["strand"]=strand
	if ForR == "F" and strand == "+":
		primers[circRNAname]["ASorSS"]="SS"
	if ForR == "F" and strand == "-":
		primers[circRNAname]["ASorSS"]="AS"
for primer in primers:
	sname=primer+"_"+primers[primer]["ASorSS"]
	#sname="_".join([primer,primers[primer]["ASorSS"])
	#print(sname)
	chrom=primers[primer]["chrom"]
	start=str(primers[primer]["coordinates"][0])
	end=str(primers[primer]["coordinates"][-1])
	strand=primers[primer]["strand"]
	print("\t".join([chrom,start,end,sname,".",strand]))
