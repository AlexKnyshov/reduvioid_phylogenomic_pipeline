from Bio import SeqIO
import sys

listfilename = sys.argv[1]
foldername = sys.argv[2]
d = {}
with open(listfilename) as listf:
	for line in listf:
		splitline = line.strip().split(",")
		transname = splitline[2]
		locname = splitline[0]
		if transname in d:
			if locname not in d[transname]:
				d[transname].add(locname)
		else:
			d[transname] = set()
			d[transname].add(locname)

for k,v in d.items():
	# print k, "a".join(v)
	with open(foldername+"/"+k+".fas") as infh:
		seqs = SeqIO.read(infh, "fasta")
		locname = "a".join(v)
		with open(locname+".fas","w") as outh:
			print >> outh, ">"+locname
			print >> outh, seqs.seq