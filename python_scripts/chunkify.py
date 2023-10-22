import sys
from Bio import SeqIO

filename = sys.argv[1]

mode = sys.argv[2]

if mode == "-pre":
	with open(filename) as seqfilehandle:
		with open(filename+".chunked.fas", "w") as outf:
			seqs = SeqIO.parse(seqfilehandle, "fasta")
			for seq in seqs:
				lnseq = len(seq.seq)
				for pos in range(0, lnseq, 1000):
					if lnseq > pos+1000:
						endB = pos+1000
					else:
						endB = lnseq
					print >> outf, ">"+seq.id+"_s"+str(pos)+"_e"+str(endB)
					print >> outf, seq.seq[pos:endB]
elif mode == "-post":
	with open(filename) as blastfilehandle:
		with open(filename+"_post", "w") as outf:
			for line in blastfilehandle:
				splitline = line.strip().split("\t")
				splitname = splitline[0].split("_")
				startB = int(splitname[-2][1:])
				newcontigname = "_".join(splitname[:-2])
				newquerystart = str(int(splitline[6])+startB)
				newqueryend = str(int(splitline[7])+startB)
				oldlinechunk1 = "\t".join(splitline[1:6])
				oldlinechunk2 = "\t".join(splitline[8:])
				print >> outf, newcontigname+"\t"+oldlinechunk1+"\t"+newquerystart+"\t"+newqueryend+"\t"+oldlinechunk2

else:
	print "error in parameter"
	sys.exit()