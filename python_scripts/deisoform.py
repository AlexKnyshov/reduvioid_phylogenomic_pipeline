from Bio import SeqIO
import sys
from itertools import combinations

infilename = sys.argv[1]
seqs = SeqIO.parse(infilename, "fasta")
d = {}
for seq in seqs:
	taxname = seq.id.split("|")[0]
	if taxname in d:
		d[taxname].append(seq.seq)
	else:
		d[taxname] = [seq.seq]

outfilename = infilename.split("/")[-1]
d2 = {}
d3 = {}
d4 = {}
d5 = {}
with open(outfilename, "w") as outhandle:
	for taxon in d:
		if len(d[taxon]) == 1:
			print >> outhandle, ">"+taxon
			print >> outhandle, d[taxon][0]
			d2[taxon] = 0
			d3[taxon] = 0
			d4[taxon] = 1
			d5[taxon] = 1
		else:
			isonumber = 0
			badseqs = set()
			for comb in combinations(range(len(d[taxon])), 2):
				seqlen = len(d[taxon][0])
				seq1 = d[taxon][comb[0]]
				seq2 = d[taxon][comb[1]]
				diffs = 0.0
				for pos in range(len(d[taxon][0])):
					if seq1[pos] != "-" and seq2[pos] != "-":
						if seq1[pos] != seq2[pos]:
							diffs += 1
					else:
						seqlen -= 1
				if seqlen > 10:
					if diffs / seqlen <= 0.01:
						if comb[0] not in badseqs and comb[1] not in badseqs:
							isonumber += 1
							if len(str(comb[0]).replace("-","")) < len(str(comb[1]).replace("-","")):
								badseqs.add(comb[0])
							else:
								badseqs.add(comb[1])
			d3[taxon] = isonumber

			chunknumber = 0
			if len(d[taxon]) - isonumber > 1:
				for seqnum in range(len(d[taxon])):
					if seqnum not in badseqs:
						# print "seq1", seqnum
						ovlp = 0
						chunk = True
						for seqnum2 in range(len(d[taxon])):
							if seqnum2 not in badseqs:
								if seqnum != seqnum2:
									# print "seq2", seqnum2
									for pos in range(len(d[taxon][0])):
										if d[taxon][seqnum][pos] != "-" and d[taxon][seqnum2][pos] != "-":
											ovlp += 1
										if ovlp > 10:
											# print "overlap"
											break
							if ovlp > 10:
								chunk = False
								break
						if chunk:
							chunknumber += 1
			d2[taxon] = chunknumber
			finalnumber=0
			for seqnum in range(len(d[taxon])):
				if seqnum not in badseqs:
					if len(d[taxon]) - len(badseqs) > 1:
						print >> outhandle, ">"+taxon+"_"+str(seqnum)
					else:
						print >> outhandle, ">"+taxon
					print >> outhandle, d[taxon][seqnum]
					finalnumber += 1
			d4[taxon] = finalnumber
			if finalnumber - chunknumber == 0:
				copynumber = 1
			else:
				copynumber = finalnumber - chunknumber
			d5[taxon] = copynumber
with open(outfilename+".txt", "w") as outhandle2:
	for tx in d:
		print >> outhandle2, tx, len(d[tx]), d2[tx], d3[tx], d4[tx], d5[tx]
#d2[tx] - extra chunk number
#d3[tx] - number of isoforms
#d4[tx] - number of seqs retained in the files
#d5[tx] - copy number