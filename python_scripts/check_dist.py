import sys
from Bio import SeqIO

infilename = sys.argv[1]

with open(infilename) as infilehandle:
	seqs = SeqIO.parse(infilehandle, "fasta")
	seqlist = []
	for seq in seqs:
		seqlist.append(str(seq.seq))

def pairwise_dist(seq1, seq2):
	seqlen = len(seq1)
	# print(seqlen)
	diffs = 0.0
	for pos in range(seqlen):
		seqset = set()
		if seq1[pos] != "-" and seq1[pos] != "N":
			seqset.add(seq1[pos])
		if seq2[pos] != "-" and seq2[pos] != "N":
			seqset.add(seq2[pos])
		# print seqset
		if len(seqset) > 1:
			diffs += 1
	if diffs/seqlen > 0.05:
		return False
	else:
		return True

def fetch_fasta(cl,lcn):
	if len(cl) == 1:
		maxlen = 0
		longestSeq = ""
		for seq in cl[0]:
			if len(seq.replace("-","").replace("N","")) > maxlen:
				longestSeq = seq
		with open(lcn+".fas", "w") as outh:
			print >> outh, ">"+lcn
			print >> outh, longestSeq.replace("-","")
	else:
		for cluster in range(len(cl)):
			maxlen = 0
			longestSeq = ""
			for seq in cl[cluster]:
				if len(seq.replace("-","").replace("N","")) > maxlen:
					longestSeq = seq
			with open(lcn+"cl"+str(cluster)+".fas", "w") as outh:
				print >> outh, ">"+lcn+"cl"+str(cluster)
				print >> outh, longestSeq.replace("-","")


seqnum = len(seqlist)
locname = infilename.split("/")[-1].split(".")[0]
if seqnum == 1:
	clusters = [[seqlist[0]]]
elif seqnum == 2:
	dists1 = pairwise_dist(seqlist[0],seqlist[1])
	if dists1 == True:
		clusters = [[seqlist[0],seqlist[1]]]
	else:
		clusters = [[seqlist[0]],[seqlist[1]]]
	
else:
	comparisons = [(0,1), (0,2), (1,2)]
	dists = []
	for comp in comparisons:
		seq1a = seqlist[comp[0]]
		seq2a = seqlist[comp[1]]
		dists.append(pairwise_dist(seq1a, seq2a))
	if dists[0] == False and dists[1] == False:
		clusters = [[seqlist[0]], [seqlist[1],seqlist[2]]]
	elif dists[1] == False and dists[2] == False:
		clusters = [[seqlist[2]], [seqlist[0],seqlist[0]]]
	elif dists[0] == False and dists[2] == False:
		clusters = [[seqlist[1]], [seqlist[0],seqlist[2]]]
	elif dists[0] == False and dists[1] == False and dists[2] == False:
		clusters = [[seqlist[0]], [seqlist[1]], [seqlist[2]]]
	else:
		clusters = [[seqlist[0], seqlist[1], seqlist[2]]]

print locname, len(clusters)
fetch_fasta(clusters, locname)
		
