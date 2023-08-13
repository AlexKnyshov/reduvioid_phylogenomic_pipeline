import sys
from Bio import AlignIO
from Bio.Align import AlignInfo

infilename = sys.argv[1]
locname = infilename.split("/")[-1]

with open(infilename) as infh:
	alignment = AlignIO.read(infh, "fasta")
	d = {}
	for seqnum in range(len(alignment)):
		taxon = alignment[seqnum].id.split("_")[0]
		if taxon in d:
			d[taxon].append(seqnum)
		else:
			d[taxon] = [seqnum]

with open(locname, "w") as outh:
	for k, v in d.items():
		if len(v) > 1:
			summary_align = AlignInfo.SummaryInfo(alignment[v[0]:v[len(v)-1]+1])
			print >> outh, ">"+k
			print >> outh, summary_align.dumb_consensus()
		else:
			print >> outh, ">"+k
			print >> outh, alignment[v[0]].seq