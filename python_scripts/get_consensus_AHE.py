import sys
from Bio import AlignIO
from Bio.Align import AlignInfo

infilename = sys.argv[1]
locname = infilename.split("/")[-1].split(".")[0]

with open(infilename) as infh:
	alignment = AlignIO.read(infh, "fasta")
	summary_align = AlignInfo.SummaryInfo(alignment)
	with open(locname+".fas", "w") as outh:
		print >> outh, ">"+locname
		print >> outh, summary_align.dumb_consensus()
