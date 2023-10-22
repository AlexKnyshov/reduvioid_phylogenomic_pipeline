import sys

tablename = sys.argv[1]
filterdir = sys.argv[2]
outdir = sys.argv[3]

print "parse filter"
samplename = ".".join(tablename.split("/")[-1].split(".")[:-1])
filtername = filterdir+"/"+samplename+".blast_decont_out.txt"
filterset = set()
with open(filtername) as filterhandle:
	for line in filterhandle:
		filterset.add(line.strip())
print "filter parsed, parse table"
with open(outdir+"/"+samplename+".hmmer", "w") as outhandle:
	with open(tablename) as tablehandle:
		for line in tablehandle:
			if line[0] != "#":
				ctgname = "_".join(line.split()[0].split("_")[:-1])
				if ctgname not in filterset:
					print >> outhandle, line.strip()
print "done"