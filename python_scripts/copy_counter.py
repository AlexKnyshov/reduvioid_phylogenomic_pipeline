import glob
import sys

infolder = sys.argv[1]

d = {}
s = set()
totalfiles = glob.glob(infolder+"/*")
print "reading loci"
total = float(len(totalfiles))
c = 0
c2 = 0
for f in totalfiles:
	fn = f.split("/")[-1][:-4]
	d[fn] = {}
	with open(f) as infh:
		for line in infh:
			splitline = line.strip().split()
			d[fn][splitline[0]] = splitline[5]
			s.add(splitline[0])
	c += 1
	pc = int(c / total *100)
	if pc - c2 == 10:
		print str(pc)+"%", c
		c2 = pc

print "writing data"
with open("output.csv","w") as outh:
	print >> outh, "loci,"+",".join(s)
	c = 0
	c2 = 0
	for loc, taxdata in d.items():
		printline = loc
		for tax in s:
			if tax in taxdata:
				printline = printline+","+taxdata[tax]
			else:
				printline = printline+",0"
		print >> outh, printline
		c += 1
		pc = int(c / total *100)
		if pc - c2 == 10:
			print str(pc)+"%", c
			c2 = pc
