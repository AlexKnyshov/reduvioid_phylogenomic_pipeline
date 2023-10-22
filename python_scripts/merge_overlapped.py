import sys

d = {}
with open (sys.argv[1]) as infile:
	for line in infile:
		splitline = line.split()
		if splitline[0] in d:
			if splitline[1] not in d[splitline[0]]:
				d[splitline[0]].add(splitline[1])
		elif splitline[1] in d:
			if splitline[0] not in d[splitline[1]]:
				d[splitline[1]].add(splitline[0])
		else:
			d[splitline[0]] = set()
			d[splitline[0]].add(splitline[1])

for k, v in d.items():
	print k, v
