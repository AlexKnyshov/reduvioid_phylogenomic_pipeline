import sys
import glob

blastfile = sys.argv[2]

ranklineagefile = sys.argv[1]

rankd = {}

with open(ranklineagefile) as ranklineagehandle:
	for row in ranklineagehandle:
		row1 = row.strip().split("|")
		#rankd[row1[0].strip()] = [i.strip() for i in row1[1:]]
		rankd[row1[0].strip()] = row1[7].strip()
print "taxdump parsed"

for bfile in glob.glob(blastfile+"/*"):
	blastd = {}
	outhandle = open(bfile+"_decont_out.txt", "w")
	with open(bfile) as blasthandle:
		for row in blasthandle:
			row1 = row.strip().split()
			if row1[1] in rankd:
				#if rankd[row1[1]][6] == "Arthropoda":
				if row1[0] not in blastd:
					blastd[row1[0]] = False
				if rankd[row1[1]] == "Arthropoda":
					#print >> outhandle, row1[0], row1[1], "Arthropoda"
					#blastd[row1[0]] = rankd[row1[1]]
					blastd[row1[0]] = True
			# 	else:
			# 		#print >> outhandle, row1[0], row1[1], rankd[row1[1]][6], "exclude"
			# 		print >> outhandle, row1[0], row1[1], rankd[row1[1]], "exclude"
			# else:
			# 	print >> outhandle, row1[0], row1[1], "NOT FOUND"

			# if row1[0] in blastd:
			# 	blastd[row1[0]].add(row1[1])
			# else:
			# 	blastd[row1[0]] = set()
			# 	blastd[row1[0]].add(row1[1])
	#print blastd

	# for k, v in rankd.items():
	# 	print k, v[6]
	for k, v in blastd.items():
		if v == False:
	 		print >> outhandle, k#, v
	#print rankd
	outhandle.close()
	print bfile, "parsed"