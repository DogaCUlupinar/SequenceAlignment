from operator import itemgetter

filename = "INS_PRAC.txt"
fileoutput = "GOOD_INS2.txt"
output = open(fileoutput,'w')
insertions = []
with open(filename) as f:
	for line in f:
		if ">" not in line:
			line = line.strip()
			ins,pos = line.split(",")
			insertions.append((ins,int(pos)))

insertions = sorted(insertions, key=itemgetter(1))


final_insertions = []
final_insertions.append(insertions[0])

for insertion in insertions:
	if insertion[1] != final_insertions[-1][1]:
		final_insertions.append(insertion)
	elif len(insertions[0]) > len(final_insertions[-1][0]):
		final_insertions[-1] = insertion

for inser in final_insertions:
	output.write(inser[0] + "," + str(inser[1]) +"\n")
output.close()

