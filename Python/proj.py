with open("/Users/allenchoy/Documents/project0.fq") as f:
	lines = []
	k = f.readline().rstrip()
	while ( k != ""):
		lines.append(k)
		k = f.readline().rstrip()
	f.close()

def average_quality_score(line:str)->float:
	s = 0
	for q in line:
		s += ord(q)-33 #ord(q): the ASCII code of q
	return(s/len(line))

new_f = open("new.fasta",'w') #new file name
for i in range(0,len(lines),4):
	name = lines[i][1:]
	GC_content = (lines[i+1].count("G")+lines[i+1].count("C"))/len(lines[i+1])*100 
	average_q_score = average_quality_score(lines[i+3])
	print("%s %s %f"%(name,GC_content,average_q_score))
	new_f.write(">%s\n"%(lines[i][1:])) #change fastq to fasta
	new_f.write(lines[i+1])
	new_f.write('\n')
new_f.close()
