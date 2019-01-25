import re
file = open('W303_reference.fasta')
line = file.readline()
stored_list = []
temp_number = file.tell()
fasta_dict = {}
while line:
	if line[0] == '>':
		temp_number = file.tell()
		stored_list.append((temp_number,re.sub(r'\W+','',line.rstrip('\n'))))
	else:
		if temp_number in fasta_dict:
			fasta_dict[temp_number] += line.rstrip('\n')
		else:
			fasta_dict[temp_number] = line.rstrip('\n')
	line = file.readline()

file.close()
'''
for fasta_l in stored_list:
	for fasta_d in fasta_dict.keys():
		if fasta_d == fasta_l[0]:
			fasta_dict[fasta_l[1]] = fasta_dict.pop(fasta_d)

for fasta_key,fasta_value in fasta_dict.items():
	f = open((fasta_key+'.fasta'), 'w')
	f.write(fasta_value)
	f.close()

'''
