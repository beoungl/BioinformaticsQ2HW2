import re


def open_file(file_name):
        file = open(file_name)
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
        return fasta_dict,stored_list


def key_change(fasta_dict,stored_list):
        for fasta_l in stored_list:
                for fasta_d in fasta_dict.keys():
                        if fasta_d == fasta_l[0]:
                                fasta_dict[fasta_l[1]] = fasta_dict.pop(fasta_d)

def make_files(fasta_dict):
        for fasta_key,fasta_value in fasta_dict.items():
                f = open((fasta_key+'.fasta'), 'w')
                f.write(fasta_value)
                f.close()

def search_contig(fasta_dict):
        contig = input('Which contig do you want to look at?:')
        if contig in fasta_dict.keys():
                print(contig,fasta_dict[contig])
if __name__ == '__main__':
        fasta_dict, fasta_list = open_file('W303_reference.fasta')
        key_change(fasta_dict,fasta_list)
        #To make files erase the #.
        #make_files(fasta_dict)
        #To search for certain contig erase the #.
        while True:
                search_contig(fasta_dict)
        
        
