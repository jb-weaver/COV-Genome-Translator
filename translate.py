genbank_file = open("MN908947.genbank", "r") #open genbank file and read lines
genbank_lines = genbank_file.readlines()
genbank_file.close()

fasta_file = open("MN908947.fasta", "r") #open fasta file and read lines
fasta_lines = fasta_file.readlines()
del fasta_lines[0] #remove title from fasta
remove_char = [fasta_lines.rstrip('\n') for fasta_lines in fasta_lines] 
string_sequence =  ''.join(remove_char)
sequence = string_sequence[:13468] + 'C' + string_sequence[13468:] #Accounts for ribosomal frameshift 


gene_library = open("Gene_Library.txt", "w")
CDS = open("DNA.fasta",'w')

coordinates = []
#Loop through each line and print lines thats matches the text pattern 'FEATURES'
for line in genbank_lines:
    line=line.strip() #remove newline character
    line=line.replace('join', '')
    line=line.replace('(','')
    line=line.replace(')','')
    line=line.replace('..',':') 
  #  print(line) #for each list in all the lines
    if line.startswith('gene'): #see if the first item in each list is equal to CDS
        line=line.replace('gene','')
        line=line.replace(' ','')
        coordinates.append(line)

dupe_names = []
dedupe_names = []

for line in genbank_lines:
    line = line.strip() 
    if line.startswith('/gene'): 
        line = line.replace('/gene=','')
        line = line.replace('"','')
        dupe_names.append('>' + line)
         
for i in dupe_names: 
    if i not in dedupe_names: 
        dedupe_names.append(i)

dictionary = dict(zip(dedupe_names,coordinates))
#print(dictionary)
data = str(dictionary)

start_str = []
stop_str = []
start_int = []
stop_int = []

#Cuts out the CDS coordinates and assigns them to start/stop values 
#Values can't be read as indices because they are currently strings, not integers 
#Values can't be converted to integers yet because they are joined by ':'
for key in dictionary: 
    unclean_coordinates = dictionary.get(key)
    unclean_coordinates = unclean_coordinates.split(":")
    start_str.append(unclean_coordinates[0]) 
    stop_str.append(unclean_coordinates[1])

#Convert start and stop coords into integers
for i in start_str: 
    i = int(i)
    if i == 266:
        i = i-1 #Python starts counting at 0, subtracting one accounts for this 
    start_int.append(i) 

for i in stop_str:
    i = int(i)
    if i != 21556:
        i +=1
    stop_int.append(int(i))

#Create dictionary with Start coords as dict Keys and Stop coords as dict Values  
clean_coordinates = dict(zip(start_int,stop_int)) #make dictionary of just int coordinates
print(clean_coordinates)
DNA = ''

#Parse DNA seq with coordinate and tag them with the appropiate gene names 
counter = 0
for i in clean_coordinates: 
    DNA = sequence[i:clean_coordinates[i]] #use key and value from clean_coordinates as string indices
    CDS.write('\n' + dedupe_names[counter] + '\n') #use counter to simulatneously loop through dedupe_names
    CDS.write(DNA)
    counter = counter + 1
    
CDS.close()
gene_library.close()

#Converting CDS DNA to RNA by converting 'T' to 'U'
RNA_file = open('RNA.fasta', 'w')

DNA_file = open('DNA.fasta','r')
DNA_lines = DNA_file.readlines()

sequence =  ''.join(DNA_lines) #lines are orignally read as lists, this makes each CDS sequence all one string 

RNA = sequence.replace("T","U")
RNA_file.write(RNA)


RNA_file.close()

RNA_file = open("RNA.fasta", "r+")
RNA_lines = RNA_file.readlines()
RNA_sequence = RNA_lines
#print(RNA_sequence)
Protein_file = open('Protein.fasta','w')

RNA_file.close()

###Translates RNA sequence into a Protein sequence
def translate(RNA):
    rna2protein = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
     "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
     "UAU":"Y", "UAC":"Y", "UAA":"", "UAG":"",
     "UGU":"C", "UGC":"C", "UGA":"", "UGG":"W",
     "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
     "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
     "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
     "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
     "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
     "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
     "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
     "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
     "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
     "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
     "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
     "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G"}

    dict(rna2protein)
    protein = ""
    if len(line)%3 == 0: #"If the sequence is evenly divisble by 3" then: 
        for i in range(0, len(line), 3): #For every three characters in a given sequence. AKA for each codon 
            codon = line[i:i+3] #Let every three character equal 'codon'
            protein += rna2protein[codon] #Then find the rna2protein value for 'codon' and add it to 'protein'
    return protein

#For each RNA sequence: 1)Write the Gene name 2)Translate the RNA sequence beneath it 3)Move on and repeat
for line in RNA_sequence:
    line = line.rstrip('\n') #get rid of new-line character so lines are more readable to 'if' statements
    if line.startswith(">"): #Find line containing Gene name
        Protein_file.write(line + '\n')
    elif translate(RNA_sequence) == "": #Prevents an empty line being written on the first line in Protein_file
        Protein_file.write("")
    else:
        Protein_file.write(translate(RNA_sequence) + '\n') #Translation function for RNA sequence 
        print(translate(RNA_sequence))


Protein_file.close()
