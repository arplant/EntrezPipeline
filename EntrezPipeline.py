# -*- coding: utf-8 -*-
"""

Information:
Entrez Pipeline
Retrieve protein fasta sequences based on search term
arplant
03/04/2023


Last update:
Added multiple search term input from file.
Inputs are separate, one search per input, no check for duplication.


Acknowledgements:
Many thanks to leadbot and the following sources:
https://biopython.org/docs/dev/api/Bio.Entrez.html
http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec145
https://github.com/gumption/Using_Biopython_Entrez/blob/master/Biopython_Tutorial_and_Cookbook_Chapter_9.ipynb
"""

#Overview
#Read query terms from file, retrieve fasta, write the fasta to an output, then process it

#from Bio import SeqIO
#from io import BytesIO
#import os
from Bio import Entrez
import time
Entrez.email = 'alastairrobertplant@gmail.com'


#Read input terms from file "input.txt" into a list; input file contains one term per line
inFile = open("input.txt","r")
inputList = inFile.read().splitlines() #Split at linebreaks, does not include newline char in item
print(inputList)

#For each input term, search the database and recover a list of IDs, then write matching data to file
for query in inputList:
    queryterm = query + "[protein] AND Bacteria[Orgn] NOT partial[Description]"
    handle = Entrez.esearch(db="protein", term=queryterm, retmax = 1, retmode="fasta")
    record = Entrez.read(handle)
    handle.close()
    print(record["IdList"])
    
    proFasta = open(query + ".fa", "a")
    #Working script to write fasta results to file iteratively
    for protID in record["IdList"]:
        handle = Entrez.efetch(db="protein", id=protID, rettype="fasta", retmode="text")
        proFasta.write(handle.read())
        handle.close()
        time.sleep(3)#Manually control query rate to avoid use restriction
    proFasta.close()
