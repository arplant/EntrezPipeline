# -*- coding: utf-8 -*-
"""

Information:
Entrez Pipeline
Retrieve protein fasta sequences based on search term
arplant
03/04/2023


Last update:
Defaulted to combining search result IDs into single set of unique IDs

To do:
Parse fasta into dictionary of ids and sequences
Alignment and partitioning

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
from Bio import SeqIO
import time
Entrez.email = 'anemail@adomain.com'
searchName = "lcw"

#Read input terms from file "input.txt" into a list
inFile = open("input.txt","r")
inputList = inFile.read().splitlines() #Split at linebreaks, does not include newline char in item
print(inputList)

#Create a set to store unique IDs
#This is used if different search terms are being used to retrieve the same type of protein
combinedIdList = set()

#Loop through the query terms, searching the database and returning IDs to the combined ID list
for query in inputList:
    queryterm = query + "[protein] NOT partial[Description]"
    handle = Entrez.esearch(db="protein", term=queryterm, retmax = 1, retmode="fasta")
    record = Entrez.read(handle)
    handle.close()
    for id in record["IdList"]:
        combinedIdList.add(id)
        print("Added " + id)

#The combined ID list can be used to retrieve sequence data, etc.
print(combinedIdList)

#Write the fasta for each ID to a file
proFasta = open(searchName + ".fa", "a")
for entry in combinedIdList:
    handle = Entrez.efetch(db="protein", id=entry, rettype="fasta", retmode="text")
    proFasta.write(handle.read())
    handle.close()
    time.sleep(0.3) #Manually control query rate
proFasta.close()

sequences = []
for seq_entry in SeqIO.parse(searchName + ".fa", "fasta"):
    sequences.append(str(seq_entry.seq))
    print(seq_entry.id)




"""
#Nested loop to find IDs based on search terms, then write the associated FASTA sequences to separate files
#For each input term, search the database and recover a list of IDs, then write matching data to file
for query in inputList:
    queryterm = query + "[protein] AND Bacteria[Orgn] NOT partial[Description]"
    handle = Entrez.esearch(db="protein", term=queryterm, retmax = 2, retmode="fasta")
    record = Entrez.read(handle)
    handle.close()
    print(record["IdList"])
    
    proFasta = open(query + ".fa", "a")
    #Working script to write fasta results to file iteratively
    for protId in record["IdList"]:
        handle = Entrez.efetch(db="protein", id=protId, rettype="fasta", retmode="text")
        proFasta.write(handle.read())
        handle.close()
        combinedIdList.add(protId)
        time.sleep(1) #Manually control query rate
    proFasta.close()


print(combinedIdList)

"""
