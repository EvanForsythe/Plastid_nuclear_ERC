# makeSeqFiles.py
# script to read in a table (.csv) of sequence names, pull out the names in each row, find those
# seqs in a fasta file, and combine the seqs in a row into a fasta file

import sys
import csv
from Bio import SeqIO

#run from directory containing files or specify paths

def main():

	#arguments
	csvFilename = sys.argv[1] #name of csv file containing list of subtrees
	# csvFilename = "min_subtrees_forERC.csv"
	fastaFilename = sys.argv[2] #big file containing all proteomes
	# fastaFilename = "all_prots.fa"
	
	#open and read csv file
	f = open(csvFilename)
	csv_f = csv.reader(f)
	next(csv_f) #skip first line because it's just headers
	
	#open and read fasta file into a dictionary
	fastadict = SeqIO.to_dict(SeqIO.parse(fastaFilename, 'fasta'))
	
	#loop through and pull out sequences names for searching purposes
	#then get sequence out of dictionary and add to file
	for row in csv_f:
		newFilename = row[0] + "_seqs.fasta" #make a new file for each subtree
		newFile = open(newFilename, "w") #overwrites 
		for i in range(1, len(row)):
			if row[i] != "NA":
				#actually get sequences
				entry = fastadict[row[i]] #find that particular entry in the dictionary
				newFile.write(">" + row[i] + "\n") #identifier of fasta
				newFile.write(str(entry.seq) + "\n") #sequence
		newFile.close()
	f.close()
	
if __name__ == "__main__":
	main()