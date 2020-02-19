#!/usr/bin/python
# Script by Kiruthiga G. Mariappan
import sys
from Bio import SeqIO
infile=sys.argv[1]
print(infile)
infile_name=infile.split(".")[0]
outfile=open(infile_name+"_size.fasta","w")
for seq_record in SeqIO.parse(open(infile),"fasta"):
	seq_id=seq_record.id
	sequence=str(seq_record.seq)
	length=len(sequence)
	size="size"+str(length)
	new_seq_name="Pver_"+seq_id+"_"+size
	outfile.write(">"+new_seq_name+"\n"+sequence+"\n")
