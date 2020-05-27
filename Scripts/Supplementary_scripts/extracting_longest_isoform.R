###############################################################################################
######## Anny's Rscript version: sorting transcripts by longest isoform. 12.03.2020 ###########
## Usage:  Rscript extracting_longest_isoform.R protein.fasta codingseq.fasta annotation.gff ##
###############################################################################################

library(phylotools)
library(data.table)
library(seqinr)
library(reshape2)

#filtering proteins
args<-commandArgs(TRUE)
prot<-phylotools::read.fasta(args[1])
prot[c('isoform', 'gene')] <-colsplit(prot$seq.name,"\t| ",c('isoform', 'gene'))
prot$gene<-gsub("split", "Pver", prot$gene)
prot$gene<-gsub("\\-.*", "", prot$gene)
prot$length<-nchar(prot$seq.text)
filtered.prot.headers=setDT(prot)[, .SD[which.max(length)], gene]
message( "Total number of isoforms: ", length(unique(prot$isoform)), "\nNumber of isoforms retained: ", length(unique(filtered.prot.headers$isoform)))
write.fasta(as.list(filtered.prot.headers$seq.text), filtered.prot.headers$isoform, "longest_isoform.faa", open = "w", nbchar = 60, as.string = FALSE)
write.table(filtered.prot.headers$isoform, "selected_transcripts.txt", quote = FALSE, row.names = FALSE)

#filtering transcripts
trans<-phylotools::read.fasta(args[2])
trans.filtered<-subset(trans, trans$seq.name %in% filtered.prot.headers$isoform)
write.fasta(as.list(trans.filtered$seq.text), trans.filtered$seq.name, "longest_isoform.fna", open = "w", nbchar = 60, as.string = FALSE)

## QC for transcripts
#refe=read.table("/Users/anny/Documents/Projects/Carol_RNASeq/Genome/Final_transcripts2.txt", header = FALSE)
#match=subset(refe, refe$V1 %in% filtered.prot.headers$isoform )
#miss=subset(refe, !refe$V1 %in% filtered.prot.headers$isoform )
# 264 correspond to equally-long isoforms 
