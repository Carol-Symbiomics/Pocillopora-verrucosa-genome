#Repeats annotation

setwd("/home/buitracn/Genomes/corals/Pocillopora.verrucosa/02_Pver_genome_repeats")

#tblastx
tblastx <- read.csv("Repeat_tblastx_annotation/Pver_HaploidGenome_repeats_filtered_stg2_thresh10_vs_Repbase.out", header = FALSE)
colnames(tblastx) <- c("qaccver", "saccver", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

repbase.annot <-  read.csv("Repeat_databases/RepBase24.02.fasta/Repbase24.02.annot.csv", header = FALSE)
colnames(repbase.annot) <- c("accessionID", "annot")

tblastx.annot <- merge(tblastx, repbase.annot, by.x=c("saccver"), by.y=c("accessionID"), sort = FALSE, all.x=TRUE)

write.csv(tblastx.annot, file="reps_tblastx_annot.csv", row.names = FALSE)

#blastx
blastx <- read.csv("Repeat_blastx_annotation/Pver_HaploidGenome_repeats_filtered_stg2_thresh10_vs_TEdb.out", header = FALSE)
colnames(blastx) <- c("qaccver", "saccver", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen")

nr.db.annot <- read.csv("Repeat_databases/nr_TE_20190922/nt_TE_20190922_annot.csv", header = FALSE)
colnames(nr.db.annot) <- c("accessionID", "annot")

blastx.annot <- merge(blastx, nr.db.annot, by.x=c("saccver"), by.y=c("accessionID"), sort = FALSE, all.x=TRUE)
write.csv(blastx.annot, file="reps_blastx_annot.csv", row.names = FALSE)
