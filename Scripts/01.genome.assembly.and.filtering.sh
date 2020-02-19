#!/bin/bash

## 1. Genome assembly of Pocillopora verrucosa

DiscovarDeNovo \
READS=./clean_bbsplit_R1.fq,./clean_bbsplit_R2.fq \
OUT_DIR=Pver_lib4_DISCOVAR_Assembly_Results \
NUM_THREADS= 40 \
MAX_MEM_GB=400 \
MEMORY_CHECK=TRUE > Pver_lib4_DISCOVAR_Assembly_Results_20191108.log

## 2. Filtering of the raw genome assembly

# a) Identify and remove circular contigs
grep "circular" a.lines.fasta | wc -l

~/useful_scripts/faSomeRecords -exclude a.lines.fasta ./01_circular_scaffolds/circular_scaffolds_raw_DISCOVAR_assembly.txt 1_Pver_lib4_nocircular.fasta
# faSomeRecords can be downloaded from http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/

# Add scaffold size to each sequence header in a fasta file
python2.7 ~/useful_scripts/calc_size_general.py 1_Pver_lib4_nocircular.fasta
rm 1_Pver_lib4_nocircular.fasta

# Identify and contigs were the begining or end are represented by string of Ns
grep -P -B1 "^N" 1_Pver_lib4_nocircular_size.fasta
grep -P -B1 "N$" 1_Pver_lib4_nocircular_size.fasta


# b) Remove overlooked vector and mitocondrial sequences
# Identify overlooked adaptor sequences identified with VecScreen (https://www.ncbi.nlm.nih.gov/tools/vecscreen/). Using same parameters as in the NCBI

blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true -evalue 700 -searchsp 1750000000000 -db ./Univec_db/UniVec_blastdb_20190717 -query 1_Pver_lib4_nocircular_size.fasta -out ./02_UniVec_blast_hits/Pver_lib4_vs_UniVecdb_evalue700.txt -num_threads 30 -outfmt "6 delim=, qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen"

# Looking at the vector db of ncbi ftp://ftp.ncbi.nlm.nih.gov/blast/db/
blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true -evalue 700 -searchsp 1750000000000 -db ./Univec_db/vector_db/vector -query 1_Pver_lib4_nocircular_size.fasta -out ./02_UniVec_blast_hits/Pver_lib4_vs_Vecdb_evalue700.txt -num_threads 30 -outfmt "6 delim=, qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen"

blastn -query 1_Pver_lib4_nocircular_size.fasta -db ./mito_seqs_db/mitochondrion.1.1.2.1.genomic_blastdb_20190717 -out ./02_mito_blast_hits/Pver_lib4_vs_mitodb_evalue1e-10_qlen.txt -evalue 1e-10  -word_size 15 -max_target_seqs 100 -gapopen 5 -gapextend 2 -penalty -3 -reward 2 -soft_masking true -num_threads 30 -outfmt "6 delim=, qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen"

# c) Identify and remove contigs containing potential contaminants that remain after the preprocessing of reads

# Blast search of conting against database symbiodinaceae, bacteria and viruses
blastn -query 2_Pver_lib4_nocircular_novec_nomito.fasta -db ./Sym_Bac_Vir_db/Symb_Bacteria_Viruses_blastdb_20190828 -out ./03_Potential_contaminat_blast_hits/Pver_nomito_novec_vs_contaminant_blastnhits.txt -evalue 1e-10 -max_target_seqs 10 -qcov_hsp_perc 50 -perc_identity 90 -num_threads 30 -outfmt "6 delim=, qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"

# print contig that hit the refence with an identity cut-off of 90% over 50% of its sequence length
awk -F "," '{print $1}' Pver_nomito_novec_vs_contaminant_blastnhits.txt | sort | uniq > Pver_blastncontaminants_contigstoremove.txt

# remove contigs with contaminant sequences
~/useful_scripts/faSomeRecords -exclude 2_Pver_lib4_nocircular_novec_nomito.fasta ./03_Potential_contaminat_blast_hits/Pver_blastncontaminants_contigstoremove.txt 3_Pver_lib4_nocircular_novec_nomito_nocontam.fasta


## 3. Reference guided scaffolding 

php csar.php -t ../3_Pver_lib4_nocircular_novec_nomito_nocontam.fasta -r ./Pdam_Scaffolds/pdam_scaffolds.fasta --pro -o 4_Pver_lib4_filtered_CSARefscaffold.fa &> Pver_lib4_FilteredDiscovar_CsarScaffolding_20191111.log

## 4. Gapfilling of scaffolded genome

perl -I Perl4-CoreLibs-0.004/lib GapFiller_v1-10_linux-x86_64/GapFiller.pl -l Pver_library4.txt -s ../4_Pver_lib4_filtered_CSARefscaffold.fa/scaffolds.pro.csar.fna -m 50 -o 3 -r 0.7 -n 10 -t 5 -d 50 -T 30 -i 3 -b Pver_lib4_filteredDiscovar_ScaffCsar_Gapfilledlib4m50 &> Pver_lib4_GapFillerm30_20191112.log

## 5. Rebuild haploid assembly by identifying oversplitted loci. Useful in highly heterozygous genomes

# a) construct a score matrix by splitting the genome fasta file in two portions with roughly 10% and 90% of the total sequences, respectively

grep ">" ../4_Pver_lib4_filteredDiscovar_ScaffCsar_Gapfilledlib4m50_size.fa | awk -F ">" '{print $2}' | head -63 > Pver_10percent_genome_seqs.txt

~/useful_scripts/faSomeRecords ../4_Pver_lib4_filteredDiscovar_ScaffCsar_Gapfilledlib4m50_size.fa Pver_10percent_genome_seqs.txt Pver_10percent_genome_seqs.fa

~/useful_scripts/faSomeRecords -exclude ../4_Pver_lib4_filteredDiscovar_ScaffCsar_Gapfilledlib4m50_size.fa Pver_10percent_genome_seqs.txt Pver_90percent_genome_seqs.fa

lastz_D_Wrapper.pl --target=Pver_10percent_genome_seqs.fa --query=Pver_90percent_genome_seqs.fa --identity=90 &>> Pver_score_matrix_20191112.log

# b) run the module A, B and D of Haplomerger2 using the custome-made script as suggested by the developers

bash run_ABD.batch # After running the pipeline A_B_D of haplomerger the resultant haploid reference genome is 380,505,698 bp

# Add scaffold size to each sequence header in a fasta file
python2.7 ./calc_size_general.py 5_Pver_lib4_filtered_HM2_ABD.fa


## 6. Assess the genome completeness using BUSCO with the metazoa_odb9 dataset 
python2.7 run_BUSCO.py -i ../5_Pver_lib4_filtered_HM2_ABD_size.fasta -o Pver_HaplogenomeHM2 -m geno -c 30 &>> 05_BUSCO_results_20191113.log



