#!/bin/bash

### 1. PASA ANNOTATION
## Clean transcript sequences with PASA script
seqclean Pver25_RL100bp_CL300bps_longest_Cnidarian_unique.fasta

## Align clean transcripts to the Genome
perl Launch_PASA_pipeline.pl -c alignAssembly.config -C -R -g 5_Pver_lib4_filtered_HM2_ABD_size.fasta -t Pver25_RL100bp_CL300bps_longest_Cnidarian_unique.fasta.clean -T -u Pver25_RL100bp_CL300bps_longest_Cnidarian_unique.fasta --ALIGNERS blat,gmap --CPU 30 --TRANSDECODER --transcribed_is_aligned_orient &> 2_launch_pasa_pipeline.log

## Call ORF
pasa_asmbls_to_training_set.dbi --pasa_transcripts_fasta Pver25_db.assemblies.fasta --pasa_transcripts_gff3 Pver25_db.pasa_assemblies.gff3 &> 3_convert_assembly_format.log

## Get a complete list of bonafide transcripts for training augustus.

# a) Identify ambiguous transcripts:
grep -P "5p2;|3p2;" Pver25_db.assemblies.fasta.transdecoder.genome.gff3 | awk -F "\t" '{print $9}' | awk -F ";" '{print $2}' | awk -F "=" '{print $2}' |sort -V | uniq > Pver25db_ambiguous_3prime_5prime_ids.txt

# b) Identify incomplete transcripts:
grep -P "\tgene\t" Pver25_db.assemblies.fasta.transdecoder.genome.gff3 | grep -v "complete" | awk -F "\t" '{print $9}' | awk -F ";" '{print $1}' | awk -F "~~" '{print $2}' | sort -V | uniq > Pver25db_incomplete_genes.txt

# c) Identify transcripts with less than 3 exons:
awk -F "\t" '{if($10<3) print $4}' Pver25_db.assemblies.fasta.transdecoder.genome.bed| awk -F ";" '{print $1}' | awk -F "=" '{print $2}' |sort -V |uniq > Pver25db_genes_with_less_than_3_exons.txt

# d) Combine 3 files:
cat Pver25db_ambiguous_3prime_5prime_ids.txt Pver25db_incomplete_genes.txt Pver25db_genes_with_less_than_3_exons.txt | sort | uniq > Pver25db_ambiguous_incomplete_less3exon.txt

# e) Filter out the transcripts identified above from the .bed file
	#create a list of the complete dataset of transcripts
grep "ID=asmbl" Pver25_db.assemblies.fasta.transdecoder.genome.bed | awk '{print $4}'| awk -F ";" '{print $1}' | sed 's/ID=//g' > Pver25db_allgenes.txt

grep -vFf Pver25db_ambiguous_incomplete_less3exon.txt Pver25db_allgenes.txt > Pver25db_unambiguous_complete_morethan3exons.txt # check that the file used to look for strings doesn't have spaces

grep -Ff Pver25db_unambiguous_complete_morethan3exons.txt Pver25_db.assemblies.fasta.transdecoder.genome.bed > Pver25db_unambiguous_complete_morethan3exons.bed

# f) Filter out redundant transcripts using the peptide sequences of the previously identified transcripts
	#create a fasta file from the .bed file to identify reduntant transcripts
bedtools getfasta -name -fi 5_Pver_lib4_filtered_HM2_ABD_size.fasta -bed Pver25db_unambiguous_complete_morethan3exons.bed -fo Pver25db_unambiguous_complete_morethan3exons.fasta
	#Translate the nucleotide sequences to amnoacids
~/useful_scripts/faTrans Pver25db_unambiguous_complete_morethan3exons.fasta Pver25db_unambiguous_complete_morethan3exons_aa.fasta
	#use the augustus perl script to BLAST all training transcript amino acid sequences against themselves and output
	#sequence ID had to be shortened in order to be compatible with BLAST
cat Pver25db_unambiguous_complete_morethan3exons_aa.fasta | awk -F ";" '{print $1}' > Pver25db_unambiguous_complete_morethan3exons_aa_shortID.fasta
	#only those protein sequences that are less than 80% redundant with any other sequence in the set
perl aa2nonred.pl --cores=30 Pver25db_unambiguous_complete_morethan3exons_aa_shortID.fasta Pver25db_unambiguous_complete_morethan3exons_aa_nonredundant.fasta

grep -c  ">" Pver25db_unambiguous_complete_morethan3exons_aa_nonredundant.fasta

# g) Filter out potential repeats within the set of transcripts previously identified
	# extract the nucleotide sequences corresponding to the non-redudant proteins using a list
	# beware that fasta headers have to match exactly, therefore, shorten the header of the nucleotide file
cat Pver25db_unambiguous_complete_morethan3exons.fasta | awk -F ";" '{print $1}' > Pver25db_unambiguous_complete_morethan3exons_shortID.fasta

~/useful_scripts/faSomeRecords Pver25db_unambiguous_complete_morethan3exons_shortID.fasta list_of_unambiguous_complete_morethan3exons_aa_nonredundant.txt Pver25db_unambiguous_complete_morethan3exons_nonredundant_nucl.fasta

blastn -query Pver25db_unambiguous_complete_morethan3exons_nonredundant_nucl.fasta -db RepeatsDB/PverRepeatScoutDB -evalue 10e-3 -outfmt 6 -num_threads 30 > Pver25db.transdecoder.genome_complete_morethan3exons_unambiguous_nonred_vs_Pverlib4RepeatScout.out

awk '{print $1}' Pver25db.transdecoder.genome_complete_morethan3exons_unambiguous_nonred_vs_Pverlib4RepeatScout.out | sort | uniq > list_of_genes_hitting_repeats.txt

~/useful_scripts/faSomeRecords -exclude Pver25db_unambiguous_complete_morethan3exons_nonredundant_nucl.fasta list_of_genes_hitting_repeats.txt Pver25db_unambiguous_complete_morethan3exons_nonred_norepeats.fasta

# h) Extract the list of bonafide transcripts for augustus training from the gff3 file
grep ">" Pver25db_unambiguous_complete_morethan3exons_nonred_norepeats.fasta | sed 's/>ID=//g' > list_of_bonafide_genes_4augustustraining.txt

grep -Ff list_of_bonafide_genes_4augustustraining.txt Pver25_db.assemblies.fasta.transdecoder.genome.gff3 > Pver25db_bonafide_transcripts.gff3


## Compute flanking regions
perl computeFlankingRegion.pl Pver25db_bonafide_transcripts.gtf > Pver25db_bonafide_transcript_flanking_regions.txt

## Convert genome file and GenomeThreader gtf training gene file to GenBank flatfile format to train AUGUSTUS
perl gff2gbSmallDNA.pl Pver25db_bonafide_transcripts.gtf 5_Pver_lib4_filtered_HM2_ABD_size.fasta 2275 Pver25db_bonafide_genes.gb


### 2. AUGUSTUS TRAINING

# Create a new species in the augustus config file
perl new_species.pl --species=pocillopora_verrucosa

etraining --species=pocillopora_verrucosa Pver25db_bonafide_genes.gb &> bonafide.out

etraining --species=pocillopora_verrucosa Pver25db_bonafide_genes.gb  2>&1 | grep "in sequence" | perl -pe 's/.*n sequence (\S+):.*/$1/' | sort -u > bad.lst # 6 genes seems to be problematic

# Filter out problematic genes
perl filterGenes.pl bad.lst Pver25db_bonafide_genes.gb > Pver_bonafide_genes_f.gb

# Split the training gene dataset into "test" and "training"
perl randomSplit.pl Pver_bonafide_genes_f.gb 200

# Re train augustus with the training set
etraining --species=pocillopora_verrucosa Pver_bonafide_genes_f.gb.train &> etrain.out

# Run Augustus with the newly trained parameters
augustus --species=pocillopora_verrucosa Pver_bonafide_genes_f.gb.test &> test.out

# Optimize the augustus parameters
perl optimize_augustus.pl --species=pocillopora_verrucosa --kfold=10 --cpus=30 Pver_bonafide_genes_f.gb.train > optimize.out

# Perform a final training with the training set using the parameters after optimization
etraining --species=pocillopora_verrucosa Pver_bonafide_genes_f.gb.train &> etrain_after_optimization.out

# Measure the accuracy (again) on the test set using the parameters after optimization
augustus --species=pocillopora_verrucosa Pver_bonafide_genes_f.gb.test &> test_after_optimization.out

# Create a hints file to run augustus

# Align ESTs/cDNAs to the softmasked genome using blat (pblat is the parallelized version of blat)
pblat -noHead  -threads=30 5_Pver_lib4_filtered_HM2_ABD_size.fasta Pver25_RL100bp_CL300bps_longest_Cnidarian_unique.fasta PvergenomeDB_vs_PverTransQuery.psl

# Filter alignments to obtain those that are potentially most useful for gene prediction
cat PvergenomeDB_vs_PverTransQuery.psl | ~/Genomes/Genome.tools/Augustus/scripts/filterPSL.pl --best --minCover=80 > PvergenomeDB_vs_PverTransQuery.est.f.psl

# Sort filtered alignments according to start position and target sequence name
cat PvergenomeDB_vs_PverTransQuery.est.f.psl | sort -n -k 16,16 | sort -s -k 14,14 > PvergenomeDB_vs_PverTransQuery.est.fs.psl

# Convert psl-file to hints
perl blat2hints.pl  --in=PvergenomeDB_vs_PverTransQuery.est.fs.psl --out=est.hints --minintronlen=45 --trunkSS

# UTR training
# Measure initial accuracy without UTR parameters
augustus --species=pocillopora_verrucosa Pver_bonafide_utr_f.gb.test > test.noUtr.out

# Optimize UTR parameters
# To notice is that the perl moduled Parallel::ForkManager was used to allow parallelization
# The afore mentioned module was installed in perl /home/linuxbrew/.linuxbrew/bin/perl
perl optimize_augustus.pl --species=pocillopora_verrucosa --rounds=3 Pver_bonafide_utr_f.gb.train --cpus=10 --kfold=10 --UTR=on --metapars=$AUGUSTUS_CONFIG_PATH/species/pocillopora_verrucosa/pocillopora_verrucosa_metapars.utr.cfg --trainOnlyUtr=1 > optimize.utr.out

# Test accuracy with UTR parameters
augustus --species=pocillopora_verrucosa Pver_bonafide_utr_f.gb.test --UTR=on --print_utr=on > test.utr.out

# Split the genome file in smaller files
mkdir split_genome
perl splitMfasta.pl 5_Pver_lib4_filtered_HM2_ABD_size.fasta --outputpath=split_genome --minsize=1000000

# Get the list of sequences in each chunk of genome to split the hints as well
mkdir split_genome_lists

for ((i=1; i<=328; i++)); do fgrep ">" split_genome/5_Pver_lib4_filtered_HM2_ABD_size.split.$i.fa | perl -pe 's/^>//' > ./split_genome_lists/part.$i.lst ;done

# Split hints
mkdir split_hints

cd split_genome_lists

for ((i=1; i<=328; i++)); do cat ../est.hints | perl getLinesMatching.pl part.$i.lst 1 > ../split_hints/est.split.$i.hints; done


# Multithreading prediction with AUGUSTUS
mkdir augustus_prediction_splitted
seq 1 328 | parallel -j 35 --bar --no-notice "nice augustus --species=pocillopora_verrucosa --strand=both --singlestrand=false --hintsfile=split_hints/est.split.{}.hints --alternatives-from-evidence=true --alternatives-from-sampling=true --sample=100 --minmeanexonintronprob=0.5 --maxtracks=2 --extrinsicCfgFile=extrinsic.pver.E.cfg --progress=true --UTR=on --print_utr=on --exonnames=on --codingseq=on split_genome/5_Pver_lib4_filtered_HM2_ABD_size.split.{}.fa > ./augustus_prediction_splitted/augustus.{}.out 2> ./augustus_prediction_splitted/augustus.{}.err"

#counting how many gens were predicted
cd augustus_prediction_splitted
for ((i=1; i<=328; i++)); do grep -c "# start gene" augustus.$i.out &>> number_of_genes_predicted.txt; done
awk '{ sum += $1 } END { print sum }' number_of_genes_predicted.txt #26,852 

# Merge AUGUSTUS output
cd ..
for ((i=1; i<=328; i++)); do cat ./augustus_prediction_splitted/augustus.$i.out ; done > concatenated_augustus_out.gff

perl join_aug_pred.pl < concatenated_augustus_out.gff > joined_augustus_Pver_predicted_genes.gff

grep -c "# start gene" joined_augustus_Pver_predicted_genes.gff #26,852

perl gtf2gff.pl < joined_augustus_Pver_predicted_genes.gff --gff3 --out=joined_augustus_Pver_predicted_genes.gff3 --printExon

grep -P -c "\tgene\t" joined_augustus_Pver_predicted_genes.gff3 #26,852

# note that AUGUSTUS also predicts alternative splice forms (ending e.g. in .t2)
# Typically one or more transcripts are transcribed from a locus. These have .t1, t2 etc. appended to the locus name e.g. a locus that expresses two alternative spliceforms might be described following the transcript IDs


### 3. PASA UPDATE

# Update the database and then the gene predictions
Load_Current_Gene_Annotations.dbi -c alignAssembly.config -g 5_Pver_lib4_filtered_HM2_ABD_size.fasta -P joined_augustus_Pver_predicted_genes.gff3

perl Launch_PASA_pipeline.pl -c alignAssembly.config -A -g 5_Pver_lib4_filtered_HM2_ABD_size.fasta -t Pver25_RL100bp_CL300bps_longest_Cnidarian_unique.fasta.clean --TRANSDECODER

# How many genes resulted after updation
grep -c -P "\tgene\t" *.gff3
#joined_augustus_Pver_predicted_genes.gff3:26852
#Pver25_db.gene_structures_post_PASA_updates.45243.gff3:27529

# Get the genes' name and isoforms list
grep ">" Pver25_db.gene_structures_post_PASA_updates.45243.gff3.prot.fasta | sed 's/>//g' > gene_names_and_isoforms.txt

# Get the number of ORFs
perl pasa_asmbls_to_training_set.extract_reference_orfs.pl Pver25_db.gene_structures_post_PASA_updates.45243.gff3 minProtLength=0 > Pver25_db.gene_structures_post_PASA_updates.45243.orfs

# how many ORFs resulted
grep -c -P "\tgene\t" *.orfs #26,634

# The prefix Pver and Pver_gene were added to the gff3 file using vim> Pver25_db.gene_structures_post_PASA_updates.45243.names.gff3

# The protein sequence was extracted from the names.gff3 file as follows
grep "#PROT" Pver25_db.gene_structures_post_PASA_updates.45243.names.gff3 | sed 's/#PROT />/g' | sed 's/\t/\n/g' > Pver25_db.gene_structures_post_PASA_updates.45243.gff3.prot.names.fasta

# The .names.gff3 file was copied into a different folder cds_mrna_exons and the following Augustus perl script was used to extract the condingseq (nt)
perl getAnnoFasta.pl Pver25_db.gene_structures_post_PASA_updates.45243.names.gff3 --seqfile=5_Pver_lib4_filtered_HM2_ABD_size.fasta
sed 's/>/>Pver_/g' Pver25_db.gene_structures_post_PASA_updates.45243.names3.codingseq > temp && mv temp Pver25_db.gene_structures_post_PASA_updates.45243.names3.codingseq


### 4. ASSESS COMPLETENESS OF THE PREDICTED GENE MODELS WITH BUSCO

# Run BUSCO in protein mode
python2.7 run_BUSCO.py -i Pver25_db.gene_structures_post_PASA_updates.45243.gff3.prot.fasta -o Pver_GeneModels_BUSCO -m proteins -c 20 &>> 12_BUSCO_results_20191203.log


### 5. CHECK RNAseq EVIDENCE FOR THE PREDICTED GENE MODELS

# Remove duplicates from .bed file
sort Pver25_db.gene_structures_post_PASA_updates.45243.bed | uniq > Pver25_db.gene_structures_post_PASA_updates.45243.bed.noDuplicates

# get nucleotides from .bed file
bedtools getfasta -fi 5_Pver_lib4_filtered_HM2_ABD_size.fasta -bed Pver25_db.gene_structures_post_PASA_updates.45243.bed.noDuplicates -name -s > Pver25_db.gene_structures_post_PASA_updates.45243.noDup.nuc.fasta

# Identify RNAseq evidence for the set of Gene models predicted
# Many RNAseq reads do not get assembled into transcript, however providing evidence to predicted gene models

# Indexing the Gene models fasta dataset
cd Pver_Gene_Models_index
bowtie2-build ../Pver25_db.gene_structures_post_PASA_updates.45243.noDup.nuc.fasta Pver_Gene_Models

# Align RNAseq read to Gene Models
cd ..
bowtie2 -p 25 -x Pver_Gene_Models_index/Pver_Gene_Models -1 M_18_3799_Pver25.R1.paired.fq.gz -2 M_18_3799_Pver25.R2.paired.fq.gz -S Pver25RNAseqTrimmedReads_vs_PverGenesModels.sam --no-unal  > Pver25RNAseqTrimmedReads_vs_PverGenesModels.log

# Covert sam to sorted bam files
samtools view -bS Pver25RNAseqTrimmedReads_vs_PverGenesModels.sam | samtools sort -o Pver25RNAseqTrimmedReads_vs_PverGenesModels.sort.bam -

# Check how many RNAseq reads aligned to each GeneModel
samtools idxstats Pver25RNAseqTrimmedReads_vs_PverGenesModels.sort.bam > Alignment_stats.20191204

#samtools idxstats produces a TAB delimited output where each column consists of ref.seq. name, ref.seq. length, number of mapped reads and number of unmapped reads.

# How many gene ID's (All the transcripts per gene were considered)
wc -l Alignment_stats.20191204
#32856 Alignment_stats.20191204 # substract the last column

# How many gene ID's have 0 reads mapping to it
awk '$3==0 {print$0}' Alignment_stats.20191204 | wc -l # 3,196  Only 9.72% of the Gene Models didn't have any RNAseq evidence

#90.27% of the Gene Models have RNAseq evidence


### 5. COMPLETE ANNOTATION BY BLASTING GENE MODELS TO Swiss, trEMBL AND nr DATASETS

# a) Build BLAST databases for swiss, TrEMBL and "nr" protein sequences

makeblastdb -in ./swiss_db/uniprot_sprot.fasta -parse_seqids -dbtype prot -out ./swiss_db/swissprot_20191205

makeblastdb -in ./trembl_db/uniprot_trembl.fasta -parse_seqids -dbtype prot -out ./trembl_db/trembl_20191101

perl update_blastdb.pl --blastdb_version 5 ./nr_db_20200121/nr --decompress

# b) BLAST against swiss 
blastp -max_target_seqs 5 -num_threads 20 -db swiss_db/swissprot_20191205 -query Pver25_db.gene_structures_post_PASA_updates.45243.gff3.prot.fasta -evalue 1e-5 -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' -out PverGeneModels_vs_sprot_1e-5_max5.out

# Get the best hit for each Gene Model (protein) swiss
cat PverGeneModels_vs_sprot_1e-5_max5.out | sort -k1,1 -k2,2 -k3,3r -k4,4r -k11,11 | awk '!seen[$1]++' > PverGeneModels_vs_sprot_1e-5_besthit.out

wc -l PverGeneModels_vs_sprot_1e-5_besthit.out #23653 

# Select the gene model proteins without hits in the swiss prot
awk '{print $1}' PverGeneModels_vs_sprot_1e-5_besthit.out > list_of_Pvergenemodelproteins_sprot.txt

~/useful_scripts/faSomeRecords -exclude Pver25_db.gene_structures_post_PASA_updates.45243.gff3.prot.fasta list_of_Pvergenemodelproteins_sprot.txt Pver25_db.gene_structures_post_PASA_updates.45243.gff3.prot4trembl.fasta

# c) BLAST the remaining protein sequences against TrEMBL

blastp -max_target_seqs 5 -num_threads 20 -db trembl_db/trembl_20191101 -query Pver25_db.gene_structures_post_PASA_updates.45243.gff3.prot4trembl.fasta -evalue 1e-5 -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' -out PverGeneModels_vs_trembl_1e-5_max5.out

# Get the best hit for each Gene Model (protein) TrEMBL
cat PverGeneModels_vs_trembl_1e-5_max5.out | sort -k1,1 -k2,2 -k3,3r -k4,4r -k11,11 | awk '!seen[$1]++' > PverGeneModels_vs_trembl_1e-5_besthit.out

wc -l PverGeneModels_vs_trembl_1e-5_besthit.out #8262 

# Select the gene model proteins without hits in TrEMBL
awk '{print $1}' PverGeneModels_vs_trembl_1e-5_besthit.out > list_of_Pvergenemodelproteins_trembl.txt

~/useful_scripts/faSomeRecords -exclude Pver25_db.gene_structures_post_PASA_updates.45243.gff3.prot4trembl.fasta list_of_Pvergenemodelproteins_trembl.txt Pver25_db.gene_structures_post_PASA_updates.45243.gff3.prot4nr.fasta

# d) BLAST the remaining protein sequences against nr

blastp -max_target_seqs 5 -num_threads 20 -db nr_db_20200121/nr -query Pver25_db.gene_structures_post_PASA_updates.45243.gff3.prot4nr.fasta -evalue 1e-5 -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' -out PverGeneModels_vs_nr_1e-5_max5.out

# Get the best hit for each Gene Model (protein) nr
cat PverGeneModels_vs_nr_1e-5_max5.out | sort -k1,1 -k2,2 -k3,3r -k4,4r -k11,11 | awk '!seen[$1]++' > PverGeneModels_vs_nr_1e-5_besthit.out

# description and GOterms were added to the list of hits using an R script and the most recent GOterms description  ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT
