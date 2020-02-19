#!/bin/bash

###

#(i) Build lmer table size 16. I've tested default l-mer size and it produces 16bp repeat frequence table
~/Genomes/Genome.tools/RepeatScout-1/build_lmer_table -l 16 -sequence 5_Pver_lib4_filtered_HM2_ABD_size.fasta -freq Pver_repeatscout.freq 

#(ii) RepeatScout using .freq file
~/Genomes/Genome.tools/RepeatScout-1/RepeatScout -sequence 5_Pver_lib4_filtered_HM2_ABD_size.fasta -output Pver_HaploidGenome_repeats.fasta -freq Pver_repeatscout.freq -l 16 &> Pver_RepeatScout20191123.log

#(iii) Filter stage 1 => remove repeats less than 50 nucleotides
cat Pver_HaploidGenome_repeats.fasta | ~/Genomes/Genome.tools/RepeatScout-1/filter-stage-1.prl > Pver_HaploidGenome_repeats_filtered_stg1.fasta

#(iv) Soft-Mask genome with repeats
~/Genomes/Genome.tools/RepeatMasker/RepeatMasker -xsmall -pa 25 -s -lib Pver_HaploidGenome_repeats_filtered_stg1.fasta 5_Pver_lib4_filtered_HM2_ABD_size.fasta -dir . &> Pver_RepeatSoftMasker20191123.log

# (v) Keep only the repeats that are occurring more than 10 times (default is 10).
cat Pver_HaploidGenome_repeats_filtered_stg1.fasta | ~/Genomes/Genome.tools/RepeatScout-1/filter-stage-2.prl --cat=5_Pver_lib4_filtered_HM2_ABD_size.fasta.out --thresh=10 > Pver_HaploidGenome_repeats_filtered_stg2_thresh10.fasta

###

# Repeats annotation - tblastx
tblastx -num_threads 25 -evalue 10e-3 -query Pver_HaploidGenome_repeats_filtered_stg2_thresh10.fasta -db ./Repeat_databases/RepBase24.02.fasta/RepbaseV24.02 -outfmt "6 delim=, qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen" -max_target_seqs 5 -out ./Repeat_tblastx_annotation/Pver_HaploidGenome_repeats_filtered_stg2_thresh10_vs_Repbase.out

# Repeats annotation blastx (custom-made database)
blastx -num_threads 25 -evalue 10e-3 -query Pver_HaploidGenome_repeats_filtered_stg2_thresh10.fasta -db ./Repeat_databases/nr_TE_20190922/TE_blastxdb_20190922 -out ./Repeat_blastx_annotation/Pver_HaploidGenome_repeats_filtered_stg2_thresh10_vs_TEdb.out -outfmt "6 delim=, qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen" -max_target_seqs 5

# Repeat annotation - Repeat masker
~/Genomes/Genome.tools/RepeatMasker/RepeatMasker -pa 20 -s -dir ./RepeatMasker_annotation -lib ./Repeat_databases/RepBase24.02.fasta/RepbaseV24.02.fasta Pver_HaploidGenome_repeats_filtered_stg2_thresh10.fasta > Repeatmasker_annotation_20191204.log

###

# Annotations present in the fasta header were added to the blast hit tables in R using the script "add.repeats.annotation.R"

# Repeats annotated using blastx
awk -F "," '{print $2}' reps_blastx_annot.csv | sort | uniq | sed 's/"//g' > listofrepeats_annotated_blastx.txt

# Repeats annotated using tblastx
awk -F "," '{print $2}' reps_tblastx_annot.csv | sort | uniq | sed 's/"//g' > listofrepeats_annotated_tblastx.txt

# The lists above are compared between them to identify the list of repeats that are uniquely identified by each strategy as well as those shared
comm -23 listofrepeats_annotated_blastx.txt listofrepeats_annotated_tblastx.txt | sed 's/^/"/g' | sed 's/$/"/g' > listofrepeats_annotated_blastx_only
comm -13 listofrepeats_annotated_blastx.txt listofrepeats_annotated_tblastx.txt | sed 's/^/"/g' | sed 's/$/"/g' > listofrepeats_annotated_tblastx_only
comm -12 listofrepeats_annotated_blastx.txt listofrepeats_annotated_tblastx.txt | sed 's/^/"/g' | sed 's/$/"/g' > listofrepeats_annotated_tblastx-blastx_shared

# Repeats in common and exclusively of blast and tblastx are identified and annotation is set apart
# blastx
grep -Ff listofrepeats_annotated_blastx_only reps_blastx_annot.csv > reps_blastx_annot_only.csv

# print best blastx hit only
cat reps_blastx_annot_only.csv | sort --field-separator=',' -k2,2 -k3,3r -k4,4r -k11,11 -k1,1 | sed 's/,/\t/g' | awk '!seen[$2]++' > reps_besthit_blastx_annot_only.txt

# tblastx
grep -Ff listofrepeats_annotated_tblastx_only reps_tblastx_annot.csv > reps_tblastx_annot_only.csv

# print best tblastx hit only
cat reps_tblastx_annot_only.csv | sort --field-separator=',' -k2,2 -k3,3r -k4,4r -k11,11 -k1,1 | sed 's/,/\t/g' | awk '!seen[$2]++' > reps_besthit_tblastx_annot_only.txt

# Repeats in common between blastx and tblastx

grep -Ff listofrepeats_annotated_tblastx-blastx_shared reps_blastx_annot.csv > reps_blastxtblastx_annot_shared.csv

# sort the hits based on pident and length and choose the best (1st per repeat)
# r in sort is used to reverse the order of the sorting to descending order

cat reps_blastxtblastx_annot_shared.csv | sort --field-separator=',' -k2,2 -k3,3r -k4,4r -k11,11 -k1,1 | sed 's/,/\t/g' | awk '!seen[$2]++' >  reps_besthit_blastx-tblasx_annot_shared.txt

# Combine the best blast hits in a single table to classiffy repeats by families. 
cat reps_besthit_blastx-tblasx_annot_shared.txt reps_besthit_blastx_annot_only.txt reps_besthit_tblastx_annot_only.txt  > best_blasthits_2874reps.txt

###

# Repeats annotated using RepeatMasker

awk '{print $5}' ./RepeatMasker_annotation/Pver_HaploidGenome_repeats_filtered_stg2_thresh10.fasta.out | sort | uniq | sed 's/\(^\|$\)/"/g' > RepeatMasker_annot_seqs
wc -l RepeatMasker_annot_seqs

# Repeats annotated as simple repeats
grep "Simple_repeat" ./RepeatMasker_annotation/Pver_HaploidGenome_repeats_filtered_stg2_thresh10.fasta.out | awk '{print $5, $11}' | sed 's/[^[:space:]]\+/"&"/g' | sort | uniq > RepeatMasker_1388Simple_Repeats_annotation.txt

awk '{print$1}' RepeatMasker_1388Simple_Repeats_annotation.txt > simple_repeats # This list is important to exclude it from the classified annotation 

awk '{print $5, $10}' ./RepeatMasker_annotation/Pver_HaploidGenome_repeats_filtered_stg2_thresh10.fasta.out | sort | uniq | awk '!seen[$1]++'| sed 's/[^[:space:]]\+/"&"/g'| grep -vFf simple_repeats > RepeatMasker_2513other_annotations.txt # print seqID and classification of RepeatMasker_annotation excluding those clasiffied as simple repeats (2 header lines removed manually)

awk '{print $2, $14}' best_blasthits_2874reps.txt | grep -vFf RepeatMasker_annot_seqs - > best_uniqblasthits_1056reps.txt # 1056 repeats annotated by blast that were not annotated with RepeatMasker

grep -c "\"\"" best_uniqblasthits_1056reps.txt # However, 42 of the uniq blast annotations didn't have a description assigned
#42

# concatenate all the Repeats annotated
cat best_uniqblasthits_1056reps.txt RepeatMasker_2513other_annotations.txt RepeatMasker_1388Simple_Repeats_annotation.txt > Pver_4958Repeats_annotated_20191213

### 

# Repeats counting

#Generate a masked file to count how many times the repeat motifs in Pver_HaploidGenome_repeats_filtered_stg2_thresh10.fasta are present in the genome

#(iv) Soft-Mask genome with repeats
~/Genomes/Genome.tools/RepeatMasker/RepeatMasker -xsmall -pa 25 -s -lib Pver_HaploidGenome_repeats_filtered_stg2_thresh10.fasta 5_Pver_lib4_filtered_HM2_ABD_size.fasta -dir ./Number_of_filtered_repeats_in_Pvergenome &> ./Number_of_filtered_repeats_in_Pvergenome/Pver_RepeatSoftMasker20200203.log

# Change NCBI keywords for specific code and remove motif ID without annotation description
sed 's/retrotransposon/OTL_NCBI_KW/g' Pver_4958Repeats_annotated_20191213 | sed 's/transposase/OTL_NCBI_KW/g' | sed 's/reverse_transcriptase/OTL_NCBI_KW/g' | sed 's/transposable_element/OTL_NCBI_KW/g' | grep -v -E '\"\"' | sort | uniq > Pver_4915Repeats_annotated_20200204 

# Some of the annotation descriptions were manually modified to facilitate the classification 

# Proceed to classify motifs using the Repbase Classification system by Bao et al 2015 in a hierarchical way

# 1. DNA transposons
06.a.find-match-DNA_transposons.sh # this classified 957 motifs into this group

# 2. LTR retrotransposon
# To prevent the incorrect classification of motifs, those previously assigned to DNA transposons will be removed from the searching list
cat ./Repbase_classif_system_Bao_et_al_2015/DNA_transposon/*.txt | grep -vwFf - Pver_4915Repeats_annotated_20200204 > Pver_3958Repeats_annotated_20200205

06.b.find-match-LTR_retrotransposons.sh # this classified 620 motifs into this group

# 3. non-LTR retrotransposon
# To prevent the incorrect classification of motifs, those previously assigned to LTR transposons will be removed from the searching list
cat Repbase_classif_system_Bao_et_al_2015/LTR_retrotransposon/*.txt | grep -vwFf - Pver_3958Repeats_annotated_20200205 > Pver_3338Repeats_annotated_20200205

06.c.find-match-non-LTR_retrotransposons.sh # this classified 849 motifs into this group

#In addition the search of SINE(s) have to be done independently due to the complication of grep to fullfil the requirement of match
grep "SINE" Pver_3338Repeats_annotated_20200205 | awk '{print $1}' |sed 's/"//g' > ./Repbase_classif_system_Bao_et_al_2015/Non-LTR_retrotransposon/SINE_repeats_list.txt
occ=$(wc -l ./Repbase_classif_system_Bao_et_al_2015/Non-LTR_retrotransposon/SINE_repeats_list.txt)
len=$(grep -wFf ./Repbase_classif_system_Bao_et_al_2015/Non-LTR_retrotransposon/SINE_repeats_list.txt Pver_HaploidGenome_repeats_filtered_stg2_thresh10_size.fasta | awk -F "size" '{ sum += $2} END {print sum}')
echo "SINE,$occ,$len" >> ./Repbase_classif_system_Bao_et_al_2015/Non-LTR_retrotransposon/Non-LTR_retrotransposon-annotated-report.csv
# 38 additional motifs were classified as SINE leaving 887 motifs assigned to this group

# 4. Simple repeats
# To prevent the incorrect classification of motifs, those previously assigned to Non-LTR transposons will be removed from the searching list
cat Repbase_classif_system_Bao_et_al_2015/Non-LTR_retrotransposon/*.txt | grep -vwFf - Pver_3338Repeats_annotated_20200205 > Pver_2451Repeats_annotated_20200205

06.d.find-match-Simple_repeats.sh # this classified 1401 motifs into this group

# 5. Multicopy genes
# To prevent the incorrect classification of motifs, those previously assigned to Simple_repeats will be removed from the searching list
cat Repbase_classif_system_Bao_et_al_2015/Simple_repeats/*.txt | grep -vwFf - Pver_2451Repeats_annotated_20200205 > Pver_1050Repeats_annotated_20200205

06.e.find-match-Multicopy_gene.sh # this classified 40 motifs into this group

# 6. Integrated_virus
# To prevent the incorrect classification of motifs, those previously assigned to Multicopy_genes will be removed from the searching list
cat Repbase_classif_system_Bao_et_al_2015/Multicopy_gene/*.txt | grep -vwFf - Pver_1050Repeats_annotated_20200205 > Pver_1010Repeats_annotated_20200205

06.f.find-match-Integrated_virus.sh # this classified 5 motifs into this group
# To prevent the incorrect classification of motifs, those previously assigned to Multicopy_genes will be removed from the searching list
cat Repbase_classif_system_Bao_et_al_2015/Integrated_virus/*.txt | grep -vwFf - Pver_1010Repeats_annotated_20200205 > Pver_1005Repeats_annotated_20200205

# 7. TEs annotated but not classified 
# Other-Non-LTR_retrotransposon >> 120 motifs
grep -iFf Repbase_classif_system_Bao_et_al_2015/Other-Non-LTR.txt Pver_1005Repeats_annotated_20200205 | awk '{print$1}' | sed 's/"//g' > Repbase_classif_system_Bao_et_al_2015/Other_TEs/Other-Non-LTR.txt

cat Repbase_classif_system_Bao_et_al_2015/Other_TEs/Other-Non-LTR.txt | grep -vwFf - Pver_1005Repeats_annotated_20200205 > Pver_885Repeats_annotated_20200205

# Other LTR_retrotransposon >> 11
grep "LTR" Pver_885Repeats_annotated_20200205|awk '{print$1}' | sed 's/"//g' > Repbase_classif_system_Bao_et_al_2015/Other_TEs/Other-LTR.txt

cat Repbase_classif_system_Bao_et_al_2015/Other_TEs/Other-LTR.txt | grep -vwFf - Pver_885Repeats_annotated_20200205 > Pver_874Repeats_annotated_20200205

# Other DNA transposons >> 223
grep -iFf Repbase_classif_system_Bao_et_al_2015/Other-DNA-transposon.txt Pver_874Repeats_annotated_20200205 | awk '{print$1}' | sed 's/"//g' > Repbase_classif_system_Bao_et_al_2015/Other_TEs/Other-DNA-transposon.txt

# Others >> 543
grep -viFf Repbase_classif_system_Bao_et_al_2015/Other_TEs/Other-DNA-transposon.txt Pver_874Repeats_annotated_20200205 | awk '{print$1}' | sed 's/"//g' > Repbase_classif_system_Bao_et_al_2015/Other_TEs/Other-TEs.txt

# Remove Empty files
find Repbase_classif_system_Bao_et_al_2015/ -size 0 -delete

# Now compile the coordinates (begining and end) of all the classified motifs to count its occurence and bps covered in the genome
# For this purpose use the scripts 07.a-g (find an example in the Supplementary scripts)

#check how many bp are covers by Unespecified/Unnanotated repeats
awk '{print$1}' Pver_4915Repeats_annotated_20200204 | sed 's/"//g' | grep -vwFf - ./Number_of_filtered_repeats_in_Pvergenome/clean_5_Pver_lib4_filtered_HM2_ABD_size.fasta.ori.out | grep "Unspecified" |awk '{print$10}' | sort | uniq | wc -l

awk '{print$1}' Pver_4915Repeats_annotated_20200204 | sed 's/"//g' | grep -vwFf - ./Number_of_filtered_repeats_in_Pvergenome/clean_5_Pver_lib4_filtered_HM2_ABD_size.fasta.ori.out | grep "Unspecified" | awk '{print$6, $7}' | sed 's/ /\t/g' > ./motif_coordinates_Unnanotated/Unnanotated_repeat_motifs.txt

# Sum up the bps covered by each of the repetitive elements class to generate a report

