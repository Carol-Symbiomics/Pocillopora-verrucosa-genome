#!/bin/bash

## 1. Transcriptome assembly of Pocillopora verrucosa
Trinity --seqType fq --max_memory 75G --left M_18_3799_Pver25_100bps.R1.fastq.gz --right M_18_3799_Pver25_100bps.R2.fastq.gz \
--output M_18_3799_Pver25_RL100bp_CL300bps_trinity_results_RF --CPU 16 --SS_lib_type RF --min_kmer_cov 5 --normalize_reads \
--min_contig_length 300 > M_18_3799_Pver25_RLmin100bps_CLmin300bp_trinity_RF.log

## 2. Filtering of the raw transcriptome assembly

# a) Select the longest isoform per transcript

# b) blast search against a holobiont custom-made database (NCBI cnidarian taxon id:6073; NCBI Symbiodiniaceae taxon id: 252141; NCBI viral taxon id: 10239 and 14,379 complete bacterial genomes)

megablast -d nt_db_Coral_Symb_Bac_Vir_blastdb_20190915 -i Pver25_RL100bp_CL300bps_longest_Trinity.fasta -e 1e-3 -m 8 -p 90 -b 1 -a 22 -o megablast_Pver25_RL100bp_CL300bps_longest_vs_AllSeqsDB_1e-3_b1.tab > megablast_LongestIsoTranscript_besthit_20190917.log

# get a list of the best megablast transcripts hit to each taxon ID
awk -F "," '{print $1,$2,$4,$5,$12,$14,$15}' megablast_transcriptbesthit_annotation_qaccsorted.csv | awk -F " " '{if ($6=="Cnidarian") print $2}' | sort | uniq   > cnidarian.bestmegablasthit.list

awk -F "," '{print $1,$2,$4,$5,$12,$14,$15}' megablast_transcriptbesthit_annotation_qaccsorted.csv | awk -F " " '{if ($6=="Symbiodiniaceae") print $2}' | sort | uniq   > symbiodiniaceae.bestmegablasthit.list

awk -F "," '{print $1,$2,$4,$5,$12,$14,$15}' megablast_transcriptbesthit_annotation_qaccsorted.csv | awk -F " " '{if ($6=="Bacteria") print $2}' | sort | uniq   > bacteria.bestmegablasthit.list

awk -F "," '{print $1,$2,$4,$5,$12,$14,$15}' megablast_transcriptbesthit_annotation_qaccsorted.csv | awk -F " " '{if ($6=="Virus") print $2}' | sort | uniq   > virus.bestmegablasthit.list

# Extract sequences from fasta file using the list of trasncript hits unique to cnidarians

~/useful_scripts/faSomeRecords Pver25_RL100bp_CL300bps_longest_Trinity.fasta cnidarian.bestmegablasthit.list Pver25_RL100bp_CL300bps_longest_Cnidarian_unique.fasta
# faSomeRecords can be downloaded from http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/

~/useful_scripts/faSomeRecords Pver25_RL100bp_CL300bps_longest_Trinity.fasta symbiodiniaceae.bestmegablasthit.list Pver25_RL100bp_CL300bps_longest_Symbiodinium_unique.fasta

## 3. Assess the cnidarian trasncriptome completeness using BUSCO with the metazoa_odb9 dataset
python2.7 /home/buitracn/Genomes/Genome.tools/busco/scripts/run_BUSCO.py -i Pver25_RL100bp_CL300bps_longest_Cnidarian_unique.fasta -o Pver25_cnidarian_clean_transcriptome_BUSCO \
-m transcriptome -c 20 &>> 01_BUSCO_results_20191203.log

## 4. Assess the Symbiodiniaceae trasncriptome completeness using BUSCO with the viridiplantae_odb10 
python2.7 /home/buitracn/Genomes/Genome.tools/busco/scripts/run_BUSCO.py -i Pver25_RL100bp_CL300bps_longest_Symbiodinium_unique.fasta -o Pver25_symbiodiniaceae_clean_transcriptome_BUSCO \
-m transcriptome -c 20 &>> 01_BUSCO_results_20200210.log