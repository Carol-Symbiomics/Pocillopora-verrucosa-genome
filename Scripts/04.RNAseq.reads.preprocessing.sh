#!/bin/bash

# Remove low quality bases and adapter sequences 
java -jar trimmomatic-0.38.jar PE M_18_3799_Pver25_D701-D508_L006_R1_001.fastq.gz M_18_3799_Pver25_D701-D508_L006_R2_001.fastq.gz \
M_18_3799_Pver25.R1.paired.fq.gz M_18_3799_Pver25.R1.unpaired.fq.gz M_18_3799_Pver25.R2.paired.fq.gz M_18_3799_Pver25.R2.unpaired.fq.gz \
ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 > M_18_3799_Pver25_trim.log

# Remove PCR duplicates
prinseq-lite.pl -fastq M_18_3799_Pver25.R1.paired.fq -fastq2 M_18_3799_Pver25.R2.paired.fq -out_format 3 -log M_18_3799_Pver25_prinseq.log \
-graph_data M_18_3799_Pver25_R1_prinseq_graph_data.txt -derep 1 -derep_min 2 -graph_stats ld,gc,qd,de > M_18_3799_Pver25_prinseq.log

# Error Correction
spades.py -1 M_18_3799_Pver25_prinseq_good_MlLX.fastq -2 M_18_3799_Pver25_prinseq_good_wk12.fastq -o M_18_3799_Pver25_corrected --only-error-correction > M_18_3799_Pver25_error_correction.log

# Filter out reads less than 100bp in length
cutadapt -m 100 -o M_18_3799_Pver25_100bps.R1.fastq.gz -p M_18_3799_Pver25_100bps.R2.fastq.gz M_18_3799_Pver25_prinseq_good_wk12.00.0_0.cor.fastq.gz M_18_3799_Pver25_prinseq_good_MlLX.00.0_0.cor.fastq.gz
