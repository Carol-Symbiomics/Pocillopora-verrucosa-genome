#!/bin/bash

# Remove low quality bases and adapter sequences 
java -jar trimmomatic-0.36.jar \
PE -threads 20 ../00_Pver_Hiseq2500_500bp_shogun_rawfiles/180625_D00658_0052_AH7KYJBCX2/Lane2/version_01/M-18_1356_Pver25-lib4-size-seleted_NEBNEXT12_L002_R1_001.fastq.gz \
1exec_trim_new.sh/M-18_1356_Pver25-lib4-size-seleted_NEBNEXT12_L002_R2_001.fastq.gz Pver_lib4.R1.paired.fq.gz Pver_lib4.R1.unpaired.fq.gz Pver_lib4.R2.paired.fq.gz Pver_lib4.R2.unpaired.fq.gz \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:15:20 LEADING:25 TRAILING:20 MINLEN:36 > Pver_lib4_trim.log

# Filter out reads less than 100bp in length
cutadapt -j 10 -m 100 -o Pver_lib4_min100bp.R1.paired.fq.gz \
-p Pver_lib4_min100bp.R2.paired.fq.gz ../0_Adapter_lowqualitybp_remove/Pver_lib4.R1.paired.fq.gz ../0_Adapter_lowqualitybp_remove/Pver_lib4.R2.paired.fq.gz

# Quality assessment of reads after filtering 
fastqc -o FASTQC_after_trimming Pver_lib4_min100bp.R1.paired.fq.gz Pver_lib4_min100bp.R2.paired.fq.gz 


# Remove PCR duplicates using clumpify.sh from the BBmap tool set
bbmap_v38.51/clumpify.sh \
in=../1Trimmomatic/Pver_lib4_min100bp.R1.paired.fq.gz \
in2=../1Trimmomatic/Pver_lib4_min100bp.R2.paired.fq.gz \
out=Pver_lib4_min100bp_dup_removed.R1.paired.fq.gz \
out2=Pver_lib4_min100bp_dup_removed.R2.paired.fq.gz dedupe=t > Pver_lib4_min100bp_dup_removed_clumpify.log 


# Filter our potential contaminants from Symbiodiniaceae, bacteria and viruses
# For this purpose I used the following databases
# 1. All viral sequences (11890 sequences) downloaded from NCBI using tax id 10239
# 2. All bacterial complete genome sequences (14379 genomes) downloaded from NCBI using tax id 2 (Complete and latest genomes)
# 3. All Symbiodiniaceae sequences and genomes (766186 sequences) downloaded from NCBI using the tax id 252141

cd contaminant_sources
bbsplit.sh \
in=../Pver_lib4_min100bp_dup_removed.R1.paired.fq.gz \
in2=../Pver_lib4_min100bp_dup_removed.R2.paired.fq.gz \
outu1=../clean_bbsplit_R1.fq \
outu2=../clean_bbsplit_R2.fq \
basename=out_%.fq \
refstats=stats_bbsplit.txt \
threads=30 \
-Xmx500g &> Pver_bbsplit_20191107.log


# before contamination removal
zcat Pver_lib4_min100bp_dup_removed.R1.paired.fq.gz | echo $((`wc -l`/4))

# after contamination removal
cat clean_bbsplit_R1.fq | echo $((`wc -l`/4))

# Check ho many reads mapped to contaminant as follows:
cat out_Symbiodiniaceae_txid252141_nucleotide-genomes_20190821.fq | echo $((`wc -l`/4))
cat out_Bacteria_txid2_completegenomes1_20190821.fq | echo $((`wc -l`/4))
cat out_Bacteria_txid2_completegenomes2_20190821.fq | echo $((`wc -l`/4))
cat out_Viruses_txid10239_allseqs_20190822.fq | echo $((`wc -l`/4))

