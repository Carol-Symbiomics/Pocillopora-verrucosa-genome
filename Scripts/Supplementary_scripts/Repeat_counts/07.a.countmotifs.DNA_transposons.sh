# Script by Larissa Morales Soto
genomeF=Number_of_filtered_repeats_in_Pvergenome/clean_5_Pver_lib4_filtered_HM2_ABD_size.fasta.ori.out
mwFile=Repbase_classif_system_Bao_et_al_2015/DNA_transposon
for f in $(ls ${mwFile})
do
	while read tid
	do
			echo "---------------------------------------------"
			echo "--> File: $f"
			echo "--> Searching motif ID $tid in the genome"
			mkdir -p motif_coordinates_DNA_transposon
			grep -w ${tid} ${genomeF} | cut -f7,8 > motif_coordinates_DNA_transposon/${tid}_${f}
			
	done < ${mwFile}/${f}
done
