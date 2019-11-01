#!/bin/bash


PHENO=/lustre/scratch115/projects/int_wes_metabol/Metabolon/CleanedData_160414/Pheno
filename=${PHENO}/AllListMet.txt
filename2=${PHENO}/ListNum_AllMet_naive.txt
awk '{print NR}' $filename > $filename2 
while read traits
do
echo 'bash /nfs/users/nfs_l/lb17/My_Script_2016/int_wes_metabol_Scripts/GzipIndexRaremetalworkerWESresults.sh' ${traits} | bsub -o /nfs/users/nfs_l/lb17/My_Script_2016/int_wes_metabol_Scripts/LogJobs/GzipIndexRaremetalworkerWESresults_${traits}.out -R "select[mem>=5000] rusage[mem=5000]" -M5000 -q normal
done < $filename2

