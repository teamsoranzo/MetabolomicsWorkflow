#!/bin/bash


PHENO=/lustre/scratch115/projects/int_wes_metabol/Metabolon/CleanedData_160414/Pheno
filename=${PHENO}/AllListMet.txt
filename2=${PHENO}/ListNum_AllMet_func.txt
awk '{print NR}' $filename > $filename2 
while read traits
do
echo 'bash /nfs/users/nfs_l/lb17/My_Script_2016/int_wes_metabol_Scripts/run_raremetal_INTERVAL_inverse_norm_functional_AllMet_1.0.sh' ${traits} | bsub -o /nfs/users/nfs_l/lb17/My_Script_2016/int_wes_metabol_Scripts/LogJobs/run_raremetal_INTERVAL_inverse_norm_functional_AllMet_1.0_${traits}_trueINVNORM.out -R "select[mem>=5000] rusage[mem=5000]" -M5000 -q long
done < $filename2
