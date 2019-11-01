#!/bin/bash

for i in {1..4} ;do
#i=5
echo 'Rscript /nfs/users/nfs_l/lb17/My_Script_2016/int_wes_metabol_Scripts/QQplot_AllMet.R --ID=${LSB_JOBINDEX}'  --Meth=${i} | bsub -J 'Approach[1-3]' -o /nfs/users/nfs_l/lb17/My_Script_2016/int_wes_metabol_Scripts/LogJobs/QQplot_AllMet_${i}_inverse_norm_Appr%I.out -R "select[mem>=40000] rusage[mem=40000]" -M40000 -q long

done
