#!/bin/bash

#for i in {1..4} ;do
i=1
echo 'Rscript /nfs/users/nfs_l/lb17/My_Script_2016/int_wes_metabol_Scripts/QQplot_AllMet_SingleVar_firsthalfMet.R --ID=${LSB_JOBINDEX}'  --Meth=${i} | bsub -J 'Approach[1]' -o /nfs/users/nfs_l/lb17/My_Script_2016/int_wes_metabol_Scripts/LogJobs/QQplot_AllMet_SingleVar_ALLMet_${i}_inverse_norm_Appr%I.out -R "select[mem>=60000] rusage[mem=60000]" -M60000 -q basement

#done
