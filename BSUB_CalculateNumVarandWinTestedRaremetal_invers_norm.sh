#!/bin/bash

echo 'bash /nfs/users/nfs_l/lb17/My_Script_2016/int_wes_metabol_Scripts/CalculateNumVarandWinTestedRaremetal_invers_norm.sh  ${LSB_JOBINDEX}'  | bsub -J "Approach[1-3]" -o /nfs/users/nfs_l/lb17/My_Script_2016/int_wes_metabol_Scripts/LogJobs/CalculateNumVarandWinTestedRaremetal_invers_norm_APPR%I.out -R "select[mem>=5000] rusage[mem=5000]" -M5000 -q normal

