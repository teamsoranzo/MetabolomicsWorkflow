#!/bin/bash

echo 'bash /nfs/users/nfs_l/lb17/My_Script_2016/int_wes_metabol_Scripts/lambdaMerge.sh'  | bsub -o /nfs/users/nfs_l/lb17/My_Script_2016/int_wes_metabol_Scripts/LogJobs/lambdaMerge_invers_norm.out -R "select[mem>=5000] rusage[mem=5000]" -M5000 -q normal
