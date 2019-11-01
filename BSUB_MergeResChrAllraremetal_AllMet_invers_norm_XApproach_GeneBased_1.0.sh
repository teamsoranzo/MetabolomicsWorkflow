#!/bin/bash

echo 'sh /nfs/users/nfs_l/lb17/My_Script_2016/int_wes_metabol_Scripts/MergeResChrAllraremetal_AllMet_invers_norm_XApproach_GeneBased_1.0.sh  ${LSB_JOBINDEX}' | bsub -J "Approach[1-3]" -o /nfs/users/nfs_l/lb17/My_Script_2016/int_wes_metabol_Scripts/LogJobs/MergeResChrAllraremetal_AllMet_invers_norm_XApproach_GeneBased_1.0_SecondRoundApproach%I.out -R "select[mem>=2000] rusage[mem=2000]" -M2000 -q long

