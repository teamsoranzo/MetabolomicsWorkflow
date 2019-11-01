#!/bin/bash
echo 'Rscript /nfs/users/nfs_l/lb17/My_Script_2016/int_wes_metabol_Scripts/MatchGenesGWAS.R --ID=${LSB_JOBINDEX}' | bsub -J 'Approach[1,2,11]' -o /nfs/users/nfs_l/lb17/My_Script_2016/int_wes_metabol_Scripts/LogJobs/MatchGenesGWAS_inverse_norm_10-4_CommonANDrare_longq_Appr%I.out -R "select[mem>=4000] rusage[mem=4000]" -M4000 -q long

