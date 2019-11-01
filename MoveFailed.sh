#!/bin/bash
pathfailed=/lustre/scratch115/projects/int_wes_metabol/raremetal/AllMetabolites/raremetal_inverse_norm_AllMet/FailedMet
pathMet=/lustre/scratch115/projects/int_wes_metabol/raremetal/AllMetabolites/raremetal_inverse_norm_AllMet/
declare -a Approaches=("naive" "functional" "LoF")
for Approach in "${Approaches[@]}"
do
     for i in $( cat ${pathfailed}/failed.txt ); do
           echo item: $i-${Approach}
           mv ${pathMet}/${Approach}/${i}  ${pathfailed}/${Approach}/          
       done
done

