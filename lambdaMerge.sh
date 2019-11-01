#!/bin/bash
raremetal=/lustre/scratch115/projects/int_wes_metabol/raremetal/AllMetabolites/raremetal_inverse_norm_AllMet



declare -a Approaches=("naive" "functional" "LoF")
declare -a Meths=("SKAT_" "VT_" "burden" "MB")


for Approach in "${Approaches[@]}"
do
        for Meth in "${Meths[@]}"
        do
		cat ${raremetal}/${Approach}/QQplot/lambda_MetOrd_${Approach}_${Meth}* | sed '/Trait/d' > ${raremetal}/${Approach}/QQplot/lambda_MetOrd_${Approach}_${Meth}_ALL.txt
		echo ${Approach}"_"${Meth}
	done 
done 


