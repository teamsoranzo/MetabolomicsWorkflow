##Start Loop
raremetal=/lustre/scratch115/projects/int_wes_metabol/raremetal/AllMetabolites/raremetal_inverse_norm_AllMet
PHENO=/lustre/scratch115/projects/int_wes_metabol/Metabolon/CleanedData_160414/Pheno

declare -a Approaches=("naive" "functional" "LoF")
declare -a Meths=("VT_" "burden")
filename=${PHENO}/AllListMet_No11Fail_No2SegFault.txt
declare -a traits
traits=( `cat "$filename"`)
j=`expr $1 - 1`
Approach=${Approaches[$j]}
echo $Approach

mkdir -p ${raremetal}/${Approach}/StatsAfterAnalysis


        for Meth in "${Meths[@]}"
        do
		echo -e "trait\tNVar" > ${raremetal}/${Approach}/StatsAfterAnalysis/Nvars${Approach}_${Meth}.txt
		echo -e "trait\tNWin" > ${raremetal}/${Approach}/StatsAfterAnalysis/NWin${Approach}_${Meth}.txt
                        for trait in "${traits[@]}"
                        do
                                echo VARINT:${Approach}_${Meth}_${trait} 
				echo ${trait} > ${raremetal}/${Approach}/StatsAfterAnalysis/tmp.txt
				echo ${trait} > ${raremetal}/${Approach}/StatsAfterAnalysis/tmp1.txt				
				awk '$8 >= 5'  ${raremetal}/${Approach}/AllMerged/Combined_ALL_INTERVAL_QCv1.0_inverse_norm_${Approach}.MAF0.001.${trait}.meta.${Meth}.results | awk '{print $9}' | sed '/VAR/d' | perl -npe 's/;/\n/g' | wc -l | paste ${raremetal}/${Approach}/StatsAfterAnalysis/tmp.txt -  >> ${raremetal}/${Approach}/StatsAfterAnalysis/Nvars${Approach}_${Meth}.txt
				awk '$8 >= 5'  ${raremetal}/${Approach}/AllMerged/Combined_ALL_INTERVAL_QCv1.0_inverse_norm_${Approach}.MAF0.001.${trait}.meta.${Meth}.results | wc -l | paste ${raremetal}/${Approach}/StatsAfterAnalysis/tmp1.txt - >> ${raremetal}/${Approach}/StatsAfterAnalysis/NWin${Approach}_${Meth}.txt

			done
done


