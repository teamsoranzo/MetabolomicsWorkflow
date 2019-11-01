#!/bin/bash
## declare an array variable
PHENO=/lustre/scratch115/projects/int_wes_metabol/Metabolon/CleanedData_160414/Pheno
## Change into working directory
raremetal=/lustre/scratch115/projects/int_wes_metabol/raremetal/AllMetabolites/raremetal_inverse_norm_AllMet
declare -a Approaches=("naive")
j=`expr $1 - 1`
Approach=${Approaches[$j]}
echo $Approach
declare -a Meths=("singlevar")
## now loop through the above array
declare -a chrs=({1..22} "X")
filename=${PHENO}/AllListMet_No11Fail_No2SegFault.txt
declare -a traits
traits=( `cat "$filename"`)
declare -a MAFs=("0.001")

        for Meth in "${Meths[@]}"
        do
			for MAF in "${MAFs[@]}"
                	do
                   		echo VARINT:${Approach}_${Meth}_${MAF}
				MetInit=M20675
				sed -n '4p' ${raremetal}/${Approach}/${MetInit}/Chr9_INTERVAL_QCv1.0_inverse_norm_${Approach}.MAF${MAF}.${MetInit}.meta.${Meth}.results | sed 's/#//' > ${raremetal}/${Approach}/header${Meth}.tmp
				for trait in "${traits[@]}"
				do
        				for chr in "${chrs[@]}"
        				do
        					echo Chr${chr}_${trait}_${Meth}
        					tail -n+5 ${raremetal}/${Approach}/${trait}/Chr${chr}_INTERVAL_QCv1.0_inverse_norm_${Approach}.MAF${MAF}.${trait}.meta.${Meth}.results | sed '/#Genomic/d'  > ${raremetal}/${Approach}/${trait}/Combined_Chr${chr}_INTERVAL_QCv1.0_inverse_norm_${Approach}.MAF${MAF}.${trait}.meta.${Meth}.results
       					done
					echo Done
					cat ${raremetal}/${Approach}/${trait}/Combined_Chr*_INTERVAL_QCv1.0_inverse_norm_${Approach}.MAF${MAF}.${trait}.meta.${Meth}.results | sort -V -k1,1 -k2,2  > ${raremetal}/${Approach}/${trait}/Combined_ALL_INTERVAL_QCv1.0_inverse_norm_${Approach}.MAF${MAF}.${trait}.meta.${Meth}.tmp
					cat ${raremetal}/${Approach}/header${Meth}.tmp ${raremetal}/${Approach}/${trait}/Combined_ALL_INTERVAL_QCv1.0_inverse_norm_${Approach}.MAF${MAF}.${trait}.meta.${Meth}.tmp  > ${raremetal}/Singlevar/Combined_ALL_INTERVAL_QCv1.0_inverse_norm_${Approach}.MAF${MAF}.${trait}.meta.${Meth}.results
					rm -f ${raremetal}/${Approach}/${trait}/Combined_ALL_INTERVAL_QCv1.0_inverse_norm_${Approach}.MAF${MAF}.${trait}.meta.${Meth}.tmp
					echo Done All ${Approach}_${trait}_${Meth}
				done
			done
	done

