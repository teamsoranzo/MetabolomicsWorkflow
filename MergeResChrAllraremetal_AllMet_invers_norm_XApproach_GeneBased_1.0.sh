#!/bin/bash
## declare an array variable
PHENO=/lustre/scratch115/projects/int_wes_metabol/Metabolon/CleanedData_160414/Pheno
## Change into working directory
raremetal=/lustre/scratch115/projects/int_wes_metabol/raremetal/AllMetabolites/raremetal_inverse_norm_AllMet
declare -a Approaches=("naive" "functional" "LoF")
j=`expr $1 - 1`
Approach=${Approaches[$j]}
echo $Approach
declare -a Meths=("SKAT_" "VT_" "burden" "MB")
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
				h1=`sed -n '4p' ${raremetal}/${Approach}/${MetInit}/Chr9_INTERVAL_QCv1.0_inverse_norm_${Approach}.MAF${MAF}.${MetInit}.meta.${Meth}.results  | cut -f2-`
				h2=`echo -e "WINDOW_NAME\tCHR\tSTART\tEND\tREGION\tGENE_NAME\tNVARWin"`
				echo -e  $h2" "$h1 > ${raremetal}/${Approach}/header${Meth}.tmp
				for trait in "${traits[@]}"
				do
					for chr in "${chrs[@]}"
					do
						echo VAR:${Approach}_${Meth}_${MAF}_${chr}_${trait}
						tail -n+5 ${raremetal}/${Approach}/${trait}/Chr${chr}_INTERVAL_QCv1.0_inverse_norm_${Approach}.MAF${MAF}.${trait}.meta.${Meth}.results | sort -k1,1 > ${raremetal}/${Approach}/${trait}/Chr${chr}_INTERVAL_QCv1.0_inverse_norm_${Approach}.MAF${MAF}.${trait}.meta.${Meth}.tmp
						sort -k5,5 /lustre/scratch115/projects/int_wes_metabol/VEP_85_GENCODE24_GFF3/Chr${chr}/${Approach}/Chr${chr}_windows_${Approach}_${MAF}.txt | sort -k5,5  > ${raremetal}/${Approach}/${trait}/2Chr${chr}_INTERVAL_QCv1.0_inverse_norm_${Approach}.MAF${MAF}.${trait}.meta.${Meth}.tmp
						join -1 5 -2 1 ${raremetal}/${Approach}/${trait}/2Chr${chr}_INTERVAL_QCv1.0_inverse_norm_${Approach}.MAF${MAF}.${trait}.meta.${Meth}.tmp  ${raremetal}/${Approach}/${trait}/Chr${chr}_INTERVAL_QCv1.0_inverse_norm_${Approach}.MAF${MAF}.${trait}.meta.${Meth}.tmp > ${raremetal}/${Approach}/${trait}/3Chr${chr}_INTERVAL_QCv1.0_inverse_norm_${Approach}.MAF${MAF}.${trait}.meta.${Meth}.tmp
						#cat  ${raremetal}/header.tmp ${raremetal}/${trait}/3.tmp > ${raremetal}/${trait}/Combined_Chr${chr}_INTERVAL_QCv1.0_naive.${trait}.meta.VT.results
						cat  ${raremetal}/${Approach}/${trait}/3Chr${chr}_INTERVAL_QCv1.0_inverse_norm_${Approach}.MAF${MAF}.${trait}.meta.${Meth}.tmp > ${raremetal}/${Approach}/${trait}/Combined_Chr${chr}_INTERVAL_QCv1.0_inverse_norm_${Approach}.MAF${MAF}.${trait}.meta.${Meth}.results
						rm -f ${raremetal}/${Approach}/${trait}/3Chr${chr}_INTERVAL_QCv1.0_inverse_norm_${Approach}.MAF${MAF}.${trait}.meta.${Meth}.tmp
						rm -f ${raremetal}/${Approach}/${trait}/Chr${chr}_INTERVAL_QCv1.0_inverse_norm_${Approach}.MAF${MAF}.${trait}.meta.${Meth}.tmp
						rm -f ${raremetal}/${Approach}/${trait}/2Chr${chr}_INTERVAL_QCv1.0_inverse_norm_${Approach}.MAF${MAF}.${trait}.meta.${Meth}.tmp
					done
					echo Done
					cat ${raremetal}/${Approach}/${trait}/Combined_Chr*_INTERVAL_QCv1.0_inverse_norm_${Approach}.MAF${MAF}.${trait}.meta.${Meth}.results | sort -n -k2 -k3  > ${raremetal}/${Approach}/${trait}/Combined_ALL_INTERVAL_QCv1.0_inverse_norm_${Approach}.MAF${MAF}.${trait}.meta.${Meth}.tmp
					cat ${raremetal}/${Approach}/header${Meth}.tmp ${raremetal}/${Approach}/${trait}/Combined_ALL_INTERVAL_QCv1.0_inverse_norm_${Approach}.MAF${MAF}.${trait}.meta.${Meth}.tmp > ${raremetal}/${Approach}/AllMerged/Combined_ALL_INTERVAL_QCv1.0_inverse_norm_${Approach}.MAF${MAF}.${trait}.meta.${Meth}.results
					rm -f ${raremetal}/${Approach}/${trait}/Combined_ALL_INTERVAL_QCv1.0_inverse_norm_${Approach}.MAF${MAF}.${trait}.meta.${Meth}.tmp
					echo Done All ${Approach}_${trait}_${Meth}
				done
			done
	done

