#!/bin/bash


PHENO=/lustre/scratch115/projects/int_wes_metabol/Metabolon/CleanedData_160414/Pheno
## Change into working directory
rareMet_S=/nfs/team151/software/raremetal_v4.14/Raremetal-master/bin/raremetal
# Change the name of teh folder based on the MAF
raremetal=/lustre/scratch115/projects/int_wes_metabol/raremetal/AllMetabolites/raremetal_inverse_norm_AllMet/naive/DropVariantInWin
rareMetW=/lustre/scratch115/projects/int_wes_metabol/raremetal/AllMetabolites/raremetalworker_inverse_norm
signres=/lustre/scratch115/projects/int_wes_metabol/raremetal/AllMetabolites/raremetal_inverse_norm_AllMet/raremetal_sign_all_var.txt


#head -n1  ${signres}  | perl -npe 's/\t/\n/g' | cat -n

mkdir -p $raremetal
cd $raremetal

declare -a traits
traits=(`tail -n+2 ${signres} |  awk -F"\t" '$25=="naive" {print $19}'  | sort | uniq`)
j=`expr $1 - 1`
#j=1
trait=${traits[$j]}
#trait=M01589
echo $trait
mkdir -p  $trait
cd $trait

wins=(`awk -F"\t" -v m=$trait '$25=="naive" && $19==m {print $6}' ${signres} | sort | uniq`)

for win in "${wins[@]}"
do
	# Set variables
	echo $win	
#	awk -F"\t" -v m=$trait -v w=$win '$25=="naive" && $19==m && $6==w {print $6,$1,$4,$7}' ${signres} | sort -g -k3,3 
	chr=`awk -F"\t" -v m=$trait -v w=$win '$25=="naive" && $19==m && $6==w {print $7}' ${signres} | sort | uniq `
        echo $chr
	sta=`awk -F"\t" -v m=$trait -v w=$win '$25=="naive" && $19==m && $6==w {print $8}' ${signres} | sort | uniq `
	echo $sta
	end=`awk -F"\t" -v m=$trait -v w=$win '$25=="naive" && $19==m && $6==w {print $9}' ${signres} | sort | uniq `
	echo $end
	bioc=`awk -F"\t" -v m=$trait -v w=$win '$25=="naive" && $19==m && $6==w {print $21}' ${signres} | sort | uniq `
	echo $bioc
	MAF="0.001"
	echo $MAF
	awk -F"\t"  -v w=$win '$1==w'  /lustre/scratch115/projects/int_wes_metabol/VEP_85_GENCODE24_GFF3/Chr${chr}/naive/*${MAF}.list > group${win}.list
        group=group${win}.list        
	echo $group
        score=${rareMetW}/Chr${chr}/Chr${chr}_INTERVAL_QCv1.0_inverse_norm.${trait}.singlevar.score.txt.list
        echo $score
        cov=${rareMetW}/Chr${chr}/Chr${chr}_INTERVAL_QCv1.0_inverse_norm.${trait}.singlevar.cov.txt.list
        echo $cov
	ipval=`awk -F"\t" -v m=$trait -v w=$win -v t=$test '$25=="naive" && $19==m && $6==w && $20==t {print $18}' ${signres} | sort | uniq `
	echo $ipval
	#GW threshold naive
	#GWTnaive=2.402737e-08
	#GW threshold functional
	#GWTfunctional=3.842134e-08
	#GW threshold naive
	#GWTnaive=1.325838e-07

	tests=(`awk -F"\t" -v m=$trait -v w=$win '$25=="naive" && $19==m && $6==w {print $20}' ${signres} | sort | uniq`)
	for test in "${tests[@]}"
	do
		echo $test
		
		len=`awk '{print NF}' $group`
		for var in $(seq 0 $len) 
		do 
		echo $var
			if [ $(($len-$var)) -gt 1 ]
                	then
				cut --complement -f$(($len-$var)) $group > group${win}.${var}.list
				group_var=group${win}.${var}.list
				droponevar=`cut -f$(($len-$var)) $group`
				# run RAREMETAL
#                        	${rareMet_S} --summaryFiles $score --covFiles $cov --groupFile $group_var --maf ${MAF} --SKAT --burden --MB --VT  --longOutput --hwe 1e-06 --tabulateHits --hitsCutoff 1e-04  --prefix Chr${chr}_INTERVAL_QCv1.0_inverse_norm_naive.MAF${MAF}.${trait}.${win}.${var}.dropone
				
				if [ "${test}" == "SKAT_" ] || [ "${test}" == "VT_" ]
                		then

				ipval=`awk -F"\t" -v m=$trait -v w=$win -v t=$test '$25=="naive" && $19==m && $6==w && $20==t {print $18}' ${signres} | sort | uniq `				
				test_vars=`awk '{print $3}' Chr${chr}_INTERVAL_QCv1.0_inverse_norm_naive.MAF${MAF}.${trait}.${win}.${var}.dropone.meta.${test}.results | tail -n1 `
				drop_effectsize=`awk '{print $10}' Chr${chr}_INTERVAL_QCv1.0_inverse_norm_naive.MAF${MAF}.${trait}.${win}.${var}.dropone.meta.${test}.results | tail -n1`
				drop_pval=`awk '{print $12}' Chr${chr}_INTERVAL_QCv1.0_inverse_norm_naive.MAF${MAF}.${trait}.${win}.${var}.dropone.meta.${test}.results | tail -n1`
				echo -e "window\tchr\tstart\tend\ttrait\tbioc\ttest\tapproach\tinitial_pval\tnumber_vars\tdropped_var\ttest_vars\tdrop_effectsize\tdrop_pval" > Drop_chr"${chr}"_INTERVAL_QCv1.0_inverse_norm_naive.MAF"${MAF}"."${trait}"."${win}"."${var}"."${test}".txt

				echo -e "$win\t$chr\t$sta\t$end\t$trait\t$bioc\t$test\tnaive\t$ipval\t$var\t$droponevar\t$test_vars\t$drop_effectsize\t$drop_pval"  >> Drop_chr"${chr}"_INTERVAL_QCv1.0_inverse_norm_naive.MAF"${MAF}"."${trait}"."${win}"."${var}"."${test}".txt

				elif [ "${test}" == "MB" ] || [ "${test}" == "burden" ]
				then
                		ipval=`awk -F"\t" -v m=$trait -v w=$win -v t=$test '$25=="naive" && $19==m && $6==w && $20==t {print $18}' ${signres} | sort | uniq `
                                test_vars=`awk '{print $3}' Chr${chr}_INTERVAL_QCv1.0_inverse_norm_naive.MAF${MAF}.${trait}.${win}.${var}.dropone.meta.${test}.results | tail -n1 `
                                drop_effectsize=`awk '{print $10}' Chr${chr}_INTERVAL_QCv1.0_inverse_norm_naive.MAF${MAF}.${trait}.${win}.${var}.dropone.meta.${test}.results | tail -n1`
                                drop_pval=`awk '{print $11}' Chr${chr}_INTERVAL_QCv1.0_inverse_norm_naive.MAF${MAF}.${trait}.${win}.${var}.dropone.meta.${test}.results | tail -n1`
                                echo -e "window\tchr\tstart\tend\ttrait\tbioc\ttest\tapproach\tinitial_pval\tnumber_vars\tdropped_var\ttest_vars\tdrop_effectsize\tdrop_pval" > Drop_chr"${chr}"_INTERVAL_QCv1.0_inverse_norm_naive.MAF"${MAF}"."${trait}"."${win}"."${var}"."${test}".txt

                                echo -e "$win\t$chr\t$sta\t$end\t$trait\t$bioc\t$test\tnaive\t$ipval\t$var\t$droponevar\t$test_vars\t$drop_effectsize\t$drop_pval"  >> Drop_chr"${chr}"_INTERVAL_QCv1.0_inverse_norm_naive.MAF"${MAF}"."${trait}"."${win}"."${var}"."${test}".txt

				fi	

			fi
		done

		echo -e "window\tchr\tstart\tend\ttrait\tbioc\ttest\tapproach\tinitial_pval\tnumber_vars\tdropped_var\ttest_vars_inWin\tdrop_effectsizes\tdrop_pvals" > Drop2m_chr"${chr}"_INTERVAL_QCv1.0_inverse_norm_naive.MAF"${MAF}"."${trait}"."${win}"."${test}".txt
		ls Drop_chr"${chr}"_INTERVAL_QCv1.0_inverse_norm_naive.MAF"${MAF}"."${trait}"."${win}".*."${test}".txt | xargs -I {}  tail -n+2 {} | tail -n1 | cut -f1-9 > tmp1

			ls Drop_chr"${chr}"_INTERVAL_QCv1.0_inverse_norm_naive.MAF"${MAF}"."${trait}"."${win}".*."${test}".txt | xargs -I {}  tail -n+2 {} | awk '{print $(NF-4)}' | perl -npe 's/\n/,/' | sed 's/.$//' > tmp2
			ls Drop_chr"${chr}"_INTERVAL_QCv1.0_inverse_norm_naive.MAF"${MAF}"."${trait}"."${win}".*."${test}".txt | xargs -I {}  tail -n+2 {} | awk '{print $(NF-3)}' | perl -npe 's/\n/,/' | sed 's/.$//'> tmp3
			ls Drop_chr"${chr}"_INTERVAL_QCv1.0_inverse_norm_naive.MAF"${MAF}"."${trait}"."${win}".*."${test}".txt | xargs -I {}  tail -n+2 {} | awk '{print $(NF-2)}' | perl -npe 's/\n/,/' | sed 's/.$//'> tmp4
			ls Drop_chr"${chr}"_INTERVAL_QCv1.0_inverse_norm_naive.MAF"${MAF}"."${trait}"."${win}".*."${test}".txt | xargs -I {}  tail -n+2 {} | awk '{print $(NF-1)}' | perl -npe 's/\n/,/' | sed 's/.$//'> tmp5
			ls Drop_chr"${chr}"_INTERVAL_QCv1.0_inverse_norm_naive.MAF"${MAF}"."${trait}"."${win}".*."${test}".txt | xargs -I {}  tail -n+2 {} | awk '{print $NF}' | perl -npe 's/\n/,/' | sed 's/.$//'> tmp6
			paste tmp* >> Drop2m_chr"${chr}"_INTERVAL_QCv1.0_inverse_norm_naive.MAF"${MAF}"."${trait}"."${win}"."${test}".txt
		
        done
done

