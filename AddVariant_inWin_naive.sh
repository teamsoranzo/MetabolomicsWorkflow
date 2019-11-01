#!/bin/bash


PHENO=/lustre/scratch115/projects/int_wes_metabol/Metabolon/CleanedData_160414/Pheno
## Change into working directory
rareMet_S=/nfs/team151/software/raremetal_v4.14/Raremetal-master/bin/raremetal
# Change the name of teh folder based on the MAF
raremetal=/lustre/scratch115/projects/int_wes_metabol/raremetal/AllMetabolites/raremetal_inverse_norm_AllMet/naive/AddVariantInWin
rareMetW=/lustre/scratch115/projects/int_wes_metabol/raremetal/AllMetabolites/raremetalworker_inverse_norm
signres=/lustre/scratch115/projects/int_wes_metabol/raremetal/AllMetabolites/raremetal_inverse_norm_AllMet/raremetal_sign_all_var_DeltaNN.txt


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
#	awk -F"\t"  -v w=$win '$1==w'  /lustre/scratch115/projects/int_wes_metabol/VEP_85_GENCODE24_GFF3/Chr${chr}/naive/*${MAF}.list > group${win}.list
#        group=group${win}.list        
#	echo $group
        score=${rareMetW}/Chr${chr}/Chr${chr}_INTERVAL_QCv1.0_inverse_norm.${trait}.singlevar.score.txt.list
        echo $score
        cov=${rareMetW}/Chr${chr}/Chr${chr}_INTERVAL_QCv1.0_inverse_norm.${trait}.singlevar.cov.txt.list
        echo $cov
	#ipval=`awk -F"\t" -v m=$trait -v w=$win -v t=$test '$25=="naive" && $19==m && $6==w && $20==t {print $18}' ${signres} | sort | uniq `
	#echo $ipval
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
	awk -F"\t" -v m=$trait -v w=$win -v t=$test '$25=="naive" && $19==m && $6==w && $20==t {print $0}' ${signres} | perl -npe 's/ /_/g'  | sort -r  -g -k26,26  |  awk '{print $1}' | perl -npe 's/\n/\t/g' | xargs -I {}  echo -e "$win\t"{} > group_all_sort.${win}.${test}.${trait}.txt 
	group_all_sort=group_all_sort.${win}.${test}.${trait}.txt
	len=`awk '{print NF}' $group_all_sort`
		for var in $(seq 2 $len) 
		do 
		echo $var
		cut -f-$var  $group_all_sort > group${win}.${test}.${trait}.${var}.list
		group_var=group${win}.${test}.${trait}.${var}.list
		#addonevar=`cut --complement -f-$var $group_all_sort | perl -npe 's/\t/;/g; s/;$//g' | cut -d";"  -f2-`
		test_vars=`cut -f-$var  $group_all_sort | cut -f2- | perl -npe 's/\t/;/g'`

				# run RAREMETAL
#                        	${rareMet_S} --summaryFiles $score --covFiles $cov --groupFile $group_var --maf ${MAF} --SKAT --burden --MB --VT  --longOutput --hwe 1e-06 --tabulateHits --hitsCutoff 1e-04  --prefix Chr${chr}_INTERVAL_QCv1.0_inverse_norm_naive.MAF${MAF}.${trait}.${win}.${var}.addMultiple
				
				if [ "${test}" == "SKAT_" ] || [ "${test}" == "VT_" ]
                		then
				
				ipval=`awk -F"\t" -v m=$trait -v w=$win -v t=$test '$25=="naive" && $19==m && $6==w && $20==t {print $18}' ${signres} | sort | uniq `				
				#test_vars=`awk '{print $3}' Chr${chr}_INTERVAL_QCv1.0_inverse_norm_naive.MAF${MAF}.${trait}.${win}.${var}.addMultiple.meta.${test}.results | tail -n1 `
				add_effectsize=`awk '{print $10}' Chr${chr}_INTERVAL_QCv1.0_inverse_norm_naive.MAF${MAF}.${trait}.${win}.${var}.addMultiple.meta.${test}.results | tail -n1`
				add_pval=`awk '{print $12}' Chr${chr}_INTERVAL_QCv1.0_inverse_norm_naive.MAF${MAF}.${trait}.${win}.${var}.addMultiple.meta.${test}.results | tail -n1`
				echo -e "window\tchr\tstart\tend\ttrait\tbioc\ttest\tapproach\tinitial_pval\tnumber_vars\ttest_vars\tadd_effectsize\tadd_pval" > AddMultiple_chr"${chr}"_INTERVAL_QCv1.0_inverse_norm_naive.MAF"${MAF}"."${trait}"."${win}"."${var}"."${test}".txt

				echo -e "$win\t$chr\t$sta\t$end\t$trait\t$bioc\t$test\tnaive\t$ipval\t$var\t$test_vars\t$add_effectsize\t$add_pval"  >> AddMultiple_chr"${chr}"_INTERVAL_QCv1.0_inverse_norm_naive.MAF"${MAF}"."${trait}"."${win}"."${var}"."${test}".txt

				elif [ "${test}" == "MB" ] || [ "${test}" == "burden" ]
				then
                		ipval=`awk -F"\t" -v m=$trait -v w=$win -v t=$test '$25=="naive" && $19==m && $6==w && $20==t {print $18}' ${signres} | sort | uniq `
                                #test_vars=`awk '{print $3}' Chr${chr}_INTERVAL_QCv1.0_inverse_norm_naive.MAF${MAF}.${trait}.${win}.${var}.addMultiple.meta.${test}.results | tail -n1 `
                                add_effectsize=`awk '{print $10}' Chr${chr}_INTERVAL_QCv1.0_inverse_norm_naive.MAF${MAF}.${trait}.${win}.${var}.addMultiple.meta.${test}.results | tail -n1`
                                add_pval=`awk '{print $11}' Chr${chr}_INTERVAL_QCv1.0_inverse_norm_naive.MAF${MAF}.${trait}.${win}.${var}.addMultiple.meta.${test}.results | tail -n1`
                                echo -e "window\tchr\tstart\tend\ttrait\tbioc\ttest\tapproach\tinitial_pval\tnumber_vars\ttest_vars\tadd_effectsize\tadd_pval" > AddMultiple_chr"${chr}"_INTERVAL_QCv1.0_inverse_norm_naive.MAF"${MAF}"."${trait}"."${win}"."${var}"."${test}".txt

                                echo -e "$win\t$chr\t$sta\t$end\t$trait\t$bioc\t$test\tnaive\t$ipval\t$var\t$test_vars\t$add_effectsize\t$add_pval"  >> AddMultiple_chr"${chr}"_INTERVAL_QCv1.0_inverse_norm_naive.MAF"${MAF}"."${trait}"."${win}"."${var}"."${test}".txt

				fi	
		done

		echo -e "window\tchr\tstart\tend\ttrait\tbioc\ttest\tapproach\tinitial_pval\tnumber_vars\ttest_vars_inWin\tadd_effectsizes\tadd_pvals" > Add2m_chr"${chr}"_INTERVAL_QCv1.0_inverse_norm_naive.MAF"${MAF}"."${trait}"."${win}"."${test}".txt
			ls AddMultiple_chr"${chr}"_INTERVAL_QCv1.0_inverse_norm_naive.MAF"${MAF}"."${trait}"."${win}".*."${test}".txt | sort -V |  xargs -I {}  tail -n+2 {} | tail -n1 | cut -f1-9 > tmp1
			ls AddMultiple_chr"${chr}"_INTERVAL_QCv1.0_inverse_norm_naive.MAF"${MAF}"."${trait}"."${win}".*."${test}".txt | sort -V | xargs -I {}  tail -n+2 {} | awk '{print $(NF-3)}' | perl -npe 's/\n/,/' | sed 's/.$//'> tmp3
			ls AddMultiple_chr"${chr}"_INTERVAL_QCv1.0_inverse_norm_naive.MAF"${MAF}"."${trait}"."${win}".*."${test}".txt | sort -V | xargs -I {}  tail -n+2 {} | awk '{print $(NF-2)}' | perl -npe 's/\n/,/' | sed 's/.$//'> tmp4
			ls AddMultiple_chr"${chr}"_INTERVAL_QCv1.0_inverse_norm_naive.MAF"${MAF}"."${trait}"."${win}".*."${test}".txt | sort -V | xargs -I {}  tail -n+2 {} | awk '{print $(NF-1)}' | perl -npe 's/\n/,/' | sed 's/.$//'> tmp5
			ls AddMultiple_chr"${chr}"_INTERVAL_QCv1.0_inverse_norm_naive.MAF"${MAF}"."${trait}"."${win}".*."${test}".txt | sort -V | xargs -I {}  tail -n+2 {} | awk '{print $NF}' | perl -npe 's/\n/,/' | sed 's/.$//'> tmp6
			paste tmp* >> Add2m_chr"${chr}"_INTERVAL_QCv1.0_inverse_norm_naive.MAF"${MAF}"."${trait}"."${win}"."${test}".txt
		
        done
done

