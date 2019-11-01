path=/lustre/scratch115/projects/int_wes_metabol/raremetal/AllMetabolites/raremetalworker_inverse_norm
PHENO=/lustre/scratch115/projects/int_wes_metabol/Metabolon/CleanedData_160414/Pheno


filename=${PHENO}/AllListMet.txt
declare -a traits
traits=( `cat "$filename"`)
j=`expr $1 - 1`
trait=${traits[$j]}
#trait=M01589
echo $trait

for i in {1..22} "X"
do
#i=1
	/usr/bin/bgzip ${path}/Chr${i}/Chr${i}_INTERVAL_QCv1.0_inverse_norm.${trait}.singlevar.score.txt
	/usr/bin/tabix -c "#" -s 1 -b 2 -e 2 ${path}/Chr${i}/Chr${i}_INTERVAL_QCv1.0_inverse_norm.${trait}.singlevar.score.txt.gz
	/usr/bin/bgzip ${path}/Chr${i}/Chr${i}_INTERVAL_QCv1.0_inverse_norm.${trait}.singlevar.cov.txt
	/usr/bin/tabix -c "#" -s 1 -b 2 -e 2 ${path}/Chr${i}/Chr${i}_INTERVAL_QCv1.0_inverse_norm.${trait}.singlevar.cov.txt.gz
	echo "Chr${i} Done for Met${trait}"
done



