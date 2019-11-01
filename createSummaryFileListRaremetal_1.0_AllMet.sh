#i=`expr $1 - 1`
#i=`expr 20 - 1`
#chroms=({1..22} "X")
#chr=${chroms[$i]}
#echo $chr
#chr=20
P=/lustre/scratch115/projects/int_wes_metabol/raremetal/AllMetabolites/raremetalworker_inverse_norm

declare -a arr=({1..22} "X")

for chr in "${arr[@]}"
do
	for group in `ls ${P}/Chr${chr}/*.score.txt.gz`
	do
		 echo $group
		 fnames=`echo $group | sed 's/.gz//'`
		 ls $group  > ${fnames}.list
	 done
	
	for group in `ls ${P}/Chr${chr}/*.cov.txt.gz`
	do
		  echo $group
		  fnamec=`echo $group | sed 's/.gz//'`
		  ls $group  > ${fnamec}.list
	 done

done
