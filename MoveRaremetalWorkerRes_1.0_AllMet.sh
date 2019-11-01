
pt=/lustre/scratch115/projects/int_wes_metabol/raremetal/AllMetabolites/raremetalworker_inverse_norm/
cd $pt
declare -a arr=({1..22} "X")

for i in "${arr[@]}"
do
   chr=Chr${i}
   echo $chr
   cd $chr
   mv listMet*/${chr}_INTERVAL_QCv1.0*.singlevar.score.txt ${pt}/${chr}/ 
   mv listMet*/${chr}_INTERVAL_QCv1.0*.singlevar.cov.txt ${pt}/${chr}/
   #mkdir -p listMet1
   #rm -f ${chr}_INTERVAL_QCv1* 
   #rename  -v 's/_inverse_norm//' Chr*
   #cd listMet3Bad3GoodComp/
   #rename  -v 's/_inverse_norm//' Chr*  
   cd ..
# or do whatever with individual element of the array
done

