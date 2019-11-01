#!/bin/bash


## echo 'bash ~/scripts/uk10k/run_raremetal_lipids_107-regions_genome_uk10k.sh ${LSB_JOBINDEX}' | bsub -J "trait[1-4]" -o /lustre/scratch113/projects/uk10k/users/kw8/lsfout4/run_raremetal_lipids_107-regions_genome_uk10k_trait%I.out -R "select[mem>=1000] rusage[mem=1000]" -M1000 -q normal 




#################################### Run RAREMETAL #####################################
PHENO=/lustre/scratch115/projects/int_wes_metabol/Metabolon/CleanedData_160414/Pheno
## Change into working directory
rareMet_S=/nfs/team151/software/raremetal_v4.14/Raremetal-master/bin/raremetal
# Change the name of teh folder based on the MAF
raremetal=/lustre/scratch115/projects/int_wes_metabol/raremetal/AllMetabolites/raremetal_inverse_norm_AllMet/LoF
rareMetW=/lustre/scratch115/projects/int_wes_metabol/raremetal/AllMetabolites/raremetalworker_inverse_norm


mkdir -p $raremetal
cd $raremetal

# Chromosome from job array
#i=`expr $1 - 1`
#i=`expr 20 - 1`
#chroms=({1..22} "X")
#chr=${chroms[$i]}
#echo $chr

## Trait index from second argument
filename=${PHENO}/AllListMet.txt
declare -a traits
traits=( `cat "$filename"`)
j=`expr $1 - 1`
trait=${traits[$j]}
echo $trait

mkdir -p  $trait

cd $trait

## run RAREMETAL

declare -a arr=("0.001")

## now loop through the above array
for MAF in "${arr[@]}"
do
	echo $MAF
	for group in `ls /lustre/scratch115/projects/int_wes_metabol/VEP_85_GENCODE24_GFF3/Chr*/LoF/*${MAF}.list`
	do 
    		chr=`ls $group |  cut -d/ -f9 | cut -d_ -f1 | cut -b4-`
    		echo $chr
    		score=${rareMetW}/Chr${chr}/Chr${chr}_INTERVAL_QCv1.0_inverse_norm.${trait}.singlevar.score.txt.list
    		echo $score

    		cov=${rareMetW}/Chr${chr}/Chr${chr}_INTERVAL_QCv1.0_inverse_norm.${trait}.singlevar.cov.txt.list
    		echo $cov

    		echo "raremetal --summaryFiles $score --covFiles $cov --groupFile $group --maf 0.001 --SKAT --burden --longOutput --tabulateHits --prefix Chr${chr}_INTERVAL_QCv1.0_LoF.${trait}"
    		${rareMet_S} --summaryFiles $score --covFiles $cov --groupFile $group --maf ${MAF} --SKAT --burden --MB --VT  --longOutput --hwe 1e-06 --tabulateHits --hitsCutoff 1e-04  --prefix Chr${chr}_INTERVAL_QCv1.0_inverse_norm_LoF.MAF${MAF}.${trait}
	done
done
