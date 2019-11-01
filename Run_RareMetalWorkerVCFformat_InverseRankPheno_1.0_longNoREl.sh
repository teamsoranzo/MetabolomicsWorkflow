#!/bin/bash

rareMetW=/lustre/scratch115/projects/int_wes_metabol/raremetal/AllMetabolites/raremetalworker_inverse_norm
CleanMET=/lustre/scratch115/projects/int_wes_metabol/Metabolon/CleanedData_160414
PHENO=/lustre/scratch115/projects/int_wes_metabol/Metabolon/CleanedData_160414/Pheno
#rareMetW_S=/nfs/team151/software/raremetal_4.13.6/raremetal/bin/raremetalworker
rareMetW_S=/nfs/team151/software/raremetal_v4.14/Raremetal-master/bin/raremetalworker
BCFTOOL=/software/vertres/bin-external/bcftools
VCFInput=/lustre/scratch115/projects/int_wes_metabol/WES_QCeddata_37

### Chromosome from job array
i=`expr $1 - 1`
#i=`expr 20 - 1`
chroms=({1..22} "X")
chr=${chroms[$i]}
echo $chr
#chr=20

## Trait index from second argument
filename=${PHENO}/ListOflistMet
declare -a ListMetFiles
ListMetFiles=( `cat "$filename"`)
j=`expr $2 - 1`
ListMetFile=${ListMetFiles[$j]}
#ListMetFile=listMet1
#ListMetFile=listMet10
#ListMetFile=listMet11
#ListMetFile=listMet6
#ListMetFile=listMet13
#ListMetFile=listMet19
echo $ListMetFile

# To be change accordingly FOLDER
InputFile=/lustre/scratch115/projects/int_wes_metabol/raremetal/AllMetabolites/data_inverse_norm
ListMetFileF=${InputFile}/${ListMetFile}


mkdir -p ${InputFile}
mkdir -p ${ListMetFileF}


mkdir -p ${rareMetW}
mkdir -p ${rareMetW}/Chr${chr}
mkdir -p ${rareMetW}/Chr${chr}/${ListMetFile}
cd ${rareMetW}/Chr${chr}/${ListMetFile}

${rareMetW_S}  --ped ${ListMetFileF}/INTERVAL_QCv1.0_${ListMetFile}.gped  --dat ${ListMetFileF}/INTERVAL_QCv1.0_${ListMetFile}.dat --vcf ${VCFInput}/Chr${chr}/Chr${chr}_INTERVAL_QCv1.0.AC-AN_SelectIndsNoREl.vcf.gz  --prefix Chr${chr}_INTERVAL_QCv1.0_inverse_norm



#--traitName M20675

#awk 'BEGIN { FS=" "; } { if (NR==FNR) { pkeynline[$1]=$0; } else if (pkeynline[$1]) { print pkeynline[$1]; } else { print; } }' ${ListMetFileF}/traits.txt ${ListMetFileF}/Chr${chr}_INTERVAL_QCv1.0_naive.ped
#bsub -M10000 -R'rusage[mem=10000] select[mem>10000]' -o /nfs/users/nfs_l/lb17/My_Script_2016/int_wes_metabol_Scripts/LogJobs/TestCreateInputforRareMetalWorkerVCFformat.out  -q normal -- /nfs/team151/software/raremetal_4.13.6/raremetal/bin/raremetalworker --ped /lustre/scratch115/projects/int_wes_metabol/raremetalTEST/data/listMet1/TestwithVCF/Chr17_INTERVAL_QCv1.0_ModList4Met.gped --dat /lustre/scratch115/projects/int_wes_metabol/raremetalTEST/data/listMet1/TestwithVCF/Chr17_INTERVAL_QCv1.0_ModList4Met.dat --vcf /lustre/scratch115/projects/int_wes_metabol/raremetalTEST/data/listMet1/TestwithVCF/Chr17_INTERVAL_QCv1.0.SelectIndsUpdateINFO.vcf.gz  --prefix /lustre/scratch115/projects/int_wes_metabol/raremetalTEST/data/listMet1/TestwithVCF/Chr17_INTERVAL_QCv1.0_ModList4Met
####
#Add command to be rearrenge
#/software/vertres/bin-external/bcftools  view -Oz -S listIndsKeep.txt Chr17_INTERVAL_QCv1.0.recode.vcf > Chr17_INTERVAL_QCv1.0.SelectInds.vcf.gz
#bgzip Chr17_INTERVAL_QCv1.0.SelectInds.vcf
#tabix -p vcf -f Chr17_INTERVAL_QCv1.0.SelectInds.vcf.gz

