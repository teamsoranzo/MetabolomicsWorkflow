#!/bin/bash

#rareMetW=/lustre/scratch115/projects/int_wes_metabol/raremetal/AllMetabolites/raremetalworker
CleanMET=/lustre/scratch115/projects/int_wes_metabol/Metabolon/CleanedData_160414
PHENO=/lustre/scratch115/projects/int_wes_metabol/Metabolon/CleanedData_160414/Pheno
#rareMetW_S=/nfs/team151/software/raremetal_4.13.6/raremetal/bin/raremetalworker
BCFTOOL=/software/vertres/bin-external/bcftools
VCFInput=/lustre/scratch115/projects/int_wes_metabol/WES_QCeddata_37

### Chromosome from job array
#i=`expr $1 - 1`
#i=`expr 20 - 1`
#chroms=({1..22} "X")
#chr=${chroms[$i]}
#echo $chr
#chr=20

## Trait index from second argument
filename=${PHENO}/ListOflistMet
declare -a ListMetFiles
ListMetFiles=( `cat "$filename"`)
j=`expr $1 - 1`
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

selT=`head -n1 ${CleanMET}/Standardised_residuals_inverse_norm_metabolon_NoREl.csv | perl -npe 's/,/\n/g' | grep -w -n -f ${PHENO}/${ListMetFile} - | cut -d: -f 1 | perl -npe 's/\n/,/' | perl -npe 's/,$/\n/'`

sex=` head -n1 ${CleanMET}/Metabolon_cleaned_data_160414NoREl.csv |  perl -npe 's/,/\n/g' | grep -w -n "sex" | cut -d: -f 1`

cut -d, -f1,${sex} ${CleanMET}/Metabolon_cleaned_data_160414NoREl.csv | perl -npe 's/,/ /g' | sed 's/Male/1/;s/Female/2/' | sort  > ${ListMetFileF}/Indssex.txt

cut -d, -f1,${selT} ${CleanMET}/Standardised_residuals_inverse_norm_metabolon_NoREl.csv | perl -npe 's/,/ /g' | head -n1 > ${ListMetFileF}/headerT_${ListMetFile}.txt

cut -d, -f1,${selT} ${CleanMET}/Standardised_residuals_inverse_norm_metabolon_NoREl.csv | perl -npe 's/,/ /g' | perl -npe 's/NA/x/g' | tail -n+2 | sort -t" " -k1,1  > ${ListMetFileF}/traits_${ListMetFile}.txt

cat ${ListMetFileF}/headerT_${ListMetFile}.txt ${ListMetFileF}/traits_${ListMetFile}.txt > ${ListMetFileF}/traitsH_${ListMetFile}.txt

#sort -t" " -k1,1 ${ListMetFileF}/Chr${chr}_INTERVAL_QCv1.0_${ListMetFile}.ped | cut -d" " -f1,7-  > ${ListMetFileF}/Chr${chr}_INTERVAL_QCv1.0_sort_${ListMetFile}.ped

### The list of inds to keep has been produce in other analysis where all the individuals with no phenotype or genetic information have been discarded total of 3926, if you only macth the common inds between genetic data and phenotype you end up 3931.   
#zcat ${VCFInput}/Chr${chr}/Chr${chr}_INTERVAL_QCv1.0.AC-AN.vcf.gz | grep -m1  "^#CHRO" | perl -npe 's/\t/\n/g' | grep "^EGA"  > ${ListMetFileF}/IndsInWES

#join -1 1 -2 1 ${ListMetFileF}/traits_${ListMetFile}.txt ${ListMetFileF}/IndsInWES | awk '{print $1,$1,"0 0 0"}' > ${ListMetFileF}/1_${ListMetFile}.tmp
#awk '{print $1}'  ${ListMetFileF}/1_${ListMetFile}.tmp > ${ListMetFileF}/listIndsKeep.txt

#join -1 1 -2 1 ${ListMetFileF}/traits_${ListMetFile}.txt ${ListMetFileF}/IndsInWES | paste -d" " ${ListMetFileF}/1_${ListMetFile}.tmp - > ${ListMetFileF}/Chr${chr}_INTERVAL_QCv1.0_${ListMetFile}.gped

#${BCFTOOL}  view -Oz -S ${ListMetFileF}/listIndsKeep.txt ${VCFInput}/Chr${chr}/Chr${chr}_INTERVAL_QCv1.0.AC-AN.vcf.gz  > ${ListMetFileF}/Chr${chr}_INTERVAL_QCv1.0.AC-AN_SelectInds.vcf.gz

### Instead using the list of inds previously ectracted.
join -1 1 -2 1 ${ListMetFileF}/Indssex.txt ${PHENO}/listIndsKeepNoREl.txt | awk '{print $1,$1,"0 0",$2}' > ${ListMetFileF}/1_${ListMetFile}.tmp

join -1 1 -2 1 ${ListMetFileF}/traits_${ListMetFile}.txt ${PHENO}/listIndsKeepNoREl.txt | cut -d" " -f2- | paste -d" " ${ListMetFileF}/1_${ListMetFile}.tmp - > ${ListMetFileF}/INTERVAL_QCv1.0_${ListMetFile}.gped

#${BCFTOOL}  view -Oz -S ${PHENO}/listIndsKeep.txt ${VCFInput}/Chr${chr}/Chr${chr}_INTERVAL_QCv1.0.AC-AN.vcf.gz  > ${ListMetFileF}/Chr${chr}_INTERVAL_QCv1.0.AC-AN_SelectInds.vcf.gz

#tabix -p vcf -f ${ListMetFileF}/Chr${chr}_INTERVAL_QCv1.0.AC-AN_SelectInds.vcf.gz

perl -npe 's/ /\n/g' ${ListMetFileF}/headerT_${ListMetFile}.txt | tail -n+2 | awk '{print "T",$1}'  > ${ListMetFileF}/INTERVAL_QCv1.0_${ListMetFile}.dat

