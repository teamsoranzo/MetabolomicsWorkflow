#!/bin/bash

PHENO=/lustre/scratch115/projects/int_wes_metabol/Metabolon/CleanedData_160414/Pheno
VCFInput=/lustre/scratch115/projects/int_wes_metabol/WES_QCeddata_37/INTERVAL_QCv1.0.AC-AN_SelectIndsNoREl.vcf.gz
carrierF=/lustre/scratch115/projects/int_wes_metabol/carriers/Significant
listofmet=/lustre/scratch115/projects/int_wes_metabol/carriers/Significant/listofmetaboliteforBoxplot.txt
listofmet1=/lustre/scratch115/projects/int_wes_metabol/carriers/Significant/listofmetaboliteforBoxplot_ln.txt
listofvar=/lustre/scratch115/projects/int_wes_metabol/carriers/Significant/listofvariantsforBoxplot.txt
CleanMET=/lustre/scratch115/projects/int_wes_metabol/Metabolon/CleanedData_160414

bcftool=/software/vertres/bin-external/bcftools 


#Vars


$bcftool view -Ov -R ${listofvar}  ${VCFInput} > ${carrierF}/geno.vcf
bgzip -f ${carrierF}/geno.vcf
/software/hgi/pkglocal/plink-1.90b3/bin/plink --threads 4  --memory 500   --vcf ${carrierF}/geno.vcf.gz  --biallelic-only strict --keep-allele-order  --recode --out ${carrierF}/genoout


# Metabolites
selT=`head -n1 ${CleanMET}/Standardised_residuals_inverse_norm_metabolon_NoREl.csv | perl -npe 's/,/\n/g' | grep -w -n -f ${listofmet} - | cut -d: -f 1 | perl -npe 's/\n/,/' | perl -npe 's/,$/\n/'`

sex=` head -n1 ${CleanMET}/Metabolon_cleaned_data_160414NoREl.csv |  perl -npe 's/,/\n/g' | grep -w -n "sex" | cut -d: -f 1`

cut -d, -f1,${sex} ${CleanMET}/Metabolon_cleaned_data_160414NoREl.csv | perl -npe 's/,/ /g' | sed 's/Male/1/;s/Female/2/' | sort  > ${carrierF}/Indssex.txt

cut -d, -f1,${selT} ${CleanMET}/Standardised_residuals_inverse_norm_metabolon_NoREl.csv | perl -npe 's/,/ /g' | head -n1 > ${carrierF}/headerT.txt

cut -d, -f1,${selT} ${CleanMET}/Standardised_residuals_inverse_norm_metabolon_NoREl.csv | perl -npe 's/,/ /g' | tail -n+2 | sort -t" " -k1,1  > ${carrierF}/traits.txt

cat ${carrierF}/headerT.txt ${carrierF}/traits.txt > ${carrierF}/traitsH.txt

join -1 1 -2 1 ${carrierF}/Indssex.txt ${PHENO}/listIndsKeepNoREl.txt | awk '{print $1,$1,"0 0",$2}' > ${carrierF}/1.tmp

join -1 1 -2 1 ${carrierF}/traits.txt ${PHENO}/listIndsKeepNoREl.txt | cut -d" " -f2- | paste -d" " ${carrierF}/1.tmp - > ${carrierF}/INTERVAL_QCv1.0.gped

perl -npe 's/ /\n/g' ${carrierF}/headerT.txt | tail -n+2 | awk '{print "T",$1}'  > ${carrierF}/INTERVAL_QCv1.0.dat


# Metabolites ln


selT=`head -n1 ${CleanMET}/Metabolon_cleaned_data_160414NoREl.csv | perl -npe 's/,/\n/g' | grep -w -n -f ${listofmet1} - | cut -d: -f 1 | perl -npe 's/\n/,/' | perl -npe 's/,$/\n/'`
sex=` head -n1 ${CleanMET}/Metabolon_cleaned_data_160414NoREl.csv |  perl -npe 's/,/\n/g' | grep -w -n "sex" | cut -d: -f 1`
cut -d, -f1,${sex} ${CleanMET}/Metabolon_cleaned_data_160414NoREl.csv | perl -npe 's/,/ /g' | sed 's/Male/1/;s/Female/2/' | sort  > ${carrierF}/Indssex.txt
cut -d, -f1,${selT} ${CleanMET}/Metabolon_cleaned_data_160414NoREl.csv | perl -npe 's/,/ /g' | head -n1 > ${carrierF}/headerT.txt
cut -d, -f1,${selT} ${CleanMET}/Metabolon_cleaned_data_160414NoREl.csv | perl -npe 's/,\./,0\./g;s/,-\./,-0\./g;/g;s/,,,,,/,NA,NA,NA,NA,/g;s/,,,,/,NA,NA,NA,/g;s/,,,/,NA,NA,/g;s/,,/,NA,/g;s/^,/NA,/;s/,$/,NA/' | perl -npe 's/,/ /g' | tail -n+2 | sort -t" " -k1,1  > ${carrierF}/traits.txt
 cat ${carrierF}/headerT.txt ${carrierF}/traits.txt > ${carrierF}/traitsH.txt
 join -1 1 -2 1 ${carrierF}/Indssex.txt ${PHENO}/listIndsKeepNoREl.txt | awk '{print $1,$1,"0 0",$2}' > ${carrierF}/1.tmp
join -1 1 -2 1 ${carrierF}/traits.txt ${PHENO}/listIndsKeepNoREl.txt | cut -d" " -f2- | paste -d" " ${carrierF}/1.tmp - > ${carrierF}/INTERVAL_QCv1_ln.0.gped
perl -npe 's/ /\n/g' ${carrierF}/headerT.txt | tail -n+2 | awk '{print "T",$1}'  > ${carrierF}/INTERVAL_QCv1_ln.0.dat










 

