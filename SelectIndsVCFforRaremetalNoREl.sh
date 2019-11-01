#!/bin/bash

#rareMetW=/lustre/scratch115/projects/int_wes_metabol/raremetal/raremetalworker_inverse_norm
CleanMET=/lustre/scratch115/projects/int_wes_metabol/Metabolon/CleanedData_160414
PHENO=/lustre/scratch115/projects/int_wes_metabol/Metabolon/CleanedData_160414/Pheno
#rareMetW_S=/nfs/team151/software/raremetal_4.13.6/raremetal/bin/raremetalworker
BCFTOOL=/software/vertres/bin-external/bcftools
VCFInput=/lustre/scratch115/projects/int_wes_metabol/WES_QCeddata_37

### Chromosome from job array
i=`expr $1 - 1`
#i=`expr 20 - 1`
chroms=({1..22} "X")
chr=${chroms[$i]}
echo $chr
#chr=20
${BCFTOOL}  view -Oz -S ${PHENO}/listIndsKeepNoREl.txt ${VCFInput}/Chr${chr}/Chr${chr}_INTERVAL_QCv1.0.AC-AN.vcf.gz  > ${VCFInput}/Chr${chr}/Chr${chr}_INTERVAL_QCv1.0.AC-AN_SelectIndsNoREl.vcf.gz

/usr/bin/tabix -p vcf -f ${VCFInput}/Chr${chr}/Chr${chr}_INTERVAL_QCv1.0.AC-AN_SelectIndsNoREl.vcf.gz


