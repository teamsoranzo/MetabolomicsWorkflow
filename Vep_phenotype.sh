#Path VARIABLE
GFF=/nfs/team151/software/VEP_NewTrialVersion_84_GRCh37/GFF/gencode.v24lift37.annotation.gff3.gz
FASTA=/lustre/scratch114/resources/ref/Homo_sapiens/1000Genomes_hs37d5/hs37d5.fa
ORIG_GEN_DATA=/lustre/scratch115/projects/int_wes_metabol/WES_QCeddata_37
VEP=/lustre/scratch115/projects/int_wes_metabol/VEP_85_GENCODE24_GFF3
plinkSKAT_functional=/lustre/scratch115/projects/int_wes_metabol/plinkSKAT/functional
GENCODE=/lustre/scratch115/projects/int_wes_metabol/GENCODE
#Path software
PLINK=/software/hgi/pkglocal/plink-1.90b3/bin/plink
VCFTOOLS=/software/hgi/pkglocal/vcftools-0.1.11/bin/vcftools
PYTHON=/software/bin/python
VEPSoft=/nfs/users/nfs_l/lb17/Softwares/variant_effect_predictor_NewTrialVersion/ensembl-vep/vep.pl
#MyChr=Chr22

#zcat /lustre/scratch115/projects/int_wes_metabol/WES_QCeddata_37/Chr22/Chr22_INTERVAL_QCv1.0.AC-AN_SelectIndsNoREl.vcf.gz |  cut -f1-8 | head -400 > ${VEP}/${MyChr}/${MyChr}_TrialINPUTPhenotype

#mkdir -p ${VEP}/Phenotype_vep_res/
#zcat /lustre/scratch115/projects/int_wes_metabol/WES_QCeddata_37/INTERVAL_QCv1.0.AC-AN.vcf.gz  |  cut -f1-8 > ${VEP}/Phenotype_vep_res/INTERVAL_QCv1.0.AC-AN_nogenotypes.vcf

#${VCFTOOLS} --vcf /lustre/scratch115/projects/int_wes_metabol/VEP_85_GENCODE24_GFF3/Phenotype_vep_res/INTERVAL_QCv1.0.AC-AN_nogenotypes.vcf  --bed /lustre/scratch115/projects/int_wes_metabol/VEP_85_GENCODE24_GFF3/Phenotype_vep_res/raremetal_sign_all_var_U.bed --out /lustre/scratch115/projects/int_wes_metabol/VEP_85_GENCODE24_GFF3/Phenotype_vep_res/INTERVAL_QCv1.0.AC-AN_nogenotypes_SelectedVars.vcf --recode --keep-INFO-all


perl $VEPSoft --force_overwrite  -gff  ${GFF} --fasta ${FASTA}  -i ${VEP}/Phenotype_vep_res/INTERVAL_QCv1.0.AC-AN_nogenotypes_SelectedVars.vcf.recode.vcf --allele_number  --plugin Phenotypes,file=/nfs/users/nfs_l/lb17/.vep/Plugins/Phenotypes.pm_homo_sapiens_84_GRCh37.bed.gz --pick  --o ${VEP}/Phenotype_vep_res/PhenotypePickUsingGFF3_GENCODE24_lift37.txt

perl $VEPSoft --force_overwrite  -gff  ${GFF} --fasta ${FASTA}  -i ${VEP}/Phenotype_vep_res/INTERVAL_QCv1.0.AC-AN_nogenotypes_SelectedVars.vcf.recode.vcf --allele_number  --plugin Phenotypes,file=/nfs/users/nfs_l/lb17/.vep/Plugins/Phenotypes.pm_homo_sapiens_84_GRCh37.bed.gz --pick  --vcf  --o ${VEP}/Phenotype_vep_res/PhenotypePickUsingGFF3_GENCODE24_lift37.vcf

perl $VEPSoft --force_overwrite  -gff  ${GFF} --fasta ${FASTA}  -i ${VEP}/Phenotype_vep_res/INTERVAL_QCv1.0.AC-AN_nogenotypes_SelectedVars.vcf.recode.vcf --allele_number  --plugin Phenotypes,file=/nfs/users/nfs_l/lb17/.vep/Plugins/Phenotypes.pm_homo_sapiens_84_GRCh37.bed.gz --pick --JSON --o ${VEP}/Phenotype_vep_res/PhenotypePickUsingGFF3_GENCODE24_lift37.json

##### LAST RUN 

perl $VEPSoft --force_overwrite  -gff  ${GFF} --fasta ${FASTA}  -i ${VEP}/Phenotype_vep_res/INTERVAL_QCv1.0.AC-AN_nogenotypes_SelectedVars_all.vcf.recode.vcf  --allele_number  --plugin Phenotypes,file=/nfs/users/nfs_l/lb17/.vep/Plugins/Phenotypes.pm_homo_sapiens_84_GRCh37.bed.gz --plugin PolyPhen_SIFT,dir=/nfs/team151/lb17t/ --plugin SameCodon --plugin LoF,human_ancestor_fa:/nfs/users/nfs_l/lb17/LoFtee/HumAncestor/human_ancestor.fa.gz,filter_position:0.05,phylocsf.sql  --canonical  --pick  --vcf --o ${VEP}/Phenotype_vep_res/PhenotypePickUsingGFF3_GENCODE24_lift37_all.vcf


-plugin dbNSFP,/lustre/scratch115/projects/int_wes_metabol/VEP_85_GENCODE24_GFF3/Phenotype_vep_res/dbNSFP.gz,MetaSVM_score,MetaSVM_rankscore,MetaSVM_pred,Reliability_index,CADD_raw,CADD_raw_rankscore,CADD_phred,phyloP46way_primate,phyloP46way_primate_rankscore,clinvar_rs,clinvar_clnsig,clinvar_trait,clinvar_golden_stars

#   --af --af_1kg --af_esp --af_exac --appris --biotype --canonical --check_existing --domains  --polyphen b --pubmed --regulatory --sift b --species homo_sapiens --symbol --tsl

#perl vep.pl --af --af_1kg --af_esp --af_exac --appris --biotype --canonical --check_existing --domains --gencode_basic --numbers --pick --plugin dbNSFP,[path_to]/ensweb-data[path_to]/dbNSFP2.9.2.txt.gz,MetaSVM_score,MetaSVM_rankscore,MetaSVM_pred,Reliability_index,CADD_raw,CADD_raw_rankscore,CADD_phred,phyloP46way_primate,phyloP46way_primate_rankscore,clinvar_rs,clinvar_clnsig,clinvar_trait,clinvar_golden_stars --plugin Condel,[path_to]/ensweb-data[path_to]/config,b --plugin Blosum62 --plugin LoFtool,[path_to]/ensweb-data[path_to]/LoFtool_scores.txt --polyphen b --pubmed --regulatory --sift b --species homo_sapiens --symbol --tsl --cache --input_file [input_data]


#Phenotypes,include_types=Variation 

#,include_sources:AMDGC & ClinVar & dbGaP & GEFOS & GIANT & HGMD-PUBLIC & MAGIC & OMIM & NHGRI-EBIGWAScatalog & Teslovich & uniprot

#present sources: ClinVar,HGMD-PUBLIC,OMIM,Uniprot

