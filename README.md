# MetabolomicsWorkflow
Stores all scripts and an Analysis workflow description for analayzing metabolomics WES in INTERVAL.

#### Analysis workflow on large scale metabolomics
###  Based on INTERVAL WES - METABOLON data
## Lorenzo Bomba 2019
#################

path for the variant group file generated using Klaudia Algorithm:
/lustre/scratch115/projects/int_wes_metabol/VEP_85_GENCODE24_GFF3/Chr*/All/ | naive/ | functional/ | LoF/
Example file:
Chr22variantsForGroupfile_naive0.000254711.txt 
#MAF=0.000254711

These are the MAF used in generating the windows for testing:
MAF<-c(1,0.05,0.01,0.005,0.001,0.000509423,0.000382067,0.000254711)
We set MAF threshold ≤ 0.001 for RVT 

Different variant selection strategies:

naive <- Exonic variants
functional <- Exonic variants with High or moderate impact annot from VEP.
LoF <- Exonic vaariants with LoF HC from Loftee VEP plugin
see VEPannot_selectVarforApproachMAF.R 


N.B.: Monomorphic sites removed.

## Run in All metabolites 17Feb2017

N.B.: this analysis is conducted in indivisulas that are not related 
(i.e. excluding any relatedness above third degree; Pi_HAT < 0.125). 
In total 3924 individuals were analysed.
All the scripts have the corresponding BSUB_* script used to submitted parallel jobs unless specify. 

##################
### Analysis #####
##################

#### Preparing files for raremetalworker

1) SelectIndsVCFforRaremetalNoREl.sh # selecting individuals from vcf files to use as input for raremetalworker

2) CreateInputforRareMetalWorkerVCFformat_InverseRankPheno_NoREl_1.0.sh # creating the gped and map file as input for raremetalworker

### Raremetalworker analysis

3) Run_RareMetalWorkerVCFformat_InverseRankPheno_1.0_BasementNoREl.sh # Run raremetalworker on the bigger chromosome see list form the BSUB_*

4) Run_RareMetalWorkerVCFformat_InverseRankPheno_1.0_longNoREl.sh # Run raremetalworker on the smaller chromosome see list form the BSUB_*

5) CheckRaremetalworkerFiles.sh #Check if the file are succesfully created - in checkfilelinesAll.txt

6) MoveRaremetalWorkerRes_1.0_AllMet.sh #Move the Score and Cov files from raremetalworker up to a parent directory - to remove the listMet* directory in subsequent processes

7) createSummaryFileListRaremetal_1.0_AllMet.sh

### Raremetal analysis

8a) BSUB_run_raremetal_INTERVAL_inverse_norm_naive_AllMet_1.0.sh
run_raremetal_INTERVAL_inverse_norm_naive_AllMet_1.0.sh

8b) BSUB_run_raremetal_INTERVAL_inverse_norm_functional_AllMet_1.0.sh
run_raremetal_INTERVAL_inverse_norm_functional_AllMet_1.0.sh

8c) BSUB_run_raremetal_INTERVAL_inverse_norm_LoF_AllMet_1.0.sh
run_raremetal_INTERVAL_inverse_norm_LoF_AllMet_1.0.sh

9) MoveFailed.sh - move files in other directory because failed due to not enough samples with phenoytpe.

### Merge and annotate results files from raremetal and raremetalworkers

10a) BSUB_MergeResChrAllraremetal_AllMet_invers_norm_XApproach_GeneBased_1.0.sh # Merge the resulting splitted Chrmomosomes in Genebased meth
MergeResChrAllraremetal_AllMet_invers_norm_XApproach_GeneBased_1.0.sh

10b) BSUB_MergeResChrAllraremetal_AllMet_invers_norm_XApproach_SingleVar_1.0.sh # Merge the resulting splitted Chrmomosomes in single var
MergeResChrAllraremetal_AllMet_invers_norm_XApproach_SingleVar_1.0.sh

### If some raremetal jobs fail this is a script to rerun selecting metabolites and chromosomes that failed:
RESUB_Exit_BSUB_run_raremetal_INTERVAL_inverse_norm_naive_AllMet_1.0.sh
RESUB_Exit_BSUB_run_raremetal_INTERVAL_inverse_norm_functional_AllMet_1.0.sh
RESUB_Exit_BSUB_run_raremetal_INTERVAL_inverse_norm_LoF_AllMet_1.0.sh


### Accessory scripts to create statistics and follow-up analysis on RVT results

11)BSUB_CalculateNumVarandWinTestedRaremetal_invers_norm.sh # simple statistic on number of variants and number of windows analysed
CalculateNumVarandWinTestedRaremetal_invers_norm.sh

12a)QQplot_AllMet.R # Calculate lambda and produce QQplots of gene-based (10 metabolites in each)
BSUB_QQplot_AllMet.sh
12b)QQplot_AllMet_SingleVar_firsthalfMet.R - QQplot_AllMet_SingleVar_secondhalfMet.R # Calculate lambda and produce QQplots of SingleVar (10 metabolites in each)
BSUB_QQplot_AllMet_SingleVar.sh

13a)BSUB_lambdaMerge.sh # Merge all the lambda and plot bowplot X approach and methods
lambdaMerge.sh
13b)BSUB_BoxPlotLambda.R # Merge all the lambda and plot bowplot X approach and methods
BoxPlotLambda.R

14)MatchGenesGWAS.R # match the Gene and Gene windows for single var hits and filter for min pval within a gene -- this has been used for another set of significant hits where no MAF filtering was done i.e. including rare variants that could be within or outside the gene (both cases are included in the final table -- see In SingleVar folder for ALL_Metabolites_*).
BSUB_MatchGenesGWAS.sh

15)QQplot_Trial_raremetalworker.R # Sanity check script for strange behaviour of SingleVAR QQPlot when using all the variant. Some Singleton form the same sample having the same identical pvalues in more than 100 variants. Seems to be fine. Replot the QQplot singleVar without the rare variants (i.e.  MAF > 0.001 ) 

16)MergeTableResults.R # Merge tables of result based on Variant level or Gene level -- Not yet unique variants or unique genes.
BSUB_MergeTableResults.sh

17)BSUB_GzipIndexRaremetalworkerWESresults.sh
GzipIndexRaremetalworkerWESresults.sh


# retrieve annotaion for gene and variants using VEP plugins and as annotation file a GFF3 of GENCODE version 24
18)Vep_phenotype.sh
19)Vep_phenotype_extra.sh

20)Select_CutoffPval_MergeMetabolites_GeneBased_AllHits.sh #This script is just combining all the metabolites without any pvalue threshold in place (all the hits included)

# ExtractPvalAllHits.sh #Extract from the output of the  preovius sh only the Pvalue and a unique ID win-met-test-approach combinig all the approaches
21)ExtractPvalFDRapp.sh
22)ExtractPvalFDRapp.sh
23)ExtractPvalFDRapptest.sh


# these scripts calclute the FDR (padj BH) using all the pvals, pvals divided by approach and pvals dividide by approach and test
24)calculateFDR.R 
25)calculateFDRapp.R
26)calculateFDRapptest.R 

# Drop one variant at the time and ricalculate the RVT for each approach
27)DropVariant_inWin_LoF.sh
28)DropVariant_inWin_functional.sh
29)DropVariant_inWin_naive.sh

# Add one variant at the time based on the magnitude of the deviation of dropone pvalue RVT and initial full model pvalue RVT  and recalculate the RVT each time and for each approach 
30)AddVariant_inWin_naive.sh
31)AddVariant_inWin_LoF.sh
32)AddVariant_inWin_functional.sh


33)extract the genotype and phenotype information to plot boxplot of carriers vs non-carriers
extractGenotypeAndMetabolite4Boxplot.sh
