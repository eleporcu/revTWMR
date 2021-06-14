# revTWMR

##########Eleonora Porcu - Differentially expressed genes reflect disease-induced rather than disease-causing changes in the transcriptome


revTWMR.R requires two input files:
- a matrix containing the univariate effect size of n SNPs on k gene expressions (these estimates come from an eQTL study) and the univariate effect sizes on the phenotype. The last three columns are: BETA_GWAS SE N
- a two columns file (delimited by space) with gene and samplesize

This folder includes the matrixes used to perform revTWMR analyses for BMI using the trans-eQTLs results downloaded from https://molgenis26.gcc.rug.nl/downloads/eqtlgen/trans-eqtl/2018-09-04-trans-eQTLsFDR-CohortInfoRemoved-BonferroniAdded.txt.gz (Vosa et al, bioRxiv 2018) and the BMI-GWAS results downloaded from http://www.nealelab.is/uk-biobank/

How to run the script:

R < revTWMR.R --no-save bmi.matrix.betaGWAS bmi genes.N 


The output files is a .txt file containing the following columns:
-gene: name of the gene tested
-alpha_original: causal effect estimated by revTWMR before applying heterogeneity test
-SE_original: standard error of alpha_original
-P_original: Pvalue calculated from alpha_original and SE_original
-N_original: number of SNPs used as instruments
-Phet_original: Pvalue of heterogenity test
-alpha: causal effect estimated by revTWMR after removing the SNPs detected as outliers by the heterogenity test
-SE: standard error of alpha
-P: Pvalue calculated from alpha and SE
-Nend: number of SNPs left after outlier removal
-Phet: Pvalue of heterogenity test after outlier removal
