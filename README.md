# GP_Sorghumhybrids
Mixed model analysis of phenotypes and heterosis in a line x tester hybrids data, heritability estimation, variance component analysis and genomic best linear unbiased prediction (GBLUP).

## data
Phenotypic data for hybrids is in "**HDP_Heterosis_ByBlock.csv**", contains hybrid values (Hyb_v), mid parent heterosi (MPH), etc.

**GRM_phased_thinned.RData** has the genomic relationship matrices (GRM) calculated from the marker data, G = inbred additive GRM, A = hybrid additive GRM, and D = hbyrid dominance matrix.

## results
Output folder has RData files that has details from model fit and variance components and BLUPs/GBLUPs.
Main results are in the parent directory of this subdirectory.

## scripts
Scripts used to run model fit, plotting, and qsub files for submission of jobs to the cluster.
