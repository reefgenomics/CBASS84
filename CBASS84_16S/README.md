#CBASS84_16S

This directory contains the scripts used to analyze and create figures for the 16S dataset.

## 16S workflow
1. The script `CBASS84_16S_QC_and_contaRemoval.R` reads Mothur output tables, identifies and removes contaminant OTUs. It also removes samples with less than 1000 counts and exports new OTU and taxonomy tables.
2. The script `CBASS84_16S_diversity.R` plots averaged and replicated barplots of the 20 most abundant bacterial families.
3. The script `CBASS84_16S_stats.R` runs PERMANOVAS to compare overall bacterial diversity across temperatures and sites
4. The script `CBASS84_16S_ordination.R` plots PCoAs of all samples
5. The script `CBASS84_16S_DESeq2.R` identifies differentially abundant OTUs between sites and temperatures
