# RNASeq_Analysis
Analysis python codes for running read processing, trimming and quality control are located in RNASeq_Analysis_Python directory amd differential expression analysis related R scripts are located in RNASeq_Analysis_R directory
### Read processing, trimming and quality control
Paired‐end Illumina reads were checked for technical artifacts using TrimGalore version 0.4.3 (Krueger 2012) following Illumina default quality filtering steps. Reads were further trimmed for low-quality ends (Phred score < 20) and cleaned up for adapter contamination with TrimGalore. 
### Genomic reference files
For sequence alignment, reference genomic gene sets for Stylophora pistillata (v1.0) (Voolstra et al. 2017) and Symbiodinium microadriaticum (v1.0) (Aranda et al. 2016) were obtained from reefgenomics.org (Liew et al. 2016) and combined to create a merged reference. 
### Transcript abundance estimation
Transcript abundance estimation was performed by using kallisto v0.44.0 (Bray et al., PMID: 27043002). 

### Abundance estimate import and differential expression analysis
Differential gene expression analysis was performed using DESeq2 package v1.22.2 (Love et al. 2014) in R after importing kallisto transcript abundance estimates with tximport package v1.14.2 (Soneson et al., PMID:26925227). To account for large dispersion with low read counts and create more accurate log2 fold change (LFC) estimates, lfcShrink function for shrinking LFC estimates was applied. Transcripts with s-values (Stephens 2017) smaller than 0.005 were defined as significantly differentially expressed. In addition, transcript quantifications in transcripts per million (TPM) were obtained using Kallisto v0.44.0 (Bray et al. 2016) 

This is a [workflowr][] project to provide reproducible search results.


[workflowr]: https://github.com/jdblischak/workflowr
