import glob, sys, os, string, datetime, re

# Data directories
run_dir = "/proj/omics4tb2/Collaborations/Vulcan"
data_dir = "%s/data/181203_K00235_0152_BHWNJ7BBXX_CBASS84_RNASeq" %run_dir
genome_dir = "%s/reference_STAR_snp" %run_dir
genome_fasta = "%s/Spis.Smic.genes.nt.longest.fa" %genome_dir

# get STAR bam alignment files
star_bams = glob.glob('%s/*/results_STAR_Spis_Smic_snp/*_star_Aligned.sortedByCoord.out.bam' %(data_dir))
print(len(star_bams))

# Function to run bcftools
def bcftools_variants():
    print
    print( "\033[34m Running bcftools Variant Calling.. \033[0m")

    # join all bam files for the command line
    joined_bams = " ".join(star_bams)

    # bcftools command
    cmd1 = '/users/sturkars/bcftools/bin/bcftools mpileup -Ou -f %s %s | /users/sturkars/bcftools/bin/bcftools call -Ou -mv | /users/sturkars/bcftools/bin/bcftools filter -s LowQual -e "QUAL<30 || DP>100" > var.flt.vcf' %(genome_fasta,joined_bams)

    print(cmd1)
    os.system(cmd1)

# run bcftools
bcftools_variants()

# Samtools indexing of BAM files if they are not indexed
for bamfile in star_bams:
    # get sample names
    sample_name = bamfile.split('/')[-1].split('_star_Aligned.sortedByCoord.out.bam')[0]
    sample_folder = "/".join(bamfile.split('/')[0:9])

    # index bam file
    cmd_index = '/users/sturkars/samtools/bin/samtools index %s' %(bamfile)
    #print('Running samtools indexing %s ' %(cmd_index))
    
    # run indexing
    #os.system(cmd_index)