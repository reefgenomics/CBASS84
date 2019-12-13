## same RNA-seq pipeline specific to using C5 reference genome fasta and gff
import glob, sys, os, string, datetime, re
now = datetime.datetime.now()
timeprint = now.strftime("%Y-%m-%d %H:%M")

# Input files
# data and results directories
run_dir = "/proj/omics4tb2/Collaborations/Vulcan"
data_dir = "%s/data/181203_K00235_0152_BHWNJ7BBXX_CBASS84_RNASeq" %run_dir
genome_dir = "%s/Genomics/Stylophora_Symbiodinium" %run_dir
genome_fasta = "%s/Spis_Smic.genome.scaffold.final.fa" %genome_dir
genome_gff = "%s/spis_smic.genome.annotation.gff3" %genome_dir
folder_name = str(sys.argv[1])
data_folder = "%s/%s" %(data_dir,folder_name)

# Programs
star_path = "/users/sturkars/STAR-2.7.0a/bin/Linux_x86_64/STAR" # path to STAR executable


####################### Create results directories ###############################
def create_dirs(data_trimmed_dir,fastqc_dir,results_dir,htseq_dir):
    dirs = [data_trimmed_dir,fastqc_dir,results_dir,htseq_dir]
    for dir in dirs:
        # create results folder
        #print(dir)
        if not os.path.exists('%s' %(dir)):
            os.makedirs('%s' %(dir))
        else:
            print('\033[31m %s directory exists. Not creating. \033[0m' %(dir))


####################### Trimgalore for quality and trimming ###############################
def trimgalore(first_pair_file,second_pair_file,folder_name,sample_id,file_ext,data_trimmed_dir,fastqc_dir):
    #print("1stpair:%s, 2ndpair:%s, folder_name:%s, sample_name:%s")%(first_pair_file,second_pair_file,folder_name,sample_name)
    print
    print ("\033[34m Running TrimGalore \033[0m")
    # create sample spepcific trimmed directory
    if not os.path.exists('%s' %(data_trimmed_dir)):
        os.makedirs('%s' %(data_trimmed_dir))
    # create sample spepcific fastqcdirectory
    if not os.path.exists('%s' %(fastqc_dir)):
        os.makedirs('%s' %(fastqc_dir))
    # run Command
    cmd = 'trim_galore --fastqc_args "--outdir %s/" --paired --output_dir %s/ %s %s' %(fastqc_dir,data_trimmed_dir,first_pair_file, second_pair_file)
    print
    print( '++++++ Trimgalore Command:', cmd)
    print
    #os.system(cmd)


####################### Collect trimmed data files ###############################
def collect_trimmed_data(data_trimmed_dir,file_ext):
    # define result files
    if file_ext == "gz":
        first_pair_trimmed = glob.glob('%s/*_val_1.fq.gz'%(data_trimmed_dir))
        second_pair_trimmed = glob.glob('%s/*_val_2.fq.gz'%(data_trimmed_dir))
    else:
        first_pair_trimmed = glob.glob('%s/*_val_1.fq'%(data_trimmed_dir))
        second_pair_trimmed = glob.glob('%s/*_val_2.fq'%(data_trimmed_dir))
    print('Trimmed Files:\n 1st:%s \n 2nd:%s' %(first_pair_trimmed,second_pair_trimmed))
    print
    first_pair_group = ' '.join(first_pair_trimmed)
    second_pair_group = ' '.join(second_pair_trimmed)
    pair_files = []
    for file in first_pair_trimmed:
        mate_file = file.replace('R1_001_val_1.fq','R2_001_val_2.fq')
        paired_mates = file + ' ' + mate_file
        pair_files.append(paired_mates)

    star_input_files = ' '.join(pair_files)

    return first_pair_group,second_pair_group, star_input_files


####################### Run STAR #####################################
def run_star(first_pair_group,second_pair_group,results_dir,star_input_files):
    print
    print('\033[33mRunning STAR! \033[0m')

    outfile_prefix = '%s/%s_star_' %(results_dir,folder_name)
    star_options ="--runThreadN 4 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --alignSJoverhangMin 5 --alignSJDBoverhangMin 3 --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.4 --outFilterScoreMinOverLread 0.4 --quantMode TranscriptomeSAM GeneCounts"

    cmd = '%s --genomeDir %s %s --readFilesIn %s %s --outFileNamePrefix %s' %(star_path, genome_dir, star_options, first_pair_group, second_pair_group, outfile_prefix)
    print('STAR run command:%s' %cmd)
    os.system(cmd)


####################### Run HTSEq Count ###############################
def run_htseq(htseq_dir,results_dir):
    print
    print('\033[33mRunning htseq-count! \033[0m')
    htseq_input = '%s/%s_star_Aligned.sortedByCoord.out.bam' %(results_dir,folder_name)
    cmd = 'htseq-count -s "reverse" -t "exon" -i "Parent" -r pos --max-reads-in-buffer 60000000 -f bam %s %s > %s/%s_htseqcounts.txt' %(htseq_input,genome_gff,htseq_dir,folder_name)
    print('htseq-count run command:%s' %cmd)
    os.system(cmd)

####################### Create STAR index ###############################
def genomeIndex():
    index_cmd = '%s --runMode genomeGenerate --runThreadN 4 --genomeDir %s --genomeFastaFiles %s --sjdbGTFfile %s --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 149 --genomeChrBinNbits 16 --limitSjdbInsertNsj 1200356' %(star_path, genome_dir, genome_fasta, genome_gff)
    print(index_cmd)

    print ("\033[34m %s Indexing genome... \033[0m") %(timeprint)
    if os.path.exists('%s/SAindex' %(genome_dir)):
        print ('Genome indexes exist. Not creating!')
    else:
        print ('Creating genome indexes')
    os.system(index_cmd)


####################### Running the Pipeline ###############################
def run_pipeline():
    folder_count = 1

    # Loop through each data folder
    #for data_folder in data_folders:
    folder_name = data_folder.split('/')[-1]
    print
    print
    print('\033[33mProcessing Folder: %s\033[0m' %(folder_name))

    # Get the list of first file names in paired end sequences
    first_pair_files = glob.glob('%s/*_R1*.fastq*' %(data_folder))
    #second_pair_files = glob.glob('%s/_R2*.fastq*' %(data_folder))

    # Program specific results directories
    data_trimmed_dir = "%s/%s/trimmed" %(data_dir,folder_name)
    fastqc_dir = "%s/%s/fastqc_results" %(data_dir,folder_name)
    results_dir = "%s/%s/results_STAR_Spis_Smic" %(data_dir,folder_name)
    htseq_dir = "%s/htseqcounts" %(results_dir)

    # Run create directories function to create directory structure
    create_dirs(data_trimmed_dir,fastqc_dir,results_dir,htseq_dir)


    # Loop through each file and create filenames
    file_count = 1
    for first_pair_file in first_pair_files:
        first_file_name_full = first_pair_file.split('/')[-1]

        second_pair_file = first_pair_file.replace('_R1', '_R2')
        second_file_name_full = second_pair_file.split('/')[-1]
        file_ext = first_pair_file.split('.')[-1]

        print ('\033[32m Processing File: %s of %s (%s)\033[0m' %(file_count, len(first_pair_files), first_file_name_full ))

        first_file_name = re.split('.fastq|.fastq.gz',first_file_name_full)[0]
        second_file_name = re.split('.fastq|.fastq.gz',second_file_name_full)[0]
        print('first_file_name:%s, second_file_name:%s' %(first_file_name,second_file_name))

        # Collect Sample attributes
        exp_name = folder_name
        print("exp_name: %s" %(exp_name))
        lane = first_file_name.split("_")[-1]
        print("Lane: %s" %(lane))
        sample_id = re.split('.fastq|.fastq.gz', first_file_name)[0]
        print("sample_id: %s" %(sample_id))
        #sys.exit()

        # Run TrimGalore
        trimgalore(first_pair_file,second_pair_file,folder_name,sample_id,file_ext,data_trimmed_dir,fastqc_dir)
        file_count = file_count + 1

        # Collect Trimmed data for input into STAR
        first_pair_group,second_pair_group,star_input_files = collect_trimmed_data(data_trimmed_dir,file_ext)

        # Run STAR
        run_star(first_pair_group,second_pair_group,results_dir,star_input_files)

        # Run HTSeq count
        run_htseq(htseq_dir,results_dir)

        folder_count = folder_count + 1
        #sys.exit()
    return data_trimmed_dir,fastqc_dir,results_dir

#genomeIndex()
data_trimmed_dir,fastqc_dir,results_dir = run_pipeline()
