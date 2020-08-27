import glob
import sys
import os
import re

# data and results directories
run_dir = "/proj/omics4tb2/Collaborations/Vulcan"
data_dir = "%s/data/181203_K00235_0152_BHWNJ7BBXX_CBASS84_RNASeq" %run_dir
genome_dir = "%s/Genomics/Stylophora_Symbiodinium" %run_dir
#transcriptome_file = "%s/GCF_002571385.1_Stylophora_pistillata_v1_rna.fna.gz" %genome_dir
transcriptome_file = "%s/Spis_Smic.genome.annotation.CDS.longest.fa.gz" %genome_dir
folder_name = str(sys.argv[1])
print(folder_name)
data_folder = "%s/%s" %(data_dir,folder_name)

############# Functions ##############

####################### Collect data files ###############################
def get_data():
    data_folders = glob.glob('%s/*' %(data_dir))
    data_folders = [element for element in data_folders if element not in ('%s/*/trimmed,%s/*/fastqc_results')%(data_dir,data_dir)]
    print('data_folders: %s, %s' %(len(data_folders),data_folders))
    return data_folders


####################### Create results directories ###############################
def create_dirs(data_trimmed_dir,fastqc_dir,results_dir): 
    dirs = [data_trimmed_dir,fastqc_dir,results_dir]
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
    
    kallisto_input_files = ' '.join(pair_files)
    
    return first_pair_group,second_pair_group, kallisto_input_files


####################### Run Salmon ###############################
def run_salmon(first_pair_group,second_pair_group,results_dir):
    print
    print('\033[33mRunning salmon! \033[0m')
    salmon_cmd = '/users/sturkars/salmon/bin/salmon quant -i %s/thaps_transcripts_index -l A -1 %s -2 %s -o %s --validateMappings --seqBias --numBootstraps 100 --gcBias' %(genome_dir, first_pair_group, second_pair_group, results_dir)
    print('Salmon run command:%s' %salmon_cmd)
    #os.system(salmon_cmd)


####################### Create Salmon index ###############################
def salmon_index():
    print
    print('\033[33mRunning salmon index! \033[0m')
    index_cmd = '/users/sturkars/salmon/bin/salmon index -t %s -i %s/thaps_transcripts_index --type quasi -k 31' %(transcriptome_file,genome_dir)
    print('salmon index command:%s' %(index_cmd))
    #os.system(index_cmd)
 
####################### Run Kalisto ###############################
def run_kallisto(first_pair_group,second_pair_group,results_dir,kallisto_input_files):
    print
    print('\033[33mRunning kallisto! \033[0m')
    kallisto_cmd = '/users/sturkars/kallisto/kallisto quant -i %s/Spis_Smic_transcripts_kallistoindex %s -o %s  -b 100 --bias -t 4 --fusion' %(genome_dir, kallisto_input_files, results_dir)
    print('Kallisto run command:%s' %kallisto_cmd)
    os.system(kallisto_cmd)
    
 ####################### Create Kallisto index ###############################
def kallisto_index():
    print
    print('\033[33mRunning kallisto index! \033[0m')
    kallistoindex_cmd = '/users/sturkars/kallisto/kallisto index -i %s/Spis_Smic_transcripts_kallistoindex %s' %(genome_dir,transcriptome_file)
    print('kallisto index command:%s' %(kallistoindex_cmd))
    os.system(kallistoindex_cmd)
   

####################### Running the Pipeline ###############################    
def run_pipeline():
    folder_count = 1
    #data_folders = get_data()

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
    results_dir = "%s/%s/results_kallisto_Spis_Smic" %(data_dir,folder_name)
    
    # Run create directories function to create directory structure
    create_dirs(data_trimmed_dir,fastqc_dir,results_dir)
 

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
        
        # 01. Run TrimGalore
        trimgalore(first_pair_file,second_pair_file,folder_name,sample_id,file_ext,data_trimmed_dir,fastqc_dir)

        file_count = file_count + 1
        
        # Run folder level salmon analysis
        first_pair_group,second_pair_group,kallisto_input_files = collect_trimmed_data(data_trimmed_dir,file_ext)
        #run_salmon(first_pair_group,second_pair_group,results_dir)
        run_kallisto(first_pair_group,second_pair_group,results_dir,kallisto_input_files)
        
        folder_count = folder_count + 1
        #sys.exit()
    return data_trimmed_dir,fastqc_dir,results_dir
    
#salmon_index()
#kallisto_index()
data_trimmed_dir,fastqc_dir,results_dir = run_pipeline()    
