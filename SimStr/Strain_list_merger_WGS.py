#coding:utf-8
from __future__ import print_function
import glob,os,re,argparse
import pandas as pd
import random

"""
#   Copyright {2020} Yuxiang Tan
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

#This script will read in the table of simulation samples‘ distribution and generate separted tables for files guiding the running of Strain_read_sim.py to simulate reads.

#This script starts from a csv table.
#The format of table is as following:
#Columns:taxa name in GG format, strain_prefix(uniq among used strains), genome location of the WGS file header(must not be gz, and must be in the following format: header_1.fastq header_2.fastq), ref_size (used to calculate read number for depth, in base unit, can just "ls -l ref.fq" to get it), read_len (used to calculate read number for depth, in base unit, can just "head WGS.fq" to get it in the length para) 
distribution in sample 1,distribution in sample 2... distribution in sample N (must be number, 1 means the minimum depth, M means M fold to 1, the relative abundance is divided by the total folds of that sample)
#Rows: 
#row1:name of each sample in the distribution columns(will be used as subfolder name, must be unique and contain no '_' character);
#row2:minimum depth(fold_coverage) of each sample in the distribution columns(must be number);
#row3 - K: K-2 strains. If a strain is not in a sample, its distribution in that column is 0.
#Note: taxa name in GG format: should be like: k__Bacteria;p__Proteobacteria;c__Epsilonproteobacteria;o__Campylobacterales;f__Campylobacteraceae;g__Campylobacter;s__jejuni;str__NCTC12851

#Output folder need to be inputted
#subfolder structure:
#1.distribution of each sample is in dist_csv folder
#2.raw simulated reads of strains are in raw_sim folder. Each sample has its prefix name as subfolder name
#3.merged fqs are in final_fq samples. Its prefix name_type as file name. errfree means without sequencing err; seqerr means with sequencing err; WGSerr means from real sequencing files. rX represent the number of the replicate sample.
#Requirements:
#seq-tk activated
#pandas install
#Minimum usage example:python Path_to_SimStr/Strain_list_merger_WGS.py -i csv_table -o outputfolder -d Path_to_SimStr
#Full usage example:python /home/yxtan/StrainSIM/Strain_list_merger_WGS.py -i csv_table -o outputfolder -r replicates  -l seq_read_len -d /home/yxtan/StrainSIM

#update log:2021-03-11 by Yuxiang Tan
##recongize the read length of each WGS data automatically (rather than default input now, which will cause mistake when WGS of different strains had different read length)
###removed read_len parameter
###add the read_len para in the input csv, this para should be manually added by user by look at the header of fastq file or doing check by themselves either by fastqc like tools or randomely pick reads to check.
###as a result, in the read_parse_csv function, one more column (num 4) added in front of c_n.
###therefore, the bench_WGS_table_converter.R was modified for the added column as well.
###the description of this script was updated as well in the Columns section of input csv.
"""

def read_parse_csv(csv_table, OutN, rep_num):
    new_data=pd.read_csv(csv_table,header=None)
    col_num= len(new_data.columns)
    row_num=len(new_data.index)
    for c_n in range(5,col_num):
        prefix_name=prefix_name=new_data.loc[0,c_n]
        min_fold_cov=float(new_data.loc[1,c_n])
        sample_data=new_data[[0,1,2,3,4,c_n]]
        sample_data.to_csv('%s/dist_csv/%s_WGS.csv' % (OutN,prefix_name),header=None, index=False)
        out_folder_n='%s/raw_sim/%s/' % (OutN,prefix_name)
        out_fq_n='%s/final_fq/%s' % (OutN,prefix_name)
        if not os.path.isdir(out_folder_n):
            os.system('mkdir %s' % (out_folder_n))
        
        os.system('echo "%s"' % (prefix_name))
        #run for each sample
        for r_n in range(2,row_num):
            fq1='%s_1.fastq' % (new_data.loc[r_n,2])
            fq2='%s_2.fastq' % (new_data.loc[r_n,2])
            #check fq files
            if not os.path.isfile(fq1):
                print('Input fq1: %s is not exist in Strain_list_merger_WGS.py; Exit now.' % (fq1))
                exit(0)
            if not os.path.isfile(fq2):
                print('Input fq2: %s is not exist in Strain_list_merger_WGS.py; Exit now.' % (fq2))
                exit(0)
            r_dist=int(new_data.loc[r_n,c_n])
            fold_cov=min_fold_cov*r_dist
            ref_size=int(new_data.loc[r_n,3])
            read_len=int(new_data.loc[r_n,4])
            read_num= int(ref_size*fold_cov/int(read_len)/2)
            STR_n=new_data.loc[r_n,1]
            sample_name='%s_f%sX' % (STR_n,fold_cov)
            a = range(1,10000)
            ran_seed_list = random.sample(a, rep_num)
            
            if fold_cov>0:
                for sam_num in range(0,rep_num):
                    ran_seed = ran_seed_list[sam_num]
                    if not os.path.isdir('%s/%s/' % (out_folder_n, sam_num)):
                        os.system('mkdir %s/%s/' % (out_folder_n, sam_num))
                    
                    os.system('seqtk sample -s %s %s %s > %s/%s/%s_WGS_1.fq' % (ran_seed, fq1, read_num, out_folder_n, sam_num, sample_name))
                    os.system('seqtk sample -s %s %s %s > %s/%s/%s_WGS_2.fq' % (ran_seed, fq2, read_num, out_folder_n, sam_num, sample_name))
                    #merge files into final_fq
                    if r_n==2:
                        os.system('cat %s/%s/%s_WGS_1.fq > %s_WGSerr_r%s_1.fq' % (out_folder_n, sam_num, sample_name, out_fq_n, sam_num))
                        os.system('cat %s/%s/%s_WGS_2.fq > %s_WGSerr_r%s_2.fq' % (out_folder_n, sam_num, sample_name, out_fq_n, sam_num))
                        
                    else:
                        os.system('cat %s/%s/%s_WGS_1.fq >> %s_WGSerr_r%s_1.fq' % (out_folder_n, sam_num, sample_name, out_fq_n, sam_num))
                        os.system('cat %s/%s/%s_WGS_2.fq >> %s_WGSerr_r%s_2.fq' % (out_folder_n, sam_num, sample_name, out_fq_n, sam_num))
                        




if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='csvt', type=str, required=True,
                        help="the path of the genome fa of the strain")
    parser.add_argument('-D', '--sp', dest='SSIM', type=str, required=True,
                        help="the path of StrainSIM scripts")
    # parser.add_argument('-l', '--length', dest='rl', type=int, required=False, default='150',
    #                     help="The midium length of WGS reads, all input samples should be the sam")
    parser.add_argument('-o', '--output', dest='OutN', type=str, required=True,
                        help="The name of outputfile")
    parser.add_argument('-r', '--rep_num', dest='rn', type=int, required=False, default='3',
                        help="The number of replicates samples")
    args = parser.parse_args()
    #print('Usage example with minimum parameters:  python /home/yxtan/StrainSIM/Strain_read_sim.py -i genome_fa --sf script_folder -f fold_cov --pf prefix -o outputfolder')
    csv_table = os.path.abspath(args.csvt)
    script_fd = os.path.abspath(args.SSIM)
    #read_len = args.rl
    OutN = os.path.abspath(args.OutN)
    rep_num = args.rn
    
    if not os.path.isfile(csv_table):
        print('Input csv table is not exist in Strain_list_merger_WGS.py; Exit now.')
        exit(0)
        
    if not os.path.isdir(script_fd):
        print('The folder of StrainSIM scripts is not exist in Strain_list_merger_WGS.py; Exit now.')
        exit(0)

    if not os.path.isdir(OutN):
        print('The output folder is not exist in Strain_list_merger_WGS.py; Exit now.')
        exit(0)
    
    #genereate output folders
    if not os.path.isdir('%s/dist_csv/' % (OutN)):
        os.system('mkdir %s/dist_csv/' % (OutN))
    
    if not os.path.isdir('%s/raw_sim/' % (OutN)):
        os.system('mkdir %s/raw_sim/' % (OutN))
    
    if not os.path.isdir('%s/final_fq/' % (OutN)):
        os.system('mkdir %s/final_fq/' % (OutN))
    
    #run the main function
    read_parse_csv(csv_table, OutN, rep_num)

    #generate the coverted csv table： _WGSerr_
    outname = csv_table.split("/")[-1].split(".csv")[0]
    os.system('Rscript %s/bench_WGS_table_converter.R sim_table=%s add_name=_WGSerr_r rep_add_time=%s name_out=%s/%s_WGSerr_all' % (script_fd, csv_table, rep_num, OutN, outname))
    
