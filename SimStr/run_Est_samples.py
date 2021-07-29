#coding:utf-8
from __future__ import print_function
import glob,os,re,argparse
import pandas as pd

"""
# Copyright {2021} Yuxiang Tan
# This file is part of SimStr. 
#
# SimStr is a pipeline to generate simulation datasets for evaluation on strain analysis from metagenomic data.
#
# SimStr is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SimStr is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with SimStr.  If not, see <https://www.gnu.org/licenses/>.


#This script will read in the table of samples fqs and run under the docker of Est, as a result, this script must be in the working folder.

#This script starts from a csv table. no header needed, column1 is fq1, column2 is fq2, column3 is the name of the sample. Only paired end is supported now.

#Output folder need to be inputted
#subfolder structure:
#1.distribution of each sample is in dist_csv folder
#2.raw simulated reads of strains are in raw_sim folder. Each sample has its prefix name as subfolder name
#3.merged fqs are in final_fq samples. Its prefix name_type as file name. errfree means without sequencing err; seqerr means with sequencing err; seqfile means from real sequencing files. rX represent the number of the replicate sample.

#Requirements:
#All installed in the docker already


#Full usage example:python run_Est_samples.py -i csv_table -p index_path -d dgrp_path -o outputfolder

Update log 2021-05-26 Yuxiang
1. add the multiprocessor option for all tools

"""

def read_parse_csv(csv_table, index_path, dgrp_path, OutN, multiprocessor):
    new_data=pd.read_csv(csv_table,header=None)
    col_num= len(new_data.columns)
    row_num=len(new_data.index)
    #run for each sample
    for r_n in range(0,row_num):
        fq1=new_data.loc[r_n,0]
        fq2=new_data.loc[r_n,1]
        sample_name=new_data.loc[r_n,2]
        print("run sample %s"% (sample_name))
        os.system('sickle pe -f %s -r %s -t sanger -o %s/intermedia_f/%s.trim1.fq -p %s/intermedia_f/%s.trim2.fq -s %s/intermedia_f/%s.singles.fq -q 20' % (fq1, fq2, OutN, sample_name, OutN, sample_name, OutN, sample_name))
        os.system('bowtie2 -p %s --very-fast --no-unal -x %s -1 %s/intermedia_f/%s.trim1.fq -2 %s/intermedia_f/%s.trim2.fq -S %s/intermedia_f/%s.sam' % (multiprocessor,index_path, OutN, sample_name, OutN, sample_name, OutN, sample_name))
        os.system('samtools view --threads %s -b %s/intermedia_f/%s.sam > %s/intermedia_f/%s.bam' % (multiprocessor,OutN, sample_name, OutN, sample_name))
        os.system('samtools sort --threads %s %s/intermedia_f/%s.bam -o %s/intermedia_f/%s.sorted.bam' % (multiprocessor,OutN, sample_name, OutN, sample_name))
        os.system('samtools index %s/intermedia_f/%s.sorted.bam' % (OutN, sample_name))
        os.system('strainest est -t %s %s %s/intermedia_f/%s.sorted.bam %s/%s' % (multiprocessor, dgrp_path, OutN, sample_name, OutN, sample_name))
        





if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='csvt', type=str, required=True,
                        help="the path of the genome fa of the strain")
    parser.add_argument('-p', '--index', dest='idp', type=str, required=True,
                        help="the path of the bowtie index of ref genome, such as P_acnes/bowtie/align")
    parser.add_argument('-d', '--dgrp', dest='dp', type=str, required=True,
                        help="the path of the dgrp of ref genome, such as P_acnes/snp_clust.dgrp")
    parser.add_argument('-o', '--output', dest='OutN', type=str, required=True,
                        help="The name of outputfile")
    parser.add_argument('-t', '--threads', dest='threads', type=int, required=False,default=1,
                        help="number of threads to use")
    args = parser.parse_args()
    csv_table = os.path.abspath(args.csvt)
    index_path = os.path.abspath(args.idp)
    dgrp_path = os.path.abspath(args.dp)
    OutN = os.path.abspath(args.OutN)
    multiprocessor=args.threads
    #需要检查输入的参数是否正确，主要是路径是否存在    
    if not os.path.isfile(csv_table):
        print('Input csv table is not exist in Strain_list_merger.py; Exit now.')
        exit(0)
    
    if not os.path.isfile("%s.1.bt2" %(index_path)):
        print('The path of the bowtie index of ref genome is not exist in Strain_list_merger.py; Exit now.')
        exit(0)
    
    if not os.path.isfile(dgrp_path):
        print('The path of the dgrp of ref genome is not exist in Strain_list_merger.py; Exit now.')
        exit(0)
    
    if not os.path.isdir(OutN):
        print('The output folder is not exist in Strain_list_merger.py; Exit now.')
        exit(0)
    
    if not os.path.isdir("%s/intermedia_f" %(OutN)):
        os.system('mkdir %s/intermedia_f' %(OutN))
    
    #run the main function
    read_parse_csv(csv_table, index_path, dgrp_path, OutN, multiprocessor)
