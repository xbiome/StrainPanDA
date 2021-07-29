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

#This script will read in the table of simulation samples‘ distribution and generate separted tables for files guiding the running of Strain_read_sim.py to simulate reads.

#This script starts from a csv table.
#The format of table is as following:
#Columns:taxa name in GG format, strain_prefix(uniq among used strains), genome location of the strain(must not be gz), distribution in sample 1,distribution in sample 2... distribution in sample N (must be number, 1 means the minimum depth, M means M fold to 1, the relative abundance is divided by the total folds of that sample)
#Rows: 
#row1:name of each sample in the distribution columns(will be used as subfolder name, must be unique and contain no '_' character);
#row2:minimum depth(fold_coverage) of each sample in the distribution columns(must be number);
#row3 - K: K-2 strains. If a strain is not in a sample, its distribution in that column is 0.
#Note: taxa name in GG format: should be like: k__Bacteria;p__Proteobacteria;c__Epsilonproteobacteria;o__Campylobacterales;f__Campylobacteraceae;g__Campylobacter;s__jejuni;str__NCTC12851

#Output folder need to be inputted
#subfolder structure:
#1.distribution of each sample is in dist_csv folder
#2.raw simulated reads of strains are in raw_sim folder. Each sample has its prefix name as subfolder name
#3.merged fqs are in final_fq samples. Its prefix name_type as file name. errfree means without sequencing err; seqerr means with sequencing err; seqfile means from real sequencing files. rX represent the number of the replicate sample.
#Requirements:
#ART installed
#samtools installed 
#pandas install
#Minimum usage example:python Path_to_SimStr/Strain_list_merger.py -i csv_table -D Path_to_SimStr --sf ART_folder -o outputfolder
#Full usage example:python Path_to_SimStr/Strain_list_merger.py -i csv_table -D Path_to_SimStr --sf ART_folder --ss seq_model -l length -m mean_frag_size -s std_frag_size -o outputfolder -r replicates
"""

def read_parse_csv(csv_table, script_fd, ART_fd, seq_md, read_len, mean_frag, std_frag, OutN, rep_num):
    new_data=pd.read_csv(csv_table,header=None)
    col_num= len(new_data.columns)
    row_num=len(new_data.index)
    for c_n in range(3,col_num):
        prefix_name=prefix_name=new_data.loc[0,c_n]
        min_fold_cov=float(new_data.loc[1,c_n])
        sample_data=new_data[[0,1,2,c_n]]
        sample_data.to_csv('%s/dist_csv/%s.csv' % (OutN,prefix_name),header=None, index=False)
        out_folder_n='%s/raw_sim/%s/' % (OutN,prefix_name)
        out_fq_n='%s/final_fq/%s' % (OutN,prefix_name)
        if not os.path.isdir(out_folder_n):
            os.system('mkdir %s' % (out_folder_n))
        
        os.system('echo "%s"' % (prefix_name))
        #run for each sample
        for r_n in range(2,row_num):
            fa=new_data.loc[r_n,2]
            r_dist=int(new_data.loc[r_n,c_n])
            fold_cov=min_fold_cov*r_dist
            STR_n=new_data.loc[r_n,1]
            sample_name='%s_f%sX' % (STR_n,fold_cov)
            if fold_cov>0:
                os.system('python %s/Strain_read_sim.py -i %s --sf %s -f %s --ss %s -l %s -m %s --std_frag %s --pf %s -o %s -r %s' % (script_fd, fa, ART_fd, fold_cov, seq_md, read_len, mean_frag, std_frag, STR_n, out_folder_n, rep_num))
                #merge files into final_fq
                for sam_num in range(0,rep_num):
                    if r_n==2:
                        os.system('cat %s/%s/%s1.fq > %s_seqerr_r%s_1.fq' % (out_folder_n, sam_num, sample_name, out_fq_n, sam_num))
                        os.system('cat %s/%s/%s2.fq > %s_seqerr_r%s_2.fq' % (out_folder_n, sam_num, sample_name, out_fq_n, sam_num))
                        os.system('cat %s/%s/%s_errFree_1.fq > %s_errfree_r%s_1.fq' % (out_folder_n, sam_num, sample_name, out_fq_n, sam_num))
                        os.system('cat %s/%s/%s_errFree_2.fq > %s_errfree_r%s_2.fq' % (out_folder_n, sam_num, sample_name, out_fq_n, sam_num))
                    
                    else:
                        os.system('cat %s/%s/%s1.fq >> %s_seqerr_r%s_1.fq' % (out_folder_n, sam_num, sample_name, out_fq_n, sam_num))
                        os.system('cat %s/%s/%s2.fq >> %s_seqerr_r%s_2.fq' % (out_folder_n, sam_num, sample_name, out_fq_n, sam_num))
                        os.system('cat %s/%s/%s_errFree_1.fq >> %s_errfree_r%s_1.fq' % (out_folder_n, sam_num, sample_name, out_fq_n, sam_num))
                        os.system('cat %s/%s/%s_errFree_2.fq >> %s_errfree_r%s_2.fq' % (out_folder_n, sam_num, sample_name, out_fq_n, sam_num))





if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='csvt', type=str, required=True,
                        help="the path of the genome fa of the strain")
    parser.add_argument('-D', '--sp', dest='SSIM', type=str, required=True,
                        help="the path of StrainSIM scripts")
    parser.add_argument('-d', '--sf', dest='sf', type=str, required=True,
                        help="the path of art_illumina scripts")
    parser.add_argument('-a', '--ss', dest='ss', type=str, required=False, default='HS25',
                        help="The sequencer modle of ART")
    parser.add_argument('-l', '--length', dest='rl', type=int, required=False, default='150',
                        help="The length of simulated reads")
    parser.add_argument('-m', '--mean_frag', dest='mf', type=int, required=False, default='270',
                        help="The expected mean of length of simulated fragment size")
    parser.add_argument('-s', '--std_frag', dest='stdF', type=int, required=False, default='27',
                        help="The expected std of length of simulated fragment size")
    parser.add_argument('-o', '--output', dest='OutN', type=str, required=True,
                        help="The name of outputfile")
    parser.add_argument('-r', '--rep_num', dest='rn', type=int, required=False, default='3',
                        help="The number of replicates samples")
    args = parser.parse_args()
    #print('Usage example with minimum parameters:  python /home/yxtan/StrainSIM/Strain_read_sim.py -i genome_fa --sf script_folder -f fold_cov --pf prefix -o outputfolder')
    csv_table = os.path.abspath(args.csvt)
    script_fd = os.path.abspath(args.SSIM)
    ART_fd = os.path.abspath(args.sf)
    seq_md = args.ss
    read_len = str(args.rl)
    mean_frag = str(args.mf)
    std_frag = str(args.stdF)
    OutN = os.path.abspath(args.OutN)
    rep_num = args.rn
  
    #check path, generate them if not exist    
    if not os.path.isfile(csv_table):
        print('Input csv table is not exist in Strain_list_merger.py; Exit now.')
        exit(0)
    
    if not os.path.isdir(script_fd):
        print('The folder of StrainSIM scripts is not exist in Strain_list_merger.py; Exit now.')
        exit(0)
    
    if not os.path.isdir(ART_fd):
        print('The folder of ART scripts is not exist in Strain_list_merger.py; Exit now.')
        exit(0)
    
    if not os.path.isdir(OutN):
        print('The output folder is not exist in Strain_list_merger.py; Exit now.')
        exit(0)
    
    #genereate output folders
    if not os.path.isdir('%s/dist_csv/' % (OutN)):
        os.system('mkdir %s/dist_csv/' % (OutN))
    
    if not os.path.isdir('%s/raw_sim/' % (OutN)):
        os.system('mkdir %s/raw_sim/' % (OutN))
    
    if not os.path.isdir('%s/final_fq/' % (OutN)):
        os.system('mkdir %s/final_fq/' % (OutN))
    
    #run the main function
    read_parse_csv(csv_table, script_fd, ART_fd, seq_md, read_len, mean_frag, std_frag, OutN, rep_num)

    #generate the coverted csv table： seqerr and errfree
    outname = csv_table.split("/")[-1].split(".csv")[0]
    os.system('Rscript %s/bench_table_converter.R sim_table=%s add_name=_seqerr_r rep_add_time=%s name_out=%s/%s_seqerr_all' % (script_fd, csv_table, rep_num, OutN, outname))
    os.system('Rscript %s/bench_table_converter.R sim_table=%s add_name=_errfree_r rep_add_time=%s name_out=%s/%s_errfree_all' % (script_fd, csv_table, rep_num, OutN, outname))
    
