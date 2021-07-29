#coding:utf-8
from __future__ import print_function
import glob,os,re,argparse

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

#This script will simulate reads (both errfree and with err) of the selectied genome of a strain with the following parameters in ART.

#This script starts from a fa (not gz) file of a strain
#Output folder need to be inputted
#Requirements:
#ART installed
#samtools installed 
#Minimum usage example:python /home/yxtan/StrainSIM/Strain_read_sim.py -i genome_fa --sf script_folder -f fold_cov --pf prefix -o outputfolder
#Full usage example:python /home/yxtan/StrainSIM/Strain_read_sim.py -i genome_fa --sf script_folder --ss seq_model -l length -f fold_cov -m mean_frag_size -s std_frag_size --pf prefix -o outputfolder -r replicates
"""

def two_steps(script_fd, seq_md, fa, read_len, fold_cov, mean_frag, std_frag, OutN, sam_num, prefix_name):
    """
    Step 1: generate reads with ART
    """    
    os.system('%s/art_illumina -ss %s -ef -i %s -p -l %s -f %s -m %s -s %s -o %s/%s/%s 1>>%s/%s/ART_run.log' % (script_fd, seq_md, fa, read_len, fold_cov, mean_frag, std_frag, OutN, sam_num, prefix_name, OutN, sam_num))
    """
    Step 2: generate errfree reads by samtools from errfree.sam
    """    
    os.system('samtools fastq -1 %s/%s/%s-errFree_1.fq -2 %s/%s/%s-errFree_2.fq %s/%s/%s_errFree.sam 2>>%s/%s/sam2fq_run.log' % (OutN, sam_num, prefix_name, OutN, sam_num, prefix_name, OutN, sam_num, prefix_name, OutN, sam_num))
    #key outputs: are prefix_name1.fq, prefix_name2.fq, prefix_name_errFree_1.fq, prefix_name_errFree_2.fq


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='fa', type=str, required=True,
                        help="the path of the genome fa of the strain")
    parser.add_argument('-d', '--sf', dest='sf', type=str, required=True,
                        help="the path of art_illumina scripts")
    parser.add_argument('-a', '--ss', dest='ss', type=str, required=False, default='MSv3',
                        help="The sequencer modle of ART")
    parser.add_argument('-l', '--length', dest='rl', type=int, required=False, default='150',
                        help="The length of simulated reads")
    parser.add_argument('-f', '--fold_cov', dest='fc', type=float, required=True,
                        help="It is the fold_coverage of the sample")
    parser.add_argument('-m', '--mean_frag', dest='mf', type=int, required=False, default='270',
                        help="The expected mean of length of simulated fragment size")
    parser.add_argument('-s', '--std_frag', dest='stdF', type=int, required=False, default='27',
                        help="The expected std of length of simulated fragment size")
    parser.add_argument('-p', '--pf', dest='pf', type=str, required=True,
                        help="The name of outputfile")
    parser.add_argument('-o', '--output', dest='OutN', type=str, required=True,
                        help="The name of outputfile. Should not contain '_' in the name ")
    parser.add_argument('-r', '--rep_num', dest='rn', type=int, required=False, default='3',
                        help="The number of replicates samples")
    args = parser.parse_args()
    #print('Usage example with minimum parameters:  python /home/yxtan/StrainSIM/Strain_read_sim.py -i genome_fa --sf script_folder -f fold_cov --pf prefix -o outputfolder')
    fa = os.path.abspath(args.fa)
    script_fd = os.path.abspath(args.sf)
    seq_md = args.ss
    read_len = str(args.rl)
    fold_cov = str(args.fc)
    mean_frag = str(args.mf)
    std_frag = str(args.stdF)
    prefix_n = args.pf
    OutN = os.path.abspath(args.OutN)
    rep_num = args.rn
    
    if not os.path.isfile(fa):
        print('Input sample fa is not exist in Strain_read_sim.py; Exit now.')
        exit(0)
    
    if not os.path.isdir(script_fd):
        print('The folder of ART scripts is not exist in Strain_read_sim.py; Exit now.')
        exit(0)
    
    if not os.path.isdir(OutN):
        print('The output folder is not exist in Strain_read_sim.py; Exit now.')
        exit(0)
    
    #run for n replicates
    prefix_name="%s_f%sX" %(prefix_n, fold_cov)
    for sam_num in range(0,rep_num):
        if not os.path.isdir('%s/%s/' % (OutN, sam_num)):
            os.system('mkdir %s/%s/' % (OutN, sam_num))
        
        two_steps(script_fd, seq_md, fa, read_len, fold_cov, mean_frag, std_frag, OutN, str(sam_num), prefix_name)
