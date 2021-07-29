#!/bin/bash

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


#This script will run the evaluation on Xtrain.
#using scripts: Xtrain_converter.R and SIM_evaluate_unifrac_xtrain.R
#This script starts from three inputs: rds from xtrain; panphlan_ref; bench_table
#output will be in the out_fd
#Library file requirement:
#NOTE1: Rscript must be activated, with tidyverse\ggplot2\ctc\pheatmap\vegan\library(phyloseq) installed

if [ $# -ne 6 ]
then
    echo ""
    echo "Usage: ref_tree_generator.sh ref_list_files ref_tb name_out StrSIM_path out_fd nproc"
    echo "Example:sh /home/yxtan/StrainSIM/ref_tree_generator.sh ref_list_files.txt Ecoli_strs.csv Ecoli_ref_xt_EST /home/yxtan/StrainSIM . 8"
    echo ""
    echo "ref_list_files - the file with list of ref strains of each methods and simulation."
    echo "ref_tb - the full strain path file for the target species."
    echo "name_out - ."
    echo "StrSIM_path - ."
    echo "out_fd - target folder, could use . for current folder."
    echo "nproc - the number of processcors to be used in parsnp."       
    exit 1
fi

#set para
ref_list_files=$1 
ref_tb=$2 
name_out=$3
StrSIM_path=$4
out_fd=$5
nproc=$6


#check files
if [ ! -s $ref_list_files ] 
then
  echo ""
  echo "Warning: The file $ref_list_files does not exist in ref_tree_generator.sh, exit."
  echo ""
  exit 1
fi

if [ ! -s $ref_tb ] 
then
  echo ""
  echo "Warning: The file $ref_tb does not exist in ref_tree_generator.sh, exit."
  echo ""
  exit 1
fi


#Check folders
if [ ! -d $StrSIM_path ] 
then
  echo ""
  echo "Warning: The directory $StrSIM_path does not exist in ref_tree_generator.sh, exit."
  echo ""
  exit 1
fi

if [ ! -d $out_fd ] 
then
  echo ""
  echo "Warning: The directory $out_fd does not exist in ref_tree_generator.sh, generate it."
  mkdir $out_fd
  echo ""

fi

name_out_fd=$out_fd"/"$name_out
echo "Step1:Run R to generate ref Strain download list "
Rscript $StrSIM_path"/"Ref_dl_list_generator.R ref_list_files=$ref_list_files ref_tb=$ref_tb name_out=$name_out_fd

ref_out_fd=$out_fd"/Refs_fa"
if [ ! -d $ref_out_fd ] 
then
  echo ""
  echo "Warning: The directory $ref_out_fd does not exist in ref_tree_generator.sh, generate it."
  mkdir $ref_out_fd
  echo ""

fi
echo "Step2:download refs and do parsnp"
wget -i $name_out_fd"_ref_dl_list.txt" -P $ref_out_fd

cd $ref_out_fd
gunzip *.gz
cd ..

time parsnp -r ! -d $ref_out_fd -p $nproc -c #-c is necessary to avoid lost of branch