#!/bin/bash

# Copyright {2021} Yuxiang Tan
# This file is part of SimStr. 
#
# SimStr is a pipeline to generate simulation datasets for evaluation on strain analysis from metagenomic data.
#
# SimStr is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# SimStr is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with SimStr.  If not, see <https://www.gnu.org/licenses/>.

#This script will run the strain selection for simulation data.
#This script starts from three info files from NCBI
#Library file requirement:
#NOTE1: Rscript must be activated, with tidyverse\ggplot2\ggpubr\pheatmap\vegan installed
#NOTE1: fastani must be activated,
if [ $# -ne 4 ]
then
    echo ""
    echo "Usage: XT_ref_str_select_ANI_mash.sh ref_str_info_csv fileout_header StrSIM_path range_type "
    echo "Example:sh /home/yxtan/StrainSIM/XT_ref_str_select_ANI_mash.sh Ecoli_strs.csv Ecoli_target_str_info /home/yxtan/StrainSIM wide "
    echo ""
    echo "ref_str_info_csv - ."
    echo "fileout_header - must be in the following format: genus-species-yearmonth(in number)."
    echo "StrSIM_path - ."
    echo "range_type - 'wide' which will get the combination with the widest range of ANIs; or 'max ANI value (e.g. 98)' which will get all the strs combinations with ANI smaller or equal to that value, the bigger this value, the more strs you will get."
       
    exit 1
fi

#set para
refstr=$1
fileout=$2
StrSIM_path=$3
range_type=$4

if [ ! -s $refstr ] 
then
  echo ""
  echo "Warning: The file $refstr does not exist in XT_ref_str_select_ANI_mash.sh, exit."
  echo ""
  exit 1
fi

#Check folders
if [ ! -d $StrSIM_path ] 
then
  echo ""
  echo "Warning: The directory $StrSIM_path does not exist in XT_ref_str_select_ANI_mash.sh, exit."
  echo ""
  exit 1
fi

#Run R to extract candidates from 3 NCBI files
Rscript $StrSIM_path"/Extract_strs_for_download.R" refstr=$refstr fileout=$fileout

#check out file from previous step
if [ ! -s $fileout"_ref_dl_list.txt" ] 
then
  echo ""
  echo "Warning: The file "$fileout"_ref_dl_list.txt does not exist in XT_ref_str_select_ANI_mash.sh; There must be no strains pass the filter, exit."
  echo ""
  exit 1
fi

if [ ! -d "strain-ref" ] 
then
  mkdir strain-ref
fi


#run fastANI
cd strain-ref
wget -i "../"$fileout"_ref_dl_list.txt"
bac_list=$fileout"_bac_list.txt"
mash sketch -o ../sketch *.fna.gz
cd ..
mash dist -t sketch.msh sketch.msh>mash.dist

#range_type="98" 
ANI_tb="mash.dist"
#get SIM strs and their GCF and GFF download list
Rscript $StrSIM_path"/Extract_strs_by_ANI_mash.R" range_type=$range_type ANI_tb=$ANI_tb ref_dl_tb=$fileout"_ref_dl_list.txt" gff_dl_tb=$fileout"_gff_dl_list.txt"

if [ ! -d "strain-pan" ] 
then
  mkdir strain-pan
fi

#download the fna and gff files #add -c option will make the download more robust
cd strain-pan
wget -c -i "../"$fileout"_ref_dl_list.txt_sel.txt"
wget -c -i "../"$fileout"_gff_dl_list.txt_sel.txt"
cd ..

#run panphlan
time panphlan_pangenome_generation.py -c $fileout --i_fna strain-pan/ --i_gff strain-pan/ -o $fileout

