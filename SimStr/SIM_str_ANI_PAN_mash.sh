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

#This script will select the candidate strains with paired WGS data for simulation.
#This script starts from three info files from NCBI
#Library file requirement:
#NOTE1: Rscript must be activated, with tidyverse\ggplot2\ggpubr\pheatmap\vegan installed
#NOTE2: fastani must be activated,
#NOTE3: panphlanv1.2.4 must be activated,
#NOTE4: sra-tookit must be activated 
#NOTE5: Must not run multiple simulation in the same folder

if [ $# -ne 7 ]
then
    echo ""
    echo "Usage: SIM_str_ANI_PAN_mash.sh sra_summary_csv sra_runinfo_csv strs.csv fileout_header Path_to_SimStr range_type min_ANI"
    echo "Example:sh /home/yxtan/StrainSIM/SIM_str_ANI_PAN_mash.sh sra_result.csv SraRunInfo.csv Ecoli_strs.csv Ecoli_target_str_info /home/yxtan/StrainSIM wide 95"
    echo ""
    echo "sra_result.csv: from https://www.ncbi.nlm.nih.gov/sra/, select the target species and select the sent to (top right) -> file -> save summary.."
    echo "SraRunInfo.csv: from https://www.ncbi.nlm.nih.gov/sra/, select the target species and select the sent to (top right) -> file -> save Runinfo."
    echo "strs.csv: from https://www.ncbi.nlm.nih.gov/genome/browse/#!/prokaryotes/, select the target species; Select all the strains or subset of the target strains; click "
    echo "fileout_header: name of output files. For example: Ecoli99_target_str_info ."
    echo "Path_to_SimStr: the full path of SimStr scripts folder ."
    echo "range_type - 'wide' which will get the combination with the widest range of ANIs; or 'max ANI value (e.g. 98)' which will get all the strs combinations with ANI smaller or equal to that value, the bigger this value, the more strs you will get."
    echo "min_ANI - the min ANI allowed for simulation, ANI>=95 is considered as within species."
   
    exit 1
fi

#set para
srafile=$1 
runinfofile=$2 
refstr=$3
fileout=$4
StrSIM_path=$5
range_type=$6
min_ANI=$7

#check files
if [ ! -s $srafile ] 
then
  echo ""
  echo "Warning: The file $srafile does not exist in SIM_str_ANI_PAN_mash.sh, exit."
  echo ""
  exit 1
fi

if [ ! -s $runinfofile ] 
then
  echo ""
  echo "Warning: The file $runinfofile does not exist in SIM_str_ANI_PAN_mash.sh, exit."
  echo ""
  exit 1
fi

if [ ! -s $refstr ] 
then
  echo ""
  echo "Warning: The file $refstr does not exist in SIM_str_ANI_PAN_mash.sh, exit."
  echo ""
  exit 1
fi

#Check folders
if [ ! -d $StrSIM_path ] 
then
  echo ""
  echo "Warning: The directory $StrSIM_path does not exist in SIM_str_ANI_PAN_mash.sh, exit."
  echo ""
  exit 1
fi

#Run R to extract candidates (with pared WGS data) from 3 NCBI files
Rscript $StrSIM_path"/Extract_candidate_strs.R" srafile=$srafile runinfofile=$runinfofile refstr=$refstr fileout=$fileout

#generate gff download list
cp $fileout"_ref_dl_list.txt" $fileout"_gff_dl_list.txt"
sed -i 's/.fna.gz/.gff.gz/g' $fileout"_gff_dl_list.txt"
#check out file from previous step
if [ ! -s $fileout"_ref_dl_list.txt" ] 
then
  echo ""
  echo "Warning: The file "$fileout"_ref_dl_list.txt does not exist in SIM_str_ANI_PAN_mash.sh; There must be no strains pass the filter, exit."
  echo ""
  exit 1
fi

if [ ! -d "strain-ref" ] 
then
  mkdir strain-ref
fi


#download the files for mash and panphlan
cd strain-ref
wget -i "../"$fileout"_ref_dl_list.txt"
wget -i "../"$fileout"_gff_dl_list.txt"
bac_list=$fileout"_bac_list.txt"

#run mash for ANI
mash sketch -o ../sketch *.fna.gz
cd ..
mash dist -t sketch.msh sketch.msh>mash.dist

#run panphlan to generate pangenome
#docker run -u $(id -u):$(id -g) -v $PWD:/data -w /data 265157181057.dkr.ecr.cn-northwest-1.amazonaws.com.cn/panphlan:v1.2.4 /bin/bash -c "panphlan_pangenome_generation.py -c target_species --i_fna strain-ref/ --i_gff strain-ref/ -o target_database/"
panphlan_pangenome_generation.py -c target_species --i_fna strain-ref/ --i_gff strain-ref/ -o target_database/

#generate stats and the list of candidate strains automatically 
ANI_tb="mash.dist"
select_info=$fileout"_sel_col.csv"
Rscript $StrSIM_path"/Extract_SIM_strs_ANI_PAN_mash.R" range_type=$range_type ANI_tb=$ANI_tb select_info=$select_info pan_ref="target_database/panphlan_target-species_pangenome.csv" min_ANI=$min_ANI

#get paired WGS data
if [ ! -d "WGS_data" ] 
then
  mkdir WGS_data
fi

cd WGS_data
wget -i ../SIM_WGS_dl_list.txt

#for loop to fastq-dump --split-files
while read LINE;do fastq-dump --skip-technical --split-3 $LINE;done  < ../SIM_WGS_fastdump_list.txt
cd ..

while read LINE;do gunzip -k $LINE;done  < SIM_ref_unzip.txt

sed -i 's/.fna.gz/.fna/g' SIM_ART_matrix.csv