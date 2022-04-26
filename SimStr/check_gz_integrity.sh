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


#This script will check the integrity of all gz files in the target folder, and remove the files with decompression errors.
#output will be the broken_files.txt in the target folder, and the broken_ID.txt can be used to grep -f and get the subset of list to re-download.

if [ $# -ne 1 ]
then
    echo ""
    echo "Usage: check_gz_integrity.sh in_folder"
    echo "Example:sh /home/yxtan/StrainSIM/check_gz_integrity.sh strain-ref  "
    echo ""
    echo "in_folder - the target folder such as strain-ref ."
        exit 1
fi

#set para
in_folder=$1 


#Check folders
if [ ! -d $in_folder ] 
then
  echo ""
  echo "Warning: The directory $in_folder does not exist in check_gz_integrity.sh, exit."
  echo ""
  exit 1
fi

cd $in_folder
for i in `find . -name "*.gz*"`; 
 do 
 gunzip -t $i 2>$i.err
done

find . -name "*err" -type f -size +0c > broken_files.txt
cut -f2 -d"/" broken_files.txt | cut -f1 -d"." | sort -u > broken_ID.txt
rm *gz*.err
sed -i "s/.err//g" broken_files.txt

for i in $(cat broken_files.txt);do 
    rm $i
done