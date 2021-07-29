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


#in the Est output folder.

#put all the abund.txt files into a single folder
mkdir abund_prof_merge
count_num=0
for file in `find . -name abund.txt`
do
let count_num=count_num+1
newfile="abund_prof_merge/profile_"$count_num".txt"
cp $file $newfile
#echo $newfile
done

#run merge_Est.R to merge all profile files into one matrix.
Rscript ./merge_Est.R