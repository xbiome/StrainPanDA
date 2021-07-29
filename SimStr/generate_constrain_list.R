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

#usage example: Rscript Path_to_SimStr/generate_constrain_list.R wfd=folder_fq_files fq_reg=recognizer of the first end fq file(for example _1.fq)
#This is for ConStrains and PStrain

for (e in commandArgs()) {
        ta = strsplit(e,"=",fixed=TRUE)
        if(! is.na(ta[[1]][2])) {
                temp = ta[[1]][2]
                if(substr(ta[[1]][1],nchar(ta[[1]][1]),nchar(ta[[1]][1])) == "I") {
                temp = as.integer(temp)
                }
        if(substr(ta[[1]][1],nchar(ta[[1]][1]),nchar(ta[[1]][1])) == "N") {
                temp = as.numeric(temp)
                }
        assign(ta[[1]][1],temp)
        } else {
        assign(ta[[1]][1],TRUE)
        }
}

#check whether the file in is exist
if (!exists("wfd")) {
    stop(paste("\nRscript generate_Est_list.R wfd=",wfd," \nWarning: Usage: the target folder is not exist, please check the path. \n\n",sep=""))
}
if (substr(wfd, 1,1)=="/"){
    f_path=wfd
}else{
    cDIR=getwd()
    f_path=paste(cDIR,"/",wfd,sep = "")
}

if (file.exists("fq_list.txt")) {
    file.remove("fq_list.txt")
}

constrain_list <- function(f_path,fq_reg){
  file_list=dir(path=f_path,pattern = paste("*",fq_reg,sep=""))
  out_list <- matrix("",nrow=length(file_list),ncol=3)
  for (fq_n in 1:length(file_list)) {
    sm_name <- strsplit(file_list[fq_n], fq_reg)[[1]]
    Ffq_path <- paste(f_path,"/",sm_name,fq_reg,sep="")
    Rfq_path <- paste(f_path,"/",sm_name,chartr(old = "1",new = "2",fq_reg),sep="")
    out_list[fq_n,1]<-Ffq_path
    out_list[fq_n,2]<-Rfq_path
    out_list[fq_n,3]<-sm_name
  cat( paste("//\nsample:",sm_name,"\n",sep=""), file="fq_list.txt" ,append=TRUE)
  cat( paste("fq1:",Ffq_path,"\n",sep=""), file="fq_list.txt",append=TRUE )
  cat( paste("fq2:",Rfq_path,"\n",sep=""), file="fq_list.txt",append=TRUE )
  }
  return(out_list)
}

outfile <- constrain_list(f_path,fq_reg)