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

#usage example: Rscript Ref_dl_list_generator.R ref_list_files=ref_list_files.txt ref_tb=Ecoli_strs.csv name_out=Ecoli_ref_xt_EST
#Inputs: 1. the list of files(the converted bench table; the output matrixs to be evaluated);2. the download list of all the candidate strs;
#the output will be a list of ref strs to be downloaded. 
#library(tidyverse)
#check arguments
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
if (!exists("ref_list_files")) {
    stop(paste("\nRscript Ref_dl_list_generator.R ref_list_files=",ref_list_files," \nWarning: Usage: The file list with lists of ref strs is not exist, please check the path. \n\n",sep=""))
}

if (!exists("ref_tb")) {
    stop(paste("\nRscript Ref_dl_list_generator.R ref_tb=",ref_tb," \nWarning: Usage: The reference table downloaded from NCBI for the target strain is not exist, please check the path. \n\n",sep=""))
}
#ref_list_files <- "ref_list_files.txt"
#ref_tb <- "Ecoli_strs.csv"
#name_out <- "Ecoli_ref_xt_EST"

ref_list_f <- read.table(ref_list_files,as.is=TRUE,stringsAsFactors = FALSE)
ref_list=""
for(ref_f in ref_list_f[,1]){
    #only csv or tab dilim files will be supoorted.
    SIM_in<-read.table(ref_f,as.is=TRUE,stringsAsFactors = FALSE)
    SIM_ref_ID <- SIM_in[,1]
    if (sum(substr(SIM_ref_ID,1,3)!="GCF")>0){
            if (sum(substr(SIM_ref_ID,1,3)!="GCA")>0){
                stop(paste("\nThe simulation table (",ref_f,") is not annotated in GCF or GCA format \n Please convert it into the right annotation, Exit! \n\n",sep=""))
            }
    }
    #because GCF is required as panphlan input, so GCF is used instead of GCA
    if (sum(substr(SIM_ref_ID,1,3)=="GCA")>0){
            SIM_ref_ID <- chartr('GCA','GCF',SIM_ref_ID)
    }
    for (i in 1:length(SIM_ref_ID)){
            SIM_ref_ID[i]<-strsplit(SIM_ref_ID[i],split='.',fixed=TRUE)[[1]][1]
    }
    if(ref_list==""){                
            ref_list<-SIM_ref_ID
    }else{
            ref_list<-union(ref_list,SIM_ref_ID)
    }
}

#GCF is wilth RefSeq.FTP
ref_tb_in <- read.csv(ref_tb)
ref_GCF <- as.matrix(ref_tb_in$RefSeq.FTP)
rownames(ref_GCF)<- ref_GCF[,1]
ref_bac <- as.matrix(ref_tb_in$RefSeq.FTP)
rownames(ref_bac)<- ref_GCF[,1]
emp_row_list<-""
for (i in 1:length(ref_GCF[,1])){
    if(ref_GCF[i,1]!=""){
        str_tmp <- strsplit(ref_GCF[i,1],split='/',fixed=TRUE)[[1]]
    }else{
        if(ref_tb_in$RefSeq.FTP[i]!=""){
                str_tmp <- strsplit(ref_tb_in$RefSeq.FTP[i],split='/',fixed=TRUE)[[1]]
        }else{
                if(emp_row_list==""){
                        emp_row_list=i
                }else{
                        emp_row_list=c(emp_row_list,i)
                }
                next
        }
    }
    str_tmp1 <-str_tmp[length(str_tmp)]
    rownames(ref_GCF)[i]<-strsplit(str_tmp1,split='.',fixed=TRUE)[[1]][1]
    ref_GCF[i,1]<-paste(ref_GCF[i,1],"/",str_tmp1,"_genomic.fna.gz",sep="")
    rownames(ref_bac)[i]<-strsplit(str_tmp1,split='.',fixed=TRUE)[[1]][1]
    ref_bac[i,1]<-paste(str_tmp1,"_genomic.fna.gz",sep="")
}
if(emp_row_list!=""){
        ref_GCF<-ref_GCF[-emp_row_list,]
        ref_bac<-ref_bac[-emp_row_list,]
}



ref_ID_intersect <- intersect(names(ref_GCF),ref_list)
miss_str <- setdiff(ref_list,ref_ID_intersect)
if(length(miss_str)>0){
        print(paste("\nWarning:the following strings (",paste(miss_str,collapse=";"),") are not in the ref_table (",ref_tb,") \n Please double check it before moving on! \n\n",sep=""))
        write.table(miss_str, file=paste(name_out, "_miss_strs.txt",sep = ""),sep=",", quote=F, row.names= FALSE, col.names= FALSE, fileEncoding="UTF-8" )
}
ref_dl_list <- ref_GCF[ref_ID_intersect]
bac_list <- ref_bac[ref_ID_intersect]
write.table(ref_dl_list, file=paste(name_out, "_ref_dl_list.txt",sep = ""),sep=",", quote=F, row.names= FALSE, col.names= FALSE, fileEncoding="UTF-8" )
write.table(bac_list, file=paste(name_out, "_bac_list.txt",sep = ""),sep=",", quote=F, row.names= FALSE, col.names= FALSE, fileEncoding="UTF-8" )

