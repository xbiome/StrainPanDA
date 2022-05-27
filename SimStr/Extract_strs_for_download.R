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


#usage example: Rscript Extract_strs_for_download.R refstr=Ecoli_strs.csv fileout=Ecoli_target_str_info
#use three input files to get candidate strs from NCBI.
#the output will be in the local folder with the header of fileout
#library(tidyverse) is requried

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

if (!exists("refstr")) {
  stop("\nRscript Extract_strs_for_download.R refstr=Ecoli_strs.csv \nWarning: Usage: refstr file is not exist, please check the path. \n\n")
}



SRA_str_WGS_pair_lv_ml <- read.csv(refstr,as.is=TRUE,stringsAsFactors = FALSE )
#out files for downstream
library(tidyverse)
SRA_str_WGS_pair_lv_ml[,"RefSeq.FTP"]<-as.character(SRA_str_WGS_pair_lv_ml[,"RefSeq.FTP"])
row_i=1
while(row_i>0){
  if (nchar(SRA_str_WGS_pair_lv_ml[row_i,"RefSeq.FTP"])>0){
    into_name <- unlist(strsplit(SRA_str_WGS_pair_lv_ml[row_i,"RefSeq.FTP"],split="/"))
    break
  }
  row_i=row_i+1
}

into_name[2] <- "empty"
FTP_sep <- separate(data = SRA_str_WGS_pair_lv_ml, col = "RefSeq.FTP", into = into_name, sep = "/")
gz_att <-unlist(rep("_genomic.fna.gz",times=dim(FTP_sep)[1]))
gff_att <-unlist(rep("_genomic.gff.gz",times=dim(FTP_sep)[1]))
gz_att2 <- unlist(rep("/",times=dim(FTP_sep)[1]))
ref_dl_list <- cbind(SRA_str_WGS_pair_lv_ml[,"RefSeq.FTP"],gz_att2,FTP_sep[,into_name[length(into_name)]],gz_att)
gff_dl_list <- cbind(SRA_str_WGS_pair_lv_ml[,"RefSeq.FTP"],gz_att2,FTP_sep[,into_name[length(into_name)]],gff_att)
#remove empty rows in RefSeq_FTP
row_have_ref_seq <- !is.na(ref_dl_list[,3])
ref_dl_list<-ref_dl_list[!is.na(ref_dl_list[,3]),]
ref_dl_list_str<-paste(ref_dl_list[,1],ref_dl_list[,2],ref_dl_list[,3],ref_dl_list[,4],sep="")
gff_dl_list<-gff_dl_list[!is.na(gff_dl_list[,3]),]
gff_dl_list_str<-paste(gff_dl_list[,1],gff_dl_list[,2],gff_dl_list[,3],gff_dl_list[,4],sep="")
bac_list<-paste(ref_dl_list[,3],ref_dl_list[,4],sep="")
write.table(ref_dl_list_str[!duplicated(ref_dl_list_str)], file=paste(fileout, "_ref_dl_list.txt",sep = ""),sep=",", quote=F, row.names= FALSE, col.names= FALSE, fileEncoding="UTF-8" )
write.table(bac_list[!duplicated(bac_list)], file=paste(fileout, "_bac_list.txt",sep = ""),sep=",", quote=F, row.names= FALSE, col.names= FALSE, fileEncoding="UTF-8" )
write.table(gff_dl_list_str[!duplicated(gff_dl_list_str)], file=paste(fileout, "_gff_dl_list.txt",sep = ""),sep=",", quote=F, row.names= FALSE, col.names= FALSE, fileEncoding="UTF-8" )



