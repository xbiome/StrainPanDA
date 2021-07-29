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



#usage example: Rscript Extract_candidate_strs.R srafile=sra_result.csv runinfofile=SraRunInfo.csv refstr=strs.csv fileout=fileout_header
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
if (!exists("srafile")) {
	stop("\nRscript Extract_candidate_strs.R srafile=sra_result.csv \nWarning: Usage: srafile para is not provied, please check the para \n\n")
}

if (!exists("runinfofile")) {
  stop("\nRscript Extract_candidate_strs.R runinfofile=SraRunInfo.csv \nWarning: Usage: runinfofile para is not provied, please check the para \n\n")
}

if (!exists("refstr")) {
  stop("\nRscript Extract_candidate_strs.R refstr=strs.csv \nWarning: Usage: strs.csv para is not provied, please check the para \n\n")
}

#check whether the file is exist
if (!file.exists(srafile)) {
  stop("\nRscript Extract_candidate_strs.R srafile=",srafile," \nWarning: Usage: srafile file is not exist, please check the path. \n\n")
}

if (!file.exists(runinfofile)) {
  stop("\nRscript Extract_candidate_strs.R runinfofile=",runinfofile," \nWarning: Usage: runinfofile file is not exist, please check the path. \n\n")
}

if (!file.exists(refstr)) {
  stop("\nRscript Extract_candidate_strs.R refstr=",refstr," \nWarning: Usage: strs.csv file is not exist, please check the path. \n\n")
}




SRA_result <- read.csv(srafile,as.is=TRUE,stringsAsFactors = FALSE )
SRA_runinfo <- read.csv(runinfofile,as.is=TRUE,stringsAsFactors = FALSE )
ref_str <- read.csv(refstr,as.is=TRUE,stringsAsFactors = FALSE )

#merge tables
SRA_merge <- merge(SRA_result,SRA_runinfo, by.x="Experiment.Accession",by.y="Experiment")
SRA_str <- merge(SRA_merge,ref_str,by.x="BioSample",by.y="BioSample")
WGS_low=FALSE
WGS_pair_low=FALSE
#extract WGS and get level and model info if there are enough strs
SRA_str_WGS <- SRA_str[which(SRA_str$LibraryStrategy=="WGS" & SRA_str$LibraryStrategy==SRA_str$Library.Strategy),]

#check whether the number of strains are enough to do the simulation
if (dim(SRA_str_WGS)[1]<=4) {
  WGS_low=TRUE
  SRA_str_WGS <- SRA_str
}

SRA_str_WGS_pair <- SRA_str_WGS[which(SRA_str_WGS$LibraryLayout=="PAIRED"),] 

#check again after filtering whether the number of strains are enough to do the simulation
if ( sum(!duplicated(SRA_str_WGS_pair[,"RefSeq.FTP"]))<=4) {
  if ( sum(!duplicated(SRA_str_WGS[,"RefSeq.FTP"]))<=4) {
    WGS_low=TRUE
    SRA_str_WGS_pair <- SRA_str
  } else {
    WGS_pair_low=TRUE
    SRA_str_WGS_pair <- SRA_str_WGS
  }
}

#will select complete is enough genomes to be select from
if(sum(SRA_str_WGS_pair$Level=="Complete")>20){
  SRA_str_WGS_pair_lv <- SRA_str_WGS_pair[which(SRA_str_WGS_pair$Level=="Complete"),]
}else{SRA_str_WGS_pair_lv <- SRA_str_WGS_pair}

#will select same sequencer as MiSeq is enough genomes to be select from
if(sum(SRA_str_WGS_pair_lv$Model=="Illumina MiSeq")>20){
  SRA_str_WGS_pair_lv_ml <- SRA_str_WGS_pair_lv[which(SRA_str_WGS_pair_lv$Model=="Illumina MiSeq"),]
}else{SRA_str_WGS_pair_lv_ml <- SRA_str_WGS_pair_lv}

#out info files
write.table(SRA_str_WGS_pair_lv_ml, file=paste(fileout, "_full.csv",sep = ""),sep=",", quote=F, row.names= FALSE, col.names= TRUE, fileEncoding="UTF-8" )
write.table(SRA_str, file=paste(fileout, "_full_nofilter.csv",sep = ""),sep=",", quote=F, row.names= FALSE, col.names= TRUE, fileEncoding="UTF-8" )

selected_col <- c("Run","BioProject.x","TaxID","Strain","Level","Size.Mb.","GC.","Scaffolds","Library.Strategy","Model","LibraryLayout","InsertSize","download_path","Assembly","GenBank.FTP","RefSeq.FTP") #重新给重要信息进行归类和排序

write.table(SRA_str[,selected_col], file=paste(fileout, "_sel_col_nofilter.csv",sep = ""),sep=",", quote=F, row.names= FALSE, col.names= TRUE, fileEncoding="UTF-8" )

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
gz_att2 <- unlist(rep("/",times=dim(FTP_sep)[1]))
ref_dl_list <- cbind(SRA_str_WGS_pair_lv_ml[,"RefSeq.FTP"],gz_att2,FTP_sep[,into_name[length(into_name)]],gz_att)
#remove empty rows in refseq_FTP
row_have_ref_seq <- !is.na(ref_dl_list[,3])
ref_dl_list<-ref_dl_list[!is.na(ref_dl_list[,3]),]
ref_dl_list_str<-paste(ref_dl_list[,1],ref_dl_list[,2],ref_dl_list[,3],ref_dl_list[,4],sep="")
bac_list<-paste(ref_dl_list[,3],ref_dl_list[,4],sep="")
write.table(ref_dl_list_str[!duplicated(ref_dl_list_str)], file=paste(fileout, "_ref_dl_list.txt",sep = ""),sep=",", quote=F, row.names= FALSE, col.names= FALSE, fileEncoding="UTF-8" )
write.table(bac_list[!duplicated(bac_list)], file=paste(fileout, "_bac_list.txt",sep = ""),sep=",", quote=F, row.names= FALSE, col.names= FALSE, fileEncoding="UTF-8" )



bac_list_ID <- bac_list
for(i in 1:length(bac_list)){
  bac_list_ID[i] <- strsplit(bac_list[i], split = ".", fixed = TRUE)[[1]][1]
}

write.table(SRA_str_WGS_pair_lv_ml[row_have_ref_seq,"download_path"], file=paste(fileout, "_WGS_dl_list.txt",sep = ""),sep=",", quote=F, row.names= FALSE, col.names= FALSE, fileEncoding="UTF-8" )

#add download and panref keys into the table
final_out <- cbind(SRA_str_WGS_pair_lv_ml[row_have_ref_seq,selected_col],ref_dl_list_str,bac_list,bac_list_ID)

write.table(final_out, file=paste(fileout, "_sel_col.csv",sep = ""),sep=",", quote=F, row.names= FALSE, col.names= TRUE, fileEncoding="UTF-8" )


if (WGS_low) {
  stop("\nThere are not enough strains data from paired WGS with 1 scaffold for WGS simulation. \n\n")
} else {
  if (WGS_pair_low) {
    stop("\nAlthough strains data from paired WGS with 1 scaffold for WGS simulation are enough, the corresponding unique ref-seqs are not enough. Therefore, the WGS simulation is not possible. \n\n")
  }
}

