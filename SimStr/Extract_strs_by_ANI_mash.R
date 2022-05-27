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

#usage example: Rscript Extract_strs_by_ANI.R range_type=wide(or specific number N to get all strs with ANI lower than N) ANI_tb=Ecoli_target_str_info-fastani-1k.out ref_dl_tb=Ecoli-complete-20200907.csv_ref_dl_list.txt
#use the ref gene profile and ANI, select target strains.
#the output will be a table of selected strs for SIM
#library(tidyverse) ,library(vegan),library(pheatmap) library(ggplot2),library(ggpubr)are requried
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
if (!exists("ANI_tb")) {
  stop(paste("\nRscript Extract_strs_by_ANI.R ANI_tb=",ANI_tb," \nWarning: Usage: ANI matrix is not exist, please check the path. \n\n",sep=""))
}

if (!exists("ref_dl_tb")) {
  stop(paste("\nRscript Extract_strs_by_ANI.R ref_dl_tb=",ref_dl_tb," \nWarning: Usage: The download list of ref is not exist, please check the path. \n\n",sep=""))
}

if (!exists("gff_dl_tb")) {
  stop(paste("\nRscript Extract_strs_by_ANI.R gff_dl_tb=",gff_dl_tb," \nWarning: Usage: The download list of gff is not exist, please check the path. \n\n",sep=""))
}


library(tidyverse)
library(vegan)
library(pheatmap)


#analyze the ANI from mash of available WGS strs
ANI_all_strain <- read.table(ANI_tb,row.names=1)
ANI_names<-rownames(ANI_all_strain)
ANI_names_list <- matrix(0, nrow=length(ANI_names),ncol=1)
for(i in 1:length(ANI_names)){
ANI_names_list[i,1] <- strsplit(ANI_names[i], split = ".", fixed = TRUE)[[1]][1]
}

rownames(ANI_all_strain) <- ANI_names_list
colnames(ANI_all_strain) <- ANI_names_list
ANI_all_strain=(1-ANI_all_strain)*100


pdf(file = "ANI_all_strain_hist.pdf" , width = 11, height= 8.5)
hist(ANI_all_strain[ANI_all_strain<100])
dev.off()

#get the strs by ANI
if(range_type=="wide") {
    print("With 'wide_range', gete the list of strs with the widest ANI range and one str in each 1 ANI break.")
    ANI_sh_strain<-ANI_all_strain
    sh_ref<-as.character(colnames(ANI_all_strain))
    
    #get the most diverse distributed str by ANI
    dis_breaks <-seq(from=round(min(ANI_sh_strain)-1), to=max(ANI_sh_strain),by=1)
    hist_list <- matrix(nrow=length(sh_ref),ncol=length(dis_breaks)-1)
    colnames(hist_list)<-dis_breaks[-1]
    rownames(hist_list)<-rownames(ANI_sh_strain)
    for (row_name in rownames(ANI_sh_strain)) {
      hist_list[row_name,] <- hist(ANI_sh_strain[row_name,ANI_sh_strain[row_name,]<100], breaks = dis_breaks, plot = FALSE)$counts
    }
    hist_list_count <- (hist_list>0)*1
    hlc <- apply(hist_list_count,1,sum)
    hist_list_max<-hist_list[which(hlc==max(hlc)),]
    hist_list_max_dist<-t(apply(hist_list_max,1,function(x) x*as.numeric(colnames(hist_list_max))))
    target_str <- rownames(hist_list_max_dist)[which(rowSums(hist_list_max_dist)==max(rowSums(hist_list_max_dist)))][1]

    #get the strs
    SIM_strs <- rep(0,sum(hist_list_count[target_str,])+1)
    names(SIM_strs) <-c(0,colnames(hist_list_count)[hist_list_count[target_str,]>0])
    SIM_strs[1] <- target_str
    for (k in 2:(sum(hist_list_count[target_str,])+1)) {
      SIM_strs[k] <- colnames(ANI_sh_strain)[(ANI_sh_strain[target_str,]<as.numeric(names(SIM_strs)[k])) & ANI_sh_strain[target_str,]>=as.numeric(names(SIM_strs)[k])-1][1]
    }
    
    #get info about target strings
    ANI_SIM_strain<-ANI_all_strain[SIM_strs,SIM_strs]
    #get the hist of share strs
    pdf(file = "ANI_SIM_strain_hist.pdf" , width = 11, height= 8.5)
    hist(ANI_SIM_strain[ANI_SIM_strain<100])
    dev.off()
    write.table(ANI_SIM_strain, file=("ANI_SIM_strain.csv"),sep=",", quote=F, row.names= TRUE, col.names= TRUE, fileEncoding="UTF-8" )
}else {
    if(is.na(as.numeric(range_type))){
        stop(paste("From input, it should be with 'lower_than_range', however the parameter range_type is ",range_type,", rather than correct numeric format",sep=""))
    } else{
    	print(paste("With 'lower_than_range', get the list of strs with ANI lower than the range_type value:",range_type,sep=""))
    	ANI_sh_strain<-ANI_all_strain
    	low_range <- as.numeric(range_type)
      #set all strain to itself is 0
    	for (cn in 1:length(rownames(ANI_sh_strain))) {
          ANI_sh_strain[cn,cn]=0
        }
    	#get the count for strains within the define range 
      sh_ref<-as.character(colnames(ANI_all_strain))
    	hist_list <- matrix(nrow=length(sh_ref),ncol=2)
    	colnames(hist_list)<-c("lower or equal","bigger")
    	rownames(hist_list)<-rownames(ANI_sh_strain)
    	for (row_name in rownames(ANI_sh_strain)) {
    		hist_list[row_name,1] <- length(which(ANI_sh_strain[row_name,]<=low_range))
    		hist_list[row_name,2] <- length(which(ANI_sh_strain[row_name,]>low_range))
    	}
      #find the strains with the most hits and start with them.
    	hist_list_max<-hist_list[which(hist_list[,1]==max(hist_list[,1])),]
    	if(is.null(dim(hist_list_max))){
        SIM_strs<-names(which(hist_list[,1]==max(hist_list[,1])))
        hist_list_sort<-hist_list[order(hist_list[,1],decreasing = TRUE),]
        left_str <- rownames(hist_list_sort)
        add_str_list <- setdiff(left_str,SIM_strs)
        
      }else{
        #get the strains to be added in order
        hist_list_sort<-hist_list[order(hist_list[,1],decreasing = TRUE),]
        left_str <- rownames(hist_list_sort)
        add_str_list <- setdiff(left_str,rownames(hist_list_max))
        seed_matr <- as.matrix(ANI_sh_strain[rownames(hist_list_max),rownames(hist_list_max)])
      
        #get the seed matrix
        name_list <- rownames(seed_matr)
        while(max(seed_matr)>low_range){
          name_max <- which(seed_matr==seed_matr[which.max(seed_matr)],arr.ind=T)
          seed_matr <- seed_matr[-name_max[1],-name_max[1]]
          name_list <- setdiff(name_list,rownames(name_max)[1])
        }
        SIM_strs <- name_list
      }
    	#add one in
    	for (add_row_name in add_str_list) {
    		tm_str<-c(SIM_strs,add_row_name)
            if (max(ANI_sh_strain[tm_str,tm_str])<=low_range){
                SIM_strs<-tm_str
            }
    	}
    	#get info about target strings
        ANI_SIM_strain<-ANI_all_strain[SIM_strs,SIM_strs]
        #get the hist of share strs
        pdf(file = "ANI_ref_strain_hist.pdf" , width = 11, height= 8.5)
        hist(ANI_SIM_strain[ANI_SIM_strain<100])
        dev.off()
        write.table(ANI_SIM_strain, file=("ANI_ref_strain.csv"),sep=",", quote=F, row.names= TRUE, col.names= TRUE, fileEncoding="UTF-8" )
    }
}

ref_dl<-read.table(ref_dl_tb)
ref_dl_list <- matrix(0, nrow=length(ref_dl[,1]),ncol=1)
for(i in 1:length(ref_dl[,1])){
   temp_list <- strsplit(as.character(ref_dl[i,1]), split = "/", fixed = TRUE)[[1]]
   ref_dl_list[i,1] <- strsplit(temp_list[length(temp_list)], split = ".", fixed = TRUE)[[1]][1]
}
sel_ref_dl <- ref_dl[ref_dl_list %in% SIM_strs,1]
write.table(sel_ref_dl, file=paste(ref_dl_tb,"_sel.txt",sep=""),sep=",", quote=F, row.names= FALSE, col.names= FALSE, fileEncoding="UTF-8" )

gff_dl<-read.table(gff_dl_tb)
gff_dl_list <- matrix(0, nrow=length(gff_dl[,1]),ncol=1)
for(i in 1:length(gff_dl[,1])){
   temp_list <- strsplit(as.character(gff_dl[i,1]), split = "/", fixed = TRUE)[[1]]
   gff_dl_list[i,1] <- strsplit(temp_list[length(temp_list)], split = ".", fixed = TRUE)[[1]][1]
}
sel_gff_dl <- gff_dl[gff_dl_list %in% SIM_strs,1]
write.table(sel_gff_dl, file=paste(gff_dl_tb,"_sel.txt",sep=""),sep=",", quote=F, row.names= FALSE, col.names= FALSE, fileEncoding="UTF-8" )


if (length(SIM_strs)<2) {
  stop(paste("\n The range_type：",range_type," is too small to generate enough candidate strs. \n\n",sep=""))
}

