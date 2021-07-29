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

#usage example: Rscript Extract_SIM_strs_ANI_PAN_mash.R range_type=wide(or specific number N to get all strs with ANI lower than N) ANI_tb=Ecoli_target_str_info-fastani-1k.out select_info=Ecoli_target_str_info_sel_col.csv pan_ref="target_database/panphlan_target-species_pangenome.csv" min_ANI=95(not required)
#use the ref ANI to select target strains, jaccard index will be calculated and compare as well.
#the output will be a table of selected strs for SIM
#library(tidyverse) ,library(vegan),library(pheatmap) library(ggplot2),library(ggpubr),library(ecodist)are requried
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
  stop(paste("\nRscript Extract_SIM_strs_ANI_PAN_mash.R ANI_tb=",ANI_tb," \nWarning: This parameter was not given. \n\n",sep=""))
}

if (!exists("select_info")) {
  stop(paste("\nRscript Extract_SIM_strs_ANI_PAN_mash.R select_info=",select_info," \nWarning: This parameter was not given. \n\n",sep=""))
}

if (!exists("pan_ref")) {
  stop(paste("\nRscript Extract_SIM_strs_ANI_PAN_mash.R pan_ref=",pan_ref," \nWarning: Usage: This parameter was not given. \n\n",sep=""))
}

#check whether the file is exist
if (!file.exists(ANI_tb)) {
  stop(paste("\nRscript Extract_SIM_strs_ANI_PAN_mash.R ANI_tb=",ANI_tb," \nWarning: Usage: ANI matrix is not exist, please check the path. \n\n",sep=""))
}

if (!file.exists(select_info)) {
  stop(paste("\nRscript Extract_SIM_strs_ANI_PAN_mash.R select_info=",select_info," \nWarning: Usage: The table list of candidate strs is not exist, please check the path. \n\n",sep=""))
}

if (!file.exists(pan_ref)) {
  stop(paste("\nRscript Extract_SIM_strs_ANI_PAN_mash.R pan_ref=",pan_ref," \nWarning: Usage: pan_ref file of candidate strs is not exist, please check the path. \n\n",sep=""))
}


library(tidyverse)
library(vegan)
library(pheatmap)


#min_ANI is not required, default is 0
if (!exists("min_ANI")) { min_ANI="0" }

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
write.table(ANI_all_strain, file=("ANI_allref_strain.csv"),sep=",", quote=F, row.names= TRUE, col.names= TRUE, fileEncoding="UTF-8" )

#get the SIM strs by ANI
if(range_type=="wide") {
    print("Simulation with 'wide_range', gete the list of strs with the widest ANI range and one str in each 1 ANI break.")
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

    #get the SIM strs
    SIM_strs <- rep(0,sum(hist_list_count[target_str,])+1)
    names(SIM_strs) <-c(0,colnames(hist_list_count)[hist_list_count[target_str,]>0])
    SIM_strs[1] <- target_str
    for (k in 2:(sum(hist_list_count[target_str,])+1)) {
      SIM_strs[k] <- colnames(ANI_sh_strain)[(ANI_sh_strain[target_str,]<as.numeric(names(SIM_strs)[k])) & ANI_sh_strain[target_str,]>=as.numeric(names(SIM_strs)[k])-1][1]
    }
    
    #get info about target strings
    ANI_SIM_strain<-as.matrix(ANI_all_strain[SIM_strs,SIM_strs])

    #check whether the min ANI is qualified
    while(min(ANI_SIM_strain)<min_ANI){
        name_min <- which(ANI_SIM_strain==ANI_SIM_strain[which.min(ANI_SIM_strain)],arr.ind=T)
        ANI_SIM_strain <- ANI_SIM_strain[-name_min[1],-name_min[1]]
    }


    #get the hist of share strs
    pdf(file = "ANI_SIM_strain_hist.pdf" , width = 11, height= 8.5)
    hist(ANI_SIM_strain[ANI_SIM_strain<100])
    dev.off()
    write.table(ANI_SIM_strain, file=("ANI_SIM_strain.csv"),sep=",", quote=F, row.names= TRUE, col.names= TRUE, fileEncoding="UTF-8" )
}else {
    if(is.na(as.numeric(range_type))){
        stop(paste("From input, it should be simulation with 'lower_than_range', however the parameter range_type is ",range_type,", rather than correct numeric format",sep=""))
    } else{
    	print(paste("Simulation with 'lower_than_range', gete the list of strs with ANI lower than the range_type value:",range_type,sep=""))
    	ANI_sh_strain<-ANI_all_strain
    	low_range <- as.numeric(range_type)
    	for (cn in 1:length(rownames(ANI_sh_strain))) {
          ANI_sh_strain[cn,cn]=0
        }
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
        ANI_SIM_strain<-as.matrix(ANI_all_strain[SIM_strs,SIM_strs])

        #check whether the min ANI is qualified
        while(min(ANI_SIM_strain)<min_ANI){
            name_min <- which(ANI_SIM_strain==ANI_SIM_strain[which.min(ANI_SIM_strain)],arr.ind=T)
            ANI_SIM_strain <- ANI_SIM_strain[-name_min[1],-name_min[1]]
        }

        #get the hist of share strs
        pdf(file = "ANI_SIM_strain_hist.pdf" , width = 11, height= 8.5)
        hist(ANI_SIM_strain[ANI_SIM_strain<100])
        dev.off()
        write.table(ANI_SIM_strain, file=("ANI_SIM_strain.csv"),sep=",", quote=F, row.names= TRUE, col.names= TRUE, fileEncoding="UTF-8" )
    }
}

if (length(SIM_strs)<2) {
  stop(paste("\n The range_type：",range_type," is too small to generate enough candidate strs. \n\n",sep=""))
}


#read in the pangenome info
refs_matrix <- read.table(pan_ref) %>% 
  select(1,3) %>% 
  reshape2::acast(V1~V3) 
refs_matrix_NC <- refs_matrix[which(rowSums(refs_matrix==0)>0),] #remove core genes
refs_matrix_NC_binary <- (refs_matrix_NC>0)*1
refs_gene_count <- apply(refs_matrix_NC_binary,2,sum)
ref_stat <- c(max(refs_matrix_NC),dim(refs_matrix)[1],dim(refs_matrix_NC)[1],dim(refs_matrix)[1]-dim(refs_matrix_NC)[1],mean(refs_gene_count),sd(refs_gene_count),max(refs_gene_count),min(refs_gene_count))
names(ref_stat) <- c("Max number of replicons","Total number of pangenes","Number of accessary genes", "Number of core genes","mean of accessary genes per str","sd of accessary genes per str","max of accessary genes per str","min of accessary genes per str")
if (dim(refs_matrix_NC_binary)[1]>60000){
  refs_matrix_NC_binary_plot<-refs_matrix_NC_binary[sample(rownames(refs_matrix_NC_binary), 60000, replace = FALSE),]
  print("WARNING:the number of genes is over the plotting limit, so 60000 rows were randomly chosed for plotting")
} else{
  refs_matrix_NC_binary_plot <-refs_matrix_NC_binary
}
pheatmap(refs_matrix_NC_binary_plot,filename = "pan_ref_accessarygenes.pdf")
write.table(ref_stat, file=("pan_ref_accessarygenes_stat.csv"),sep=",", quote=F, row.names= TRUE, col.names= FALSE, fileEncoding="UTF-8" )

strain_num <- dim(refs_matrix_NC_binary)[2]
refs_matrix_NC_binary_dis <- matrix(nrow=strain_num,ncol=strain_num)
colnames(refs_matrix_NC_binary_dis)<-colnames(refs_matrix_NC_binary)
rownames(refs_matrix_NC_binary_dis)<-colnames(refs_matrix_NC_binary)
for (col_name in colnames(refs_matrix_NC_binary_dis)) {
  for (row_name in rownames(refs_matrix_NC_binary_dis)) {
    refs_matrix_NC_binary_dis[row_name,col_name]<-sum(refs_matrix_NC_binary[,row_name]!=refs_matrix_NC_binary[,col_name])
  }
}
#convert the names into the same format
pan_names <- colnames(refs_matrix_NC_binary_dis)
for(i in 1:length(pan_names)){
  pan_names[i] <- strsplit(pan_names[i], split = ".", fixed = TRUE)[[1]][1]
}
colnames(refs_matrix_NC_binary_dis)<-pan_names
rownames(refs_matrix_NC_binary_dis)<-pan_names
pdf(file = "pan_ref_accessarygenes_diff_hist.pdf" , width = 11, height= 8.5)
hist(unlist(refs_matrix_NC_binary_dis[which(refs_matrix_NC_binary_dis>0)]))
dev.off()

#calculate Jaccard Index (similarity) on all genes
library(ecodist)
#The size of matrix could not exceed 20000*325=6500000; As a result, cut the matrix into equal pieces and calculate means.
if(prod(dim(refs_matrix))>6500000){
    group_n <-as.integer(prod(dim(refs_matrix))/6500000)+1
    all_row_n <- dim(refs_matrix)[1]
    sub_row_n <- round(all_row_n/group_n)
    refs_matrix_JI <- matrix(0,nrow=dim(refs_matrix)[2],ncol=dim(refs_matrix)[2])
    for(g_n in 1:group_n){
        refs_matrix_JI_temp<-1-as.matrix(distance(t(refs_matrix[((g_n-1)*sub_row_n+1):min(g_n*sub_row_n,all_row_n),]), method = "jaccard"))
        refs_matrix_JI=refs_matrix_JI+refs_matrix_JI_temp    
    }
    refs_matrix_JI<-refs_matrix_JI/group_n
    JI_names <- colnames(refs_matrix_JI_temp)
}else{
    refs_matrix_JI<-1-as.matrix(distance(t(refs_matrix), method = "jaccard"))
    JI_names <- colnames(refs_matrix_JI)
}
for(i in 1:length(JI_names)){
  JI_names[i] <- strsplit(JI_names[i], split = ".", fixed = TRUE)[[1]][1]
}
colnames(refs_matrix_JI)<-JI_names
rownames(refs_matrix_JI)<-JI_names
pdf(file = "pan_ref_JI_hist.pdf" , width = 11, height= 8.5)
hist(unlist(refs_matrix_JI[which(refs_matrix_JI<1)]))
dev.off()
write.table(refs_matrix_JI[SIM_strs,SIM_strs], file=("AJI_SIM_strain.csv"),sep=",", quote=F, row.names= TRUE, col.names= TRUE, fileEncoding="UTF-8" )
write.table(refs_matrix_JI, file=("AJI_allref_strain.csv"),sep=",", quote=F, row.names= TRUE, col.names= TRUE, fileEncoding="UTF-8" )


#calculate Jaccard Index (similarity) on accessary genes
if(prod(dim(refs_matrix_NC_binary))>6500000){
    group_n <-as.integer(prod(dim(refs_matrix_NC_binary))/6500000)+1
    all_row_n <- dim(refs_matrix_NC_binary)[1]
    sub_row_n <- round(all_row_n/group_n)
    refs_matrix_NC_binary_JI <- matrix(0,nrow=dim(refs_matrix_NC_binary)[2],ncol=dim(refs_matrix_NC_binary)[2])
    for(g_n in 1:group_n){
        refs_matrix_NC_binary_JI_temp<-1-as.matrix(distance(t(refs_matrix_NC_binary[((g_n-1)*sub_row_n+1):min(g_n*sub_row_n,all_row_n),]), method = "jaccard"))
        refs_matrix_NC_binary_JI=refs_matrix_JI+refs_matrix_JI_temp    
    }
    refs_matrix_NC_binary_JI<-refs_matrix_NC_binary_JI/group_n
    JI_names <- colnames(refs_matrix_NC_binary_JI_temp)
}else{
    refs_matrix_NC_binary_JI<-1-as.matrix(distance(t(refs_matrix_NC_binary), method = "jaccard"))
    JI_names <- colnames(refs_matrix_NC_binary_JI)
}
for(i in 1:length(JI_names)){
  JI_names[i] <- strsplit(JI_names[i], split = ".", fixed = TRUE)[[1]][1]
}
colnames(refs_matrix_NC_binary_JI)<-JI_names
rownames(refs_matrix_NC_binary_JI)<-JI_names
pdf(file = "pan_ref_accessarygenes_JI_hist.pdf" , width = 11, height= 8.5)
hist(unlist(refs_matrix_NC_binary_JI[which(refs_matrix_NC_binary_JI<1)]))
dev.off()

#get the interscet betweet these two sets, in case some ref was lost in downloading.
sh_ref <- intersect(colnames(refs_matrix_NC_binary_dis),as.character(colnames(ANI_all_strain)))

ANI_sh_strain<-ANI_all_strain[sh_ref,sh_ref]
refs_matrix_NC_sh_dis <- refs_matrix_NC_binary_dis[sh_ref,sh_ref]
refs_matrix_NC_sh_JI <- refs_matrix_NC_binary_JI[sh_ref,sh_ref]
refs_matrix_sh_JI <- refs_matrix_JI[sh_ref,sh_ref]

#plot correlation
library(ggplot2)
library(ggpubr)
dim_matrix <- length(sh_ref)
ANI_dis_sh <- matrix(nrow=dim_matrix*dim_matrix,ncol=2)
ANI_dis_sh[,1] <- unlist(ANI_sh_strain)
ANI_dis_sh[,2] <- unlist(refs_matrix_NC_sh_dis)
colnames(ANI_dis_sh)<-c("ANI","gene_diff")
ggplot(data=data.frame(ANI_dis_sh), aes(x=ANI,y=gene_diff))+geom_point(color="red")+stat_smooth(method="lm",se=FALSE)+stat_cor(data=data.frame(ANI_dis_sh), method = "pearson")
ggsave(file="corr_on_dis_ANI_of_share_str.pdf")

ANI_dis_sh <- matrix(nrow=dim_matrix*dim_matrix,ncol=2)
ANI_dis_sh[,1] <- unlist(ANI_sh_strain)
ANI_dis_sh[,2] <- unlist(refs_matrix_sh_JI)
colnames(ANI_dis_sh)<-c("ANI","Jaccard_Index")
ggplot(data=data.frame(ANI_dis_sh), aes(x=ANI,y=Jaccard_Index))+geom_point(color="red")+stat_smooth(method="lm",se=FALSE)+stat_cor(data=data.frame(ANI_dis_sh), method = "pearson")
ggsave(file="corr_on_AJI_ANI_of_share_str.pdf")

ANI_dis_sh <- matrix(nrow=dim_matrix*dim_matrix,ncol=2)
ANI_dis_sh[,1] <- unlist(ANI_sh_strain)
ANI_dis_sh[,2] <- unlist(refs_matrix_NC_sh_JI)
colnames(ANI_dis_sh)<-c("ANI","Jaccard_Index")
ggplot(data=data.frame(ANI_dis_sh), aes(x=ANI,y=Jaccard_Index))+geom_point(color="red")+stat_smooth(method="lm",se=FALSE)+stat_cor(data=data.frame(ANI_dis_sh), method = "pearson")
ggsave(file="corr_on_JI_ANI_of_share_str.pdf")

ANI_dis_sh <- matrix(nrow=dim_matrix*dim_matrix,ncol=2)
ANI_dis_sh[,1] <- unlist(refs_matrix_NC_sh_dis)
ANI_dis_sh[,2] <- unlist(refs_matrix_NC_sh_JI)
colnames(ANI_dis_sh)<-c("gene_diff","Jaccard_Index")
ggplot(data=data.frame(ANI_dis_sh), aes(x=gene_diff,y=Jaccard_Index))+geom_point(color="red")+stat_smooth(method="lm",se=FALSE)+stat_cor(data=data.frame(ANI_dis_sh), method = "pearson")
ggsave(file="corr_on_dis_JI_of_share_str.pdf")

ANI_dis_sh <- matrix(nrow=dim_matrix*dim_matrix,ncol=2)
ANI_dis_sh[,1] <- unlist(refs_matrix_NC_sh_dis)
ANI_dis_sh[,2] <- unlist(refs_matrix_sh_JI)
colnames(ANI_dis_sh)<-c("gene_diff","Jaccard_Index")
ggplot(data=data.frame(ANI_dis_sh), aes(x=gene_diff,y=Jaccard_Index))+geom_point(color="red")+stat_smooth(method="lm",se=FALSE)+stat_cor(data=data.frame(ANI_dis_sh), method = "pearson")
ggsave(file="corr_on_dis_AJI_of_share_str.pdf")


#get WGS dl for SIM strains
SRA_table <- read.csv(select_info)
col_n <- dim(SRA_table)[2]
SRA_table_dedup <- SRA_table[!duplicated(SRA_table[,col_n]),]
rownames(SRA_table_dedup)<-SRA_table_dedup[,col_n]
SRA_table_SIM <- SRA_table_dedup[SIM_strs,]
SRA_table_SIM$Size.Mb. <- round(SRA_table_SIM[,"Size.Mb."]*1012515)
write.table(SRA_table_SIM, file=("SIM_WGS_matrix.csv"),sep=",", quote=F, row.names= FALSE, col.names= TRUE, fileEncoding="UTF-8" )
write.table(SRA_table_SIM[,"download_path"], file=("SIM_WGS_dl_list.txt"),sep=",", quote=F, row.names= FALSE, col.names= FALSE, fileEncoding="UTF-8" )
fastq_dump_list <- as.character(SRA_table_SIM[,"download_path"])
for(i in 1:length(fastq_dump_list)){
  fastq_dump_list[i] <- paste(getwd(),"WGS_data",tail(strsplit(fastq_dump_list[i], split = "/", fixed = TRUE)[[1]],n=1),sep="/")
}
write.table(fastq_dump_list, file=("SIM_WGS_fastdump_list.txt"),sep=",", quote=F, row.names= FALSE, col.names= FALSE, fileEncoding="UTF-8" )
SRA_table_SIM[,"download_path"] <- fastq_dump_list
#write.table(SRA_table_SIM[,c("bac_list_ID","download_path","Size.Mb.")], file=("SIM_WGS_dl_matrix.csv"),sep=",", quote=F, row.names= FALSE, col.names= FALSE, fileEncoding="UTF-8" )

#get str_ref path for SIM strains
SRA_table_SIM$ART<-paste(getwd(),"/strain-ref/",SRA_table_SIM[,"bac_list"],sep="")
write.table(SRA_table_SIM[,"ART"], file=("SIM_ref_unzip.txt"),sep=",", quote=F, row.names= FALSE, col.names= FALSE, fileEncoding="UTF-8" )
write.table(SRA_table_SIM[,c("bac_list_ID","ART")], file=("SIM_ART_matrix.csv"),sep=",", quote=F, row.names= FALSE, col.names= FALSE, fileEncoding="UTF-8" )
write.table(SRA_table_SIM[,c("bac_list_ID","download_path","Size.Mb.")], file=("SIM_WGS_matrix.csv"),sep=",", quote=F, row.names= FALSE, col.names= FALSE, fileEncoding="UTF-8" )
