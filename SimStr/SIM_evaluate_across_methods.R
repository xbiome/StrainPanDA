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


#Inputs: the converted output of methods, benchtable, metadate of samples; distance matrixs; statistical matrixs
#Note: the 2nd column in meta must be character and not numeric

#update log:20200928
#add JSD plotting
#add abundance plot using the same structure of ranked abudance plot

#update log:20201017
#correct the ggplot boxplot colour code to make it easier to read ( color = "variable"）
#Change abundance plot to report only the bench strains and all other FPs will be grouped. This will make the plot easier to read. This is done by add the section: keep only the bench strains and merge all the FPs into all_FPs
#In this case the ranked_abund sections are not that meaningful to use.

#update log:20201119
#Add the file path check before running. (Original, it had only the parameter check)

#usage example: Rscript /home/yxtan/StrainSIM/SIM_evaluate_across_methods.R stat_all_list_f=stat_all_list_f metadata_f=metadata_f dis_all_list_f=dis_all_list_f RAD_all_list_f=RAD_all_list_f RA_all_list_f=RA_all_list_f out_name=WGS_Xtrain_EST
#The Xtrain output wes mereged and annotated by gene profile neighbors.
#the output will be the distance of each sample and plots
#library(vegan) library(ggplot2) library(phyloseq) library(reshape2) library(Rmisc) are requried
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

#check whether the para is exist
if (!exists("stat_all_list_f")) {
    stop(paste("\nRscript SIM_evaluate_across_methods.R stat_all_list_f=",stat_all_list_f," \nWarning: Usage: This parameter was not given. \n\n",sep=""))
}

if (!exists("metadata_f")) {
    stop(paste("\nRscript SIM_evaluate_across_methods.R metadata_f=",metadata_f," \nWarning: Usage: This parameter was not given. \n\n",sep=""))
}

if (!exists("dis_all_list_f")) {
    stop(paste("\nRscript SIM_evaluate_across_methods.R dis_all_list_f=",dis_all_list_f," \nWarning: Usage: This parameter was not given. \n\n",sep=""))
}

if (!exists("RAD_all_list_f")) {
    stop(paste("\nRscript SIM_evaluate_across_methods.R RAD_all_list_f=",RAD_all_list_f," \nWarning: Usage: This parameter was not given. \n\n",sep=""))
}

if (!exists("RA_all_list_f")) {
    stop(paste("\nRscript SIM_evaluate_across_methods.R RA_all_list_f=",RA_all_list_f," \nWarning: Usage: This parameter was not given. \n\n",sep=""))
}

if (!exists("AD_all_list_f")) {
    stop(paste("\nRscript SIM_evaluate_across_methods.R AD_all_list_f=",AD_all_list_f," \nWarning: Usage: This parameter was not given. \n\n",sep=""))
}

if (!exists("Ab_all_list_f")) {
    stop(paste("\nRscript SIM_evaluate_across_methods.R Ab_all_list_f=",Ab_all_list_f," \nWarning: Usage: This parameter was not given. \n\n",sep=""))
}

#check whether the file is exist
if (!file.exists(stat_all_list_f)) {
    stop(paste("\nRscript SIM_evaluate_across_methods.R stat_all_list_f=",stat_all_list_f," \nWarning: Usage: simulation bench file is not exist, please check the path. \n\n",sep=""))
}

if (!file.exists(metadata_f)) {
    stop(paste("\nRscript SIM_evaluate_across_methods.R metadata_f=",metadata_f," \nWarning: Usage: Xtrain annotated output file is not exist, please check the path. \n\n",sep=""))
}

if (!file.exists(dis_all_list_f)) {
    stop(paste("\nRscript SIM_evaluate_across_methods.R dis_all_list_f=",dis_all_list_f," \nWarning: Usage: Xtrain merged unannotated output file is not exist, please check the path. \n\n",sep=""))
}

if (!file.exists(RAD_all_list_f)) {
    stop(paste("\nRscript SIM_evaluate_across_methods.R RAD_all_list_f=",RAD_all_list_f," \nWarning: Usage: the distance file of all strains is not exist, please check the path. \n\n",sep=""))
}

if (!file.exists(RA_all_list_f)) {
    stop(paste("\nRscript SIM_evaluate_across_methods.R RA_all_list_f=",RA_all_list_f," \nWarning: Usage: The tree file of all strains is not exist, please check the path. \n\n",sep=""))
}

if (!file.exists(AD_all_list_f)) {
    stop(paste("\nRscript SIM_evaluate_across_methods.R AD_all_list_f=",AD_all_list_f," \nWarning: Usage: the distance file of all strains is not exist, please check the path. \n\n",sep=""))
}

if (!file.exists(Ab_all_list_f)) {
    stop(paste("\nRscript SIM_evaluate_across_methods.R Ab_all_list_f=",Ab_all_list_f," \nWarning: Usage: The tree file of all strains is not exist, please check the path. \n\n",sep=""))
}

boxplot_meta= function(data_in,out_name){
    library(reshape2)
    library(ggplot2)
    data_ml <- melt(data_in)
    #p <- ggplot(data_ml,aes(V2,value, color=V2,fill=variable)) +
    p <- ggplot(data_ml,aes(V2,value, color=variable)) +
    geom_boxplot(alpha=0.4,outlier.shape= "diamond", outlier.colour= "black") + #outliers are the one over 1.58 IQR/sqrt
    #geom_boxplot(alpha=0.4) + geom_point(alpha=0.4) +
    #facet_grid( ~ V2) +
    xlab("") +
    ylab(out_name) +
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour="black"))+ #remove bg and grid
    labs(title=out_name,colour ="Methods")
    #labs(title=out_name,colour = "Groups", fill="Methods")
    ggsave(file=paste(out_name,"_boxplot.pdf",sep = "")) 

}
violinplot_meta= function(data_in,out_name){
    library(reshape2)
    library(ggplot2)
    data_ml <- melt(data_in)
    p <- ggplot(data_ml,aes(V2,value, color=V2,fill=variable)) +
    geom_violin(trim=TRUE,color="white") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(以下设为背景为白色，其实表示不要轮廓线)
    #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
    geom_boxplot(width=0.2,position=position_dodge(0.9))+ #绘制箱线图
    theme_bw()+ #背景变为白色
    #geom_boxplot(alpha=0.4,outlier.shape= "diamond", outlier.colour= "black") + geom_point(alpha=0.4) + #outliers are the one over 1.58 IQR/sqrt
    #geom_boxplot(alpha=0.4) + geom_point(alpha=0.4) +
    xlab("") +
    ylab(out_name) +
    theme_set(theme_bw()) + #white background
    theme(panel.grid.major=element_line(colour=NA)) + #remove grid
    labs(title=out_name,colour = "Groups", fill="Methods")
    ggsave(file=paste(out_name,"_violin.pdf",sep = "")) 
}

#general function to compute paired p-value by using wilcoxon Ab_all_ml_list_dc
pvalue_meta= function(metadata,eudisk_RAD_mg,RAD_all_list,out_name,att_name){
    methods <- RAD_all_list[,2]
    if(length(methods)>1){
        for(uniq_group in unique(metadata[,2])){
            tmp_data <- eudisk_RAD_mg[eudisk_RAD_mg$V2==uniq_group,]
            methods <- RAD_all_list[,2]
            combn_list <- combn(length(methods),2)
            outp_matrix <- matrix("NA",nrow=length(methods)-1,ncol=length(methods)-1)
            rownames(outp_matrix) <- methods[1:(length(methods)-1)]
            colnames(outp_matrix) <- methods[2:length(methods)]
            for(i in 1:dim(combn_list)[2]){
                tem_vec <- cbind(c(tmp_data[,methods[combn_list[1,i]]],tmp_data[,methods[combn_list[2,i]]]),c(rep(combn_list[1,i],dim(tmp_data)[1]),rep(combn_list[2,i],dim(tmp_data)[1])))
                colnames(tem_vec) <- c("am","bm")
                outp_matrix[methods[combn_list[1,i]],methods[combn_list[2,i]]]<-wilcox.test(am~bm,data = tem_vec)$p.value
            }
            write.table(outp_matrix,file=paste(out_name,"_",uniq_group,att_name,sep = ""), sep='\t',quote=FALSE,row.names=TRUE, fileEncoding="UTF-8")
        }
    }
}
#Step 1: read in the stat all matrixs and do cross method comparison
#stat_all_list_f <- "stat_all_file_list.txt"
stat_all_list <- read.table(stat_all_list_f,as.is=TRUE,stringsAsFactors = FALSE)
for(i in 1:dim(stat_all_list)[1]){
    stat_all_tmp<-read.csv(stat_all_list[i,1],as.is=TRUE,stringsAsFactors = FALSE)
    assign(stat_all_list[i,2],stat_all_tmp)
    if(i==1){
        sample_names<-colnames(get(stat_all_list[i,2]))
    }else {
        sample_names<-union(sample_names,colnames(get(stat_all_list[i,2])))
    }
}
#get(stat_all_list[i,2])
sen_matr <- matrix(0, nrow = length(sample_names), ncol = dim(stat_all_list)[1])
precision_matr <- matrix(0, nrow = length(sample_names), ncol = dim(stat_all_list)[1])
f1_matr <- matrix(0, nrow = length(sample_names), ncol = dim(stat_all_list)[1])
MCC_matr <- matrix(0, nrow = length(sample_names), ncol = dim(stat_all_list)[1])
rownames(sen_matr) <- sample_names
colnames(sen_matr) <- stat_all_list[,2]
rownames(precision_matr) <- sample_names
colnames(precision_matr) <- stat_all_list[,2]
rownames(f1_matr) <- sample_names
colnames(f1_matr) <- stat_all_list[,2]
rownames(MCC_matr) <- sample_names
colnames(MCC_matr) <- stat_all_list[,2]

for(i in 1:dim(stat_all_list)[1]){
    sen_matr[colnames(get(stat_all_list[i,2])),i]<-unlist(get(stat_all_list[i,2])["Sensitivity",colnames(get(stat_all_list[i,2]))])
    precision_matr[colnames(get(stat_all_list[i,2])),i]<-unlist(get(stat_all_list[i,2])["Precision",colnames(get(stat_all_list[i,2]))])
    f1_matr[colnames(get(stat_all_list[i,2])),i]<-unlist(get(stat_all_list[i,2])["F-measure",colnames(get(stat_all_list[i,2]))])
    MCC_matr[colnames(get(stat_all_list[i,2])),i]<-unlist(get(stat_all_list[i,2])["MCC",colnames(get(stat_all_list[i,2]))])
}

#control the sample to be plotted in metadata
metadata <- read.table(metadata_f,as.is=TRUE,stringsAsFactors = FALSE)
sen_mg = merge(sen_matr,metadata,by.x="row.names",by.y="V1",all.y= T)
precision_mg = merge(precision_matr,metadata,by.x="row.names",by.y="V1",all.y= T)
f1_mg = merge(f1_matr,metadata,by.x="row.names",by.y="V1",all.y= T)
MCC_mg = merge(MCC_matr,metadata,by.x="row.names",by.y="V1",all.y= T)

write.table(sen_mg,file=paste(out_name,"_Sensitivity_matr.txt",sep = ""), sep='\t',quote=FALSE,row.names=FALSE, fileEncoding="UTF-8")
write.table(precision_mg,file=paste(out_name,"_Precision_matr.txt",sep = ""), sep='\t',quote=FALSE,row.names=FALSE, fileEncoding="UTF-8")
write.table(f1_mg,file=paste(out_name,"_F-measure_matr.txt",sep = ""), sep='\t',quote=FALSE,row.names=FALSE, fileEncoding="UTF-8")
write.table(MCC_mg,file=paste(out_name,"_MCC_matr.txt",sep = ""), sep='\t',quote=FALSE,row.names=FALSE, fileEncoding="UTF-8")

boxplot_meta(sen_mg,paste(out_name,"_Sensitivity",sep = ""))
boxplot_meta(precision_mg,paste(out_name,"_Precision",sep = ""))
boxplot_meta(f1_mg,paste(out_name,"_F-measure",sep = ""))
boxplot_meta(MCC_mg,paste(out_name,"_MCC",sep = ""))
violinplot_meta(sen_mg,paste(out_name,"_Sensitivity",sep = ""))
violinplot_meta(precision_mg,paste(out_name,"_Precision",sep = ""))
violinplot_meta(f1_mg,paste(out_name,"_F-measure",sep = ""))
violinplot_meta(MCC_mg,paste(out_name,"_MCC",sep = ""))

pvalue_meta(metadata,sen_mg,stat_all_list,out_name,"_Sensitivity_pvalue.txt")
pvalue_meta(metadata,precision_mg,stat_all_list,out_name,"_Precision_pvalue.txt")
pvalue_meta(metadata,f1_mg,stat_all_list,out_name,"_F-measure_pvalue.txt")
pvalue_meta(metadata,MCC_mg,stat_all_list,out_name,"_MCC_pvalue.txt")


#Step2: similar to step1 for distance value 
dis_all_list <- read.table(dis_all_list_f,as.is=TRUE,stringsAsFactors = FALSE)
for(i in 1:dim(dis_all_list)[1]){
    stat_all_tmp<-read.csv(dis_all_list[i,1],as.is=TRUE,stringsAsFactors = FALSE)
    assign(dis_all_list[i,2],stat_all_tmp)
    if(i==1){
        sample_names<-colnames(get(dis_all_list[i,2]))
    }else {
        sample_names<-union(sample_names,colnames(get(dis_all_list[i,2])))
    }
}
#get(dis_all_list[i,2])
UnW_matr <- matrix(0, nrow = length(sample_names), ncol = dim(dis_all_list)[1])
Weight_matr <- matrix(0, nrow = length(sample_names), ncol = dim(dis_all_list)[1])
BrayD_matr <- matrix(0, nrow = length(sample_names), ncol = dim(dis_all_list)[1])
JCD_matr <- matrix(0, nrow = length(sample_names), ncol = dim(dis_all_list)[1])
JSD_matr <- matrix(0, nrow = length(sample_names), ncol = dim(dis_all_list)[1])
rownames(UnW_matr) <- sample_names
colnames(UnW_matr) <- dis_all_list[,2]
rownames(Weight_matr) <- sample_names
colnames(Weight_matr) <- dis_all_list[,2]
rownames(BrayD_matr) <- sample_names
colnames(BrayD_matr) <- dis_all_list[,2]
rownames(JCD_matr) <- sample_names
colnames(JCD_matr) <- dis_all_list[,2]
rownames(JSD_matr) <- sample_names
colnames(JSD_matr) <- dis_all_list[,2]

for(i in 1:dim(dis_all_list)[1]){
    UnW_matr[colnames(get(dis_all_list[i,2])),i]<-unlist(get(dis_all_list[i,2])["UnW_D",colnames(get(dis_all_list[i,2]))])
    Weight_matr[colnames(get(dis_all_list[i,2])),i]<-unlist(get(dis_all_list[i,2])["Weight_D",colnames(get(dis_all_list[i,2]))])
    BrayD_matr[colnames(get(dis_all_list[i,2])),i]<-unlist(get(dis_all_list[i,2])["Bray_D",colnames(get(dis_all_list[i,2]))])
    JCD_matr[colnames(get(dis_all_list[i,2])),i]<-unlist(get(dis_all_list[i,2])["JC_D",colnames(get(dis_all_list[i,2]))])
    JSD_matr[colnames(get(dis_all_list[i,2])),i]<-unlist(get(dis_all_list[i,2])["JSD",colnames(get(dis_all_list[i,2]))])
}

UnW_mg = merge(UnW_matr,metadata,by.x="row.names",by.y="V1",all.y= T)
Weight_mg = merge(Weight_matr,metadata,by.x="row.names",by.y="V1",all.y= T)
BrayD_mg = merge(BrayD_matr,metadata,by.x="row.names",by.y="V1",all.y= T)
JCD_mg = merge(JCD_matr,metadata,by.x="row.names",by.y="V1",all.y= T)
JSD_mg = merge(JSD_matr,metadata,by.x="row.names",by.y="V1",all.y= T)

write.table(UnW_mg,file=paste(out_name,"_UnW_UniFrac_matr.txt",sep = ""), sep='\t',quote=FALSE,row.names=FALSE, fileEncoding="UTF-8")
write.table(Weight_mg,file=paste(out_name,"_Weighted_UniFrac_matr.txt",sep = ""), sep='\t',quote=FALSE,row.names=FALSE, fileEncoding="UTF-8")
write.table(BrayD_mg,file=paste(out_name,"_Bray_curtis_matr.txt",sep = ""), sep='\t',quote=FALSE,row.names=FALSE, fileEncoding="UTF-8")
write.table(JCD_mg,file=paste(out_name,"_Jaccard_matr.txt",sep = ""), sep='\t',quote=FALSE,row.names=FALSE, fileEncoding="UTF-8")
write.table(JSD_mg,file=paste(out_name,"_JSD_matr.txt",sep = ""), sep='\t',quote=FALSE,row.names=FALSE, fileEncoding="UTF-8")

boxplot_meta(UnW_mg,paste(out_name,"_UnW_UniFrac",sep = ""))
boxplot_meta(Weight_mg,paste(out_name,"_Weighted_UniFrac",sep = ""))
boxplot_meta(BrayD_mg,paste(out_name,"_Bray_curtis",sep = ""))
boxplot_meta(JCD_mg,paste(out_name,"_Jaccard",sep = ""))
boxplot_meta(JSD_mg,paste(out_name,"_JSD",sep = ""))
violinplot_meta(UnW_mg,paste(out_name,"_UnW_UniFrac",sep = ""))
violinplot_meta(Weight_mg,paste(out_name,"_Weighted_UniFrac",sep = ""))
violinplot_meta(BrayD_mg,paste(out_name,"_Bray_curtis",sep = ""))
violinplot_meta(JCD_mg,paste(out_name,"_Jaccard",sep = ""))
violinplot_meta(JSD_mg,paste(out_name,"_JSD",sep = ""))

pvalue_meta(metadata,UnW_mg,dis_all_list,out_name,"_UnW_UniFrac_pvalue.txt")
pvalue_meta(metadata,Weight_mg,dis_all_list,out_name,"_Weighted_UniFrac_pvalue.txt")
pvalue_meta(metadata,BrayD_mg,dis_all_list,out_name,"_Bray_curtis_pvalue.txt")
pvalue_meta(metadata,JCD_mg,dis_all_list,out_name,"_Jaccard_pvalue.txt")
pvalue_meta(metadata,JSD_mg,dis_all_list,out_name,"_JSD_pvalue.txt")

#Step3: abundance
#3.1 ranked diff of abundance
RAD_all_list <- read.table(RAD_all_list_f,as.is=TRUE,stringsAsFactors = FALSE)
for(i in 1:dim(RAD_all_list)[1]){
    RAD_all_tmp<-read.csv(RAD_all_list[i,1],as.is=TRUE,stringsAsFactors = FALSE)
    assign(RAD_all_list[i,2],RAD_all_tmp)
    if(i==1){
        sample_names<-colnames(get(RAD_all_list[i,2]))
    }else {
        sample_names<-union(sample_names,colnames(get(RAD_all_list[i,2])))
    }
}
#get(RAD_all_list[i,2])
eu_dist <- function(x){
  eudist <- sqrt(sum(x^2))
  return(eudist)
}
eudisk_RAD <- matrix(0, nrow = length(sample_names), ncol = dim(RAD_all_list)[1])
rownames(eudisk_RAD) <- sample_names
colnames(eudisk_RAD) <- RAD_all_list[,2] 
for(i in 1:dim(RAD_all_list)[1]){        
    eudisk_RAD[colnames(get(RAD_all_list[i,2])),i]<-apply(get(RAD_all_list[i,2])[,colnames(get(RAD_all_list[i,2]))],2,eu_dist)   
}
eudisk_RAD_mg = merge(eudisk_RAD,metadata,by.x="row.names",by.y="V1",all.y= T)
write.table(eudisk_RAD_mg,file=paste(out_name,"_ranked_abun_diff_EUdist.txt",sep = ""), sep='\t',quote=FALSE,row.names=FALSE, fileEncoding="UTF-8")
boxplot_meta(eudisk_RAD_mg,paste(out_name,"_ranked_abun_diff_EUdist",sep = ""))
violinplot_meta(eudisk_RAD_mg,paste(out_name,"_ranked_abun_diff_EUdist",sep = ""))

#general function to compute paired p-value
pvalue_meta(metadata,eudisk_RAD_mg,RAD_all_list,out_name,"_ranked_abun_diff_EUdist_pvalue.txt")

rank_num<-dim(get(RAD_all_list[1,2]))[1]
for(j in 1:rank_num){
    rank_tmp <- matrix(0, nrow = length(sample_names), ncol = dim(RAD_all_list)[1])
    rownames(rank_tmp) <- sample_names
    colnames(rank_tmp) <- RAD_all_list[,2]    
    for(i in 1:dim(RAD_all_list)[1]){
        rank_tmp[colnames(get(RAD_all_list[i,2])),i]<-unlist(abs(get(RAD_all_list[i,2])[paste("Rank",j,sep=""),colnames(get(RAD_all_list[i,2]))]))
    }
    rank_tmp_mg = merge(rank_tmp,metadata,by.x="row.names",by.y="V1",all.y= T)
    rank_tmp_mg$Rank <- paste("Rank",j,sep="")
    library(reshape2)
    
    rank_tmp_ml <- melt(rank_tmp_mg)
    if(j==1){
        RAD_all_ml_list <- rank_tmp_ml
    } else{
        RAD_all_ml_list <- rbind(RAD_all_ml_list,rank_tmp_ml)
    }
}
#in the plot, all the diff of abun on ranks were stat together
write.table(RAD_all_ml_list,file=paste(out_name,"_ranked_abun_diff_melt.txt",sep = ""), sep='\t',quote=FALSE,row.names=FALSE, fileEncoding="UTF-8")
library(ggplot2)
#p <- ggplot(RAD_all_ml_list,aes(V2,value, color=V2,fill=variable)) +
p <- ggplot(RAD_all_ml_list,aes(V2,value, color=variable)) +
geom_boxplot(alpha=0.4,outlier.shape= "diamond", outlier.colour= "black") + #outliers are the one over 1.58 IQR/sqrt
#geom_boxplot(alpha=0.4) + geom_point(alpha=0.4) +
xlab("") +
ylab("Abs_of_Ranked_Abundance_Diff_Compared_to_Bench") +
#labs(title=paste(out_name,"_ranked_abun_diff",sep = ""),colour = "Groups", fill="Methods")
labs(title=paste(out_name,"_ranked_abun_diff",sep = ""),colour ="Methods")
ggsave(file=paste(out_name,"_ranked_abun_diff_boxplot.pdf",sep = "")) 

p <- ggplot(RAD_all_ml_list,aes(V2,value, color=V2,fill=variable)) +
geom_violin(trim=TRUE,color="white") + 
geom_boxplot(width=0.2,position=position_dodge(0.9))+ 
theme_bw()+ 
xlab("") +
ylab("Abs_of_Ranked_Abundance_Diff_Compared_to_Bench") +
labs(title=paste(out_name,"_ranked_abun_diff",sep = ""),colour = "Groups", fill="Methods")
ggsave(file=paste(out_name,"_ranked_abun_diff_violin.pdf",sep = "")) 

#3.2 ranked abundance
#Ecoli98_WGS_all_Xtrain_rank_abun_all.csv
#RA_all_list_f <- "prof_all_file_list.txt"
RA_all_list <- read.table(RA_all_list_f,as.is=TRUE,stringsAsFactors = FALSE)
for(i in 1:dim(RA_all_list)[1]){
    RA_all_tmp<-read.csv(RA_all_list[i,1],as.is=TRUE,stringsAsFactors = FALSE)
    assign(RA_all_list[i,2],RA_all_tmp)
    if(i==1){
        sample_names<-colnames(get(RA_all_list[i,2]))
    }else {
        sample_names<-union(sample_names,colnames(get(RA_all_list[i,2])))
    }
}
#get(RA_all_list[i,2])
rank_num<-dim(get(RA_all_list[1,2]))[1]/2

for(j in 1:rank_num){
    rank_tmp <- matrix(0, nrow = length(sample_names), ncol = (dim(RA_all_list)[1]+1))
    rownames(rank_tmp) <- sample_names
    colnames(rank_tmp) <- c(RA_all_list[,2],"bench")
    for(i in 1:dim(RA_all_list)[1]){
        rank_tmp[colnames(get(RA_all_list[i,2])),i]<-unlist(get(RA_all_list[i,2])[paste("predict_rank",j,sep=""),colnames(get(RA_all_list[i,2]))])
    }
    rank_tmp[colnames(get(RA_all_list[i,2])),(dim(RA_all_list)[1]+1)]<-unlist(get(RA_all_list[i,2])[paste("bench_rank",j,sep=""),colnames(get(RA_all_list[i,2]))])
    rank_tmp_mg = merge(rank_tmp,metadata,by.x="row.names",by.y="V1",all.y= T)
    rank_tmp_mg$Rank <- paste("Rank",j,sep="")
    library(reshape2)    
    rank_tmp_ml <- melt(rank_tmp_mg)
    if(j==1){
        RA_all_ml_list <- rank_tmp_ml
    } else{
        RA_all_ml_list <- rbind(RA_all_ml_list,rank_tmp_ml)
    }
}
write.table(RA_all_ml_list,file=paste(out_name,"_ranked_abun_melt.txt",sep = ""), sep='\t',quote=FALSE,row.names=FALSE, fileEncoding="UTF-8")

#plot each group separately
library(ggplot2)
for (group_name in unique(RA_all_ml_list$V2)){
    
    library(Rmisc)
    sub_RA_all_ml_list <-  summarySE(RA_all_ml_list[RA_all_ml_list$V2==group_name,],measurevar="value", groupvars=c("Rank","variable"))
    sub_RA_all_ml_list_tmp <- sub_RA_all_ml_list[sub_RA_all_ml_list$variable=="bench",]
    rank_keep <- sub_RA_all_ml_list_tmp[sub_RA_all_ml_list[sub_RA_all_ml_list$variable=="bench",]$value>0,]$Rank
    sub_RA_all_ml_list_ranked <- sub_RA_all_ml_list[sub_RA_all_ml_list$Rank %in% rank_keep,]
    ggplot(data=sub_RA_all_ml_list_ranked, mapping=aes(variable,value, color=Rank,fill=Rank)) +
    geom_bar(stat="identity",width=0.5,position='dodge')+
    geom_errorbar(aes(ymin=value-se, ymax=value+se), width = 0.5, position = 'dodge') + #errorbar is SE
    xlab("") +
    ylab("Ranked_Relative_Abundance") +
    labs(title=paste(out_name,"_",group_name,sep = ""), fill="Ranks")+
    theme_minimal()
    ggsave(file=paste(out_name,"_",group_name,"_ranked_abun_barplot.pdf",sep = "")) 
}
RA_all_ml_list_stat<-summarySE(RA_all_ml_list,measurevar="value", groupvars=c("Rank","variable","V2"))
write.table(RA_all_ml_list_stat,file=paste(out_name,"_ranked_abun_stat_grouped.txt",sep = ""), sep='\t',quote=FALSE,row.names=FALSE, fileEncoding="UTF-8")

ggplot(RA_all_ml_list_stat, aes(x = V2, y = value, fill = Rank)) + 
geom_bar(stat="identity",position="stack") +
facet_wrap(~variable) + 
theme_bw() + 
labs(x=NULL, y="relative abundances")+
theme(axis.text.x = element_text(angle=45, vjust=1,hjust=1))
ggsave(file=paste(out_name, "_ranked_abun_stackplot_by_method.pdf", sep=""))

ggplot(RA_all_ml_list_stat, aes(x = variable, y = value, fill = Rank)) + 
geom_bar(stat="identity") +
facet_wrap(~V2) + 
theme_bw() + 
labs(x=NULL, y="relative abundances")+
theme(axis.text.x = element_text(angle=45, vjust=1,hjust=1))
ggsave(file=paste(out_name, "_ranked_abun_stackplot_by_group.pdf", sep=""))

#order the V2(groups) by the most abun str in bench
bench_matr <- subset(RA_all_ml_list_stat, RA_all_ml_list_stat$variable=="bench")
max_value <- max(bench_matr$value)
max_strain <- bench_matr$Rank[which(bench_matr$value==max_value)[1]]
test_matr <- subset(bench_matr, bench_matr$Rank==max_strain)
sample_order <- test_matr[order(test_matr$value, decreasing = TRUE), ]$V2
RA_all_ml_list_stat_ordered <- RA_all_ml_list_stat
RA_all_ml_list_stat_ordered$V2 <- factor(RA_all_ml_list_stat_ordered$V2,levels = sample_order,ordered = TRUE)
ggplot(RA_all_ml_list_stat_ordered, aes(x = V2, y = value, fill = Rank)) + 
geom_bar(stat="identity",position="stack") +
facet_wrap(~variable) + 
theme_bw() + 
labs(x=NULL, y="relative abundances")+
theme(axis.text.x = element_text(angle=45, vjust=1,hjust=1))
ggsave(file=paste(out_name, "_ranked_abun_stackplot_by_method_ordered.pdf", sep=""))


######################################################################################
#Step4: abundance
#4.1 abundance difference
#AD_all_list_f <- "abun_diff_all_file_list.txt"
AD_all_list <- read.table(AD_all_list_f,as.is=TRUE,stringsAsFactors = FALSE)
for(i in 1:dim(AD_all_list)[1]){
    AD_all_tmp<-read.csv(AD_all_list[i,1],as.is=TRUE,stringsAsFactors = FALSE)
    assign(AD_all_list[i,2],AD_all_tmp)
    if(i==1){
        sample_names<-colnames(get(AD_all_list[i,2]))
    }else {
        sample_names<-union(sample_names,colnames(get(AD_all_list[i,2])))
    }
}
#get(AD_all_list[i,2])
eu_dist <- function(x){
  eudist <- sqrt(sum(x^2))
  return(eudist)
}
eudisk_AD <- matrix(0, nrow = length(sample_names), ncol = dim(AD_all_list)[1])
rownames(eudisk_AD) <- sample_names
colnames(eudisk_AD) <- AD_all_list[,2] 
for(i in 1:dim(AD_all_list)[1]){        
    eudisk_AD[colnames(get(AD_all_list[i,2])),i]<-apply(get(AD_all_list[i,2])[,colnames(get(AD_all_list[i,2]))],2,eu_dist)   
}
eudisk_AD_mg = merge(eudisk_AD,metadata,by.x="row.names",by.y="V1",all.y= T)
write.table(eudisk_AD_mg,file=paste(out_name,"_abun_diff_EUdist.txt",sep = ""), sep='\t',quote=FALSE,row.names=FALSE, fileEncoding="UTF-8")
boxplot_meta(eudisk_AD_mg,paste(out_name,"_abun_diff_EUdist",sep = ""))
violinplot_meta(eudisk_AD_mg,paste(out_name,"_abun_diff_EUdist",sep = ""))

#general function to compute paired p-value
pvalue_meta(metadata,eudisk_AD_mg,AD_all_list,out_name,"_abun_diff_EUdist_pvalue.txt")

str_names <-rownames(get(AD_all_list[1,2]))
if(dim(AD_all_list)[1]>1){
    for(i in 2:dim(AD_all_list)[1]){
        str_names <-union(str_names,rownames(get(AD_all_list[i,2])))
    }
}
rank_num<-length(str_names)
for(j in 1:rank_num){
    rank_tmp <- matrix(0, nrow = length(sample_names), ncol = dim(AD_all_list)[1])
    rownames(rank_tmp) <- sample_names
    colnames(rank_tmp) <- AD_all_list[,2]    
    for(i in 1:dim(AD_all_list)[1]){
        if(str_names[j]%in%rownames(get(AD_all_list[i,2]))){
            rank_tmp[colnames(get(AD_all_list[i,2])),i]<-unlist(abs(get(AD_all_list[i,2])[str_names[j],colnames(get(AD_all_list[i,2]))]))
        }
    }
    rank_tmp_mg = merge(rank_tmp,metadata,by.x="row.names",by.y="V1",all.y= T)
    rank_tmp_mg$Rank <- str_names[j]
    library(reshape2)
    
    rank_tmp_ml <- melt(rank_tmp_mg)
    if(j==1){
        AD_all_ml_list <- rank_tmp_ml
    } else{
        AD_all_ml_list <- rbind(AD_all_ml_list,rank_tmp_ml)
    }
}
#in the plot, all the diff of abun on ranks were stat together
write.table(AD_all_ml_list,file=paste(out_name,"_abun_diff_melt.txt",sep = ""), sep='\t',quote=FALSE,row.names=FALSE, fileEncoding="UTF-8")
library(ggplot2)
p <- ggplot(AD_all_ml_list,aes(V2,value, color=variable)) +
geom_boxplot(alpha=0.4,outlier.shape= "diamond", outlier.colour= "black") + #outliers are the one over 1.58 IQR/sqrt
#geom_boxplot(alpha=0.4) + geom_point(alpha=0.4) +
xlab("") +
ylab("Abs_of_Abundance_Diff_Compared_to_Bench") +
labs(title=paste(out_name,"_abun_diff",sep = ""),colour ="Methods")
ggsave(file=paste(out_name,"_abun_diff_boxplot.pdf",sep = "")) 

p <- ggplot(AD_all_ml_list,aes(V2,value, color=V2,fill=variable)) +
geom_violin(trim=TRUE,color="white") + 
geom_boxplot(width=0.2,position=position_dodge(0.9))+ #绘制箱线图
theme_bw()+ 
xlab("") +
ylab("Abs_of_Abundance_Diff_Compared_to_Bench") +
labs(title=paste(out_name,"_abun_diff",sep = ""),colour = "Groups", fill="Methods")
ggsave(file=paste(out_name,"_abun_diff_violin.pdf",sep = "")) 

#4.2 abundance
#Ecoli98_WGS_all_Xtrain_rank_abun_all.csv
Ab_all_list <- read.table(Ab_all_list_f,as.is=TRUE,stringsAsFactors = FALSE)
for(i in 1:dim(Ab_all_list)[1]){
    Ab_all_tmp<-read.csv(Ab_all_list[i,1],as.is=TRUE,stringsAsFactors = FALSE)
    assign(Ab_all_list[i,2],Ab_all_tmp)
    if(i==1){
        sample_names<-colnames(get(Ab_all_list[i,2]))
    }else {
        sample_names<-union(sample_names,colnames(get(Ab_all_list[i,2])))
    }
}

rank_num<-rank_num

for(j in 1:rank_num){
    rank_tmp <- matrix(0, nrow = length(sample_names), ncol = (dim(Ab_all_list)[1]+1))
    rownames(rank_tmp) <- sample_names
    colnames(rank_tmp) <- c(Ab_all_list[,2],"bench")
    for(i in 1:dim(Ab_all_list)[1]){
        if(paste("predict_",str_names[j],sep="")%in%rownames(get(Ab_all_list[i,2]))){
            rank_tmp[colnames(get(Ab_all_list[i,2])),i]<-unlist(abs(get(Ab_all_list[i,2])[paste("predict_",str_names[j],sep=""),colnames(get(Ab_all_list[i,2]))]))
            rank_tmp[colnames(get(Ab_all_list[i,2])),(dim(Ab_all_list)[1]+1)]<-unlist(get(Ab_all_list[i,2])[paste("bench_",str_names[j],sep=""),colnames(get(Ab_all_list[i,2]))])
        }

        #rank_tmp[colnames(get(Ab_all_list[i,2])),i]<-unlist(get(Ab_all_list[i,2])[paste("predict_",str_names,sep=""),colnames(get(Ab_all_list[i,2]))])
    }
    #rank_tmp[colnames(get(Ab_all_list[i,2])),(dim(Ab_all_list)[1]+1)]<-unlist(get(Ab_all_list[i,2])[paste("bench_",str_names,sep=""),colnames(get(Ab_all_list[i,2]))])
    rank_tmp_mg = merge(rank_tmp,metadata,by.x="row.names",by.y="V1",all.y= T)
    rank_tmp_mg$Rank <- str_names[j]
    library(reshape2)    
    rank_tmp_ml <- melt(rank_tmp_mg)
    if(j==1){
        Ab_all_ml_list <- rank_tmp_ml
    } else{
        Ab_all_ml_list <- rbind(Ab_all_ml_list,rank_tmp_ml)
    }
}
write.table(Ab_all_ml_list,file=paste(out_name,"_abun_melt.txt",sep = ""), sep='\t',quote=FALSE,row.names=FALSE, fileEncoding="UTF-8")

#keep only the bench strains and merge all the FPs into all_FPs, using dcast(https://www.jianshu.com/p/1c9ec35bef42)
Ab_all_ml_list_dc = dcast(Ab_all_ml_list,Rank + variable ~ Row.names)
Ab_all_ml_list_dc_bench <- Ab_all_ml_list_dc[Ab_all_ml_list_dc$variable=="bench",]
rownames(Ab_all_ml_list_dc_bench)<-Ab_all_ml_list_dc_bench$Rank
Ab_all_ml_list_dc_bench_an <- Ab_all_ml_list_dc_bench[,3:dim(Ab_all_ml_list_dc_bench)[2]]
bench_names <- names(which(rowSums(Ab_all_ml_list_dc_bench_an)>0))
Ab_all_ml_list_dc_select <- Ab_all_ml_list_dc[Ab_all_ml_list_dc$Rank %in% bench_names,]
other_FPs <- Ab_all_ml_list_dc[1:length(levels(Ab_all_ml_list_dc$variable)),]
other_FPs$Rank <- "all_FPs"
for (var_name in levels(Ab_all_ml_list_dc$variable)){
    other_FPs[other_FPs$variable==var_name,3:dim(Ab_all_ml_list_dc_bench)[2]]<- 1-colSums(Ab_all_ml_list_dc_select[Ab_all_ml_list_dc_select$variable==var_name,3:dim(Ab_all_ml_list_dc_bench)[2]])
}
Ab_all_list_select <- rbind(Ab_all_ml_list_dc_select,other_FPs)
Ab_all_ml_list_select <- melt(Ab_all_list_select)
colnames(Ab_all_ml_list_select)[3] <- "sample" 
Ab_all_ml_list_select_mg <- merge(Ab_all_ml_list_select,metadata,by.x="sample",by.y="V1",all.y= T)

#plot each group separately
library(ggplot2)
for (group_name in unique(Ab_all_ml_list_select_mg$V2)){
    library(Rmisc)
    sub_Ab_all_ml_list <-  summarySE(Ab_all_ml_list_select_mg[Ab_all_ml_list_select_mg$V2==group_name,],measurevar="value", groupvars=c("Rank","variable"))
    sub_Ab_all_ml_list_tmp <- sub_Ab_all_ml_list[sub_Ab_all_ml_list$variable=="bench",]
    rank_keep <- sub_Ab_all_ml_list_tmp[sub_Ab_all_ml_list[sub_Ab_all_ml_list$variable=="bench",]$value>=0,]$Rank
    sub_Ab_all_ml_list_ranked <- sub_Ab_all_ml_list[sub_Ab_all_ml_list$Rank %in% rank_keep,]
    ggplot(data=sub_Ab_all_ml_list_ranked, mapping=aes(variable,value, color=Rank,fill=Rank)) +
    geom_bar(stat="identity",width=0.5,position='dodge')+
    geom_errorbar(aes(ymin=value-se, ymax=value+se), width = 0.5, position = 'dodge') + #errorbar is SE
    xlab("") +
    ylab("Relative_Abundance") +
    labs(title=paste(out_name,"_",group_name,sep = ""), fill="Strains")+
    theme_minimal()
    ggsave(file=paste(out_name,"_",group_name,"_abun_barplot.pdf",sep = "")) 
}
Ab_all_ml_list_stat<-summarySE(Ab_all_ml_list_select_mg,measurevar="value", groupvars=c("Rank","variable","V2"))
write.table(Ab_all_ml_list_stat,file=paste(out_name,"_abun_stat_grouped.txt",sep = ""), sep='\t',quote=FALSE,row.names=FALSE, fileEncoding="UTF-8")

ggplot(Ab_all_ml_list_stat, aes(x = V2, y = value, fill = Rank)) + 
geom_bar(stat="identity",position="stack") +
facet_wrap(~variable) + 
theme_bw() + 
labs(x=NULL, y="relative abundances")+
theme(axis.text.x = element_text(angle=45, vjust=1,hjust=1))
ggsave(file=paste(out_name, "_abun_stackplot_by_method.pdf", sep=""))

ggplot(Ab_all_ml_list_stat, aes(x = variable, y = value, fill = Rank)) + 
geom_bar(stat="identity") +
facet_wrap(~V2) + 
theme_bw() + 
labs(x=NULL, y="relative abundances")+
theme(axis.text.x = element_text(angle=45, vjust=1,hjust=1))
ggsave(file=paste(out_name, "_abun_stackplot_by_group.pdf", sep=""))

#order the V2(groups) by the most abun str in bench
bench_matr <- subset(Ab_all_ml_list_stat, Ab_all_ml_list_stat$variable=="bench")
max_value <- max(bench_matr$value)
max_strain <- bench_matr$Rank[which(bench_matr$value==max_value)[1]]
test_matr <- subset(bench_matr, bench_matr$Rank==max_strain)
sample_order <- test_matr[order(test_matr$value, decreasing = TRUE), ]$V2
Ab_all_ml_list_stat_ordered <- Ab_all_ml_list_stat
Ab_all_ml_list_stat_ordered$V2 <- factor(Ab_all_ml_list_stat_ordered$V2,levels = sample_order,ordered = TRUE)
ggplot(Ab_all_ml_list_stat_ordered, aes(x = V2, y = value, fill = Rank)) + 
geom_bar(stat="identity",position="stack") +
facet_wrap(~variable) + 
theme_bw() + 
labs(x=NULL, y="relative abundances")+
theme(axis.text.x = element_text(angle=45, vjust=1,hjust=1))
ggsave(file=paste(out_name, "_abun_stackplot_by_method_ordered.pdf", sep=""))
