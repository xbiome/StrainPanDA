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


#created 2021-0604
#modified from SIM_evaluate_across_methods.R, use only the sort abundance and JSD.
#注：meta里面的分组名称不能 用数值，因为会被当成新的一组结果来识别，而不是注释。


#usage example: Rscript /home/tanyuxiang/StrainSIM/SIM_evaluate_across_methods_sJSD.R sJSD_all_list_f=sJSD_all_list_f metadata_f=metadata_f sAb_all_list_f=sAb_all_list_f out_name=WGS_Xtrain_EST
#the output will be the sJSD files and plots with stackplot of sort abundance
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
if (!exists("sJSD_all_list_f")) {
    stop(paste("\nRscript SIM_evaluate_across_methods_sJSD.R sJSD_all_list_f=",sJSD_all_list_f," \nWarning: Usage: This parameter was not given. \n\n",sep=""))
}

if (!exists("sAb_all_list_f")) {
    stop(paste("\nRscript SIM_evaluate_across_methods_sJSD.R sAb_all_list_f=",sAb_all_list_f," \nWarning: Usage: This parameter was not given. \n\n",sep=""))
}

if (!exists("metadata_f")) {
    stop(paste("\nRscript SIM_evaluate_across_methods_sJSD.R metadata_f=",metadata_f," \nWarning: Usage: This parameter was not given. \n\n",sep=""))
}


#check whether the file is exist
if (!file.exists(sJSD_all_list_f)) {
    stop(paste("\nRscript SIM_evaluate_across_methods_sJSD.R sJSD_all_list_f=",sJSD_all_list_f," \nWarning: Usage: sort JSD file is not exist, please check the path. \n\n",sep=""))
}

if (!file.exists(sAb_all_list_f)) {
    stop(paste("\nRscript SIM_evaluate_across_methods_sJSD.R sAb_all_list_f=",sAb_all_list_f," \nWarning: Usage: sort abundance file is not exist, please check the path. \n\n",sep=""))
}


if (!file.exists(metadata_f)) {
    stop(paste("\nRscript SIM_evaluate_across_methods_sJSD.R metadata_f=",metadata_f," \nWarning: Usage: metadata file is not exist, please check the path. \n\n",sep=""))
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
pvalue_meta= function(metadata,eudisk_RAD_mg,RssAb_all_list,out_name,att_name){
    methods <- RssAb_all_list[,2]
    if(length(methods)>1){
        for(uniq_group in unique(metadata[,2])){
            tmp_data <- eudisk_RAD_mg[eudisk_RAD_mg$V2==uniq_group,]
            methods <- RssAb_all_list[,2]
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
#sJSD_all_list_f <- "stat_all_file_list.txt"
sJSD_all_list <- read.table(sJSD_all_list_f,as.is=TRUE,stringsAsFactors = FALSE)
for(i in 1:dim(sJSD_all_list)[1]){
    sJSD_all_tmp<-read.csv(sJSD_all_list[i,1],as.is=TRUE,stringsAsFactors = FALSE)
    assign(sJSD_all_list[i,2],sJSD_all_tmp)
    if(i==1){
        sample_names<-colnames(get(sJSD_all_list[i,2]))
    }else {
        sample_names<-union(sample_names,colnames(get(sJSD_all_list[i,2])))
    }
}
#get(sJSD_all_list[i,2])
sJSD_matr <- matrix(0, nrow = length(sample_names), ncol = dim(sJSD_all_list)[1])
rownames(sJSD_matr) <- sample_names
colnames(sJSD_matr) <- sJSD_all_list[,2]

for(i in 1:dim(sJSD_all_list)[1]){
    sJSD_matr[colnames(get(sJSD_all_list[i,2])),i]<-unlist(get(sJSD_all_list[i,2])["sorted_JSD",colnames(get(sJSD_all_list[i,2]))])
    
}

#use meta data to control the sample to be plot
metadata <- read.table(metadata_f,as.is=TRUE,stringsAsFactors = FALSE)
sJSD_mg = merge(sJSD_matr,metadata,by.x="row.names",by.y="V1",all.y= T)

write.table(sJSD_mg,file=paste(out_name,"_sJSD_matr.txt",sep = ""), sep='\t',quote=FALSE,row.names=FALSE, fileEncoding="UTF-8")

boxplot_meta(sJSD_mg,paste(out_name,"_sJSD",sep = ""))
violinplot_meta(sJSD_mg,paste(out_name,"_sJSD",sep = ""))
pvalue_meta(metadata,sJSD_mg,sJSD_all_list,out_name,"_sJSD_pvalue.txt")



######################################################################################
#Step4: sorted abundance
#sAb_all_list_f <- "sAb_all_file_list.txt"
sAb_all_list <- read.table(sAb_all_list_f,as.is=TRUE,stringsAsFactors = FALSE)
for(i in 1:dim(sAb_all_list)[1]){
    Ab_all_tmp<-read.csv(sAb_all_list[i,1],as.is=TRUE,stringsAsFactors = FALSE)
    assign(sAb_all_list[i,2],Ab_all_tmp)
    if(i==1){
        sample_names<-colnames(get(sAb_all_list[i,2]))
    }else {
        sample_names<-union(sample_names,colnames(get(sAb_all_list[i,2])))
    }
}

str_names <-rownames(get(sAb_all_list[1,2]))
if(dim(sAb_all_list)[1]>1){
    for(i in 2:dim(sAb_all_list)[1]){
        str_names <-union(str_names,rownames(get(sAb_all_list[i,2])))
    }
}

for (str_i in 1:length(str_names)){
    str_names_tmp <- strsplit(str_names[str_i],split='_',fixed=TRUE)[[1]]
    str_names[str_i]<- paste(str_names_tmp[2:length(str_names_tmp)],collapse="_")
}

str_names <- unique(str_names)
rank_num<-length(str_names)

for(j in 1:rank_num){
    rank_tmp <- matrix(0, nrow = length(sample_names), ncol = (dim(sAb_all_list)[1]+1))
    rownames(rank_tmp) <- sample_names
    colnames(rank_tmp) <- c(sAb_all_list[,2],"bench")
    for(i in 1:dim(sAb_all_list)[1]){
        if(paste("predict_",str_names[j],sep="")%in%rownames(get(sAb_all_list[i,2]))){
            rank_tmp[colnames(get(sAb_all_list[i,2])),i]<-unlist(abs(get(sAb_all_list[i,2])[paste("predict_",str_names[j],sep=""),colnames(get(sAb_all_list[i,2]))]))
            rank_tmp[colnames(get(sAb_all_list[i,2])),(dim(sAb_all_list)[1]+1)]<-unlist(get(sAb_all_list[i,2])[paste("bench_",str_names[j],sep=""),colnames(get(sAb_all_list[i,2]))])
        }

    }
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
write.table(Ab_all_ml_list,file=paste(out_name,"_sAB_melt.txt",sep = ""), sep='\t',quote=FALSE,row.names=FALSE, fileEncoding="UTF-8")

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
sAb_all_list_select <- rbind(Ab_all_ml_list_dc_select,other_FPs)
Ab_all_ml_list_select <- melt(sAb_all_list_select)
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
    ggsave(file=paste(out_name,"_",group_name,"_sAB_barplot.pdf",sep = "")) 
}
Ab_all_ml_list_stat<-summarySE(Ab_all_ml_list_select_mg,measurevar="value", groupvars=c("Rank","variable","V2"))
write.table(Ab_all_ml_list_stat,file=paste(out_name,"_sAB_stat_grouped.txt",sep = ""), sep='\t',quote=FALSE,row.names=FALSE, fileEncoding="UTF-8")

ggplot(Ab_all_ml_list_stat, aes(x = V2, y = value, fill = Rank)) + 
geom_bar(stat="identity",position="stack") +
facet_wrap(~variable) + 
theme_bw() + 
labs(x=NULL, y="relative abundances")+
theme(axis.text.x = element_text(angle=45, vjust=1,hjust=1))
ggsave(file=paste(out_name, "_sAB_stackplot_by_method.pdf", sep=""))

ggplot(Ab_all_ml_list_stat, aes(x = variable, y = value, fill = Rank)) + 
geom_bar(stat="identity") +
facet_wrap(~V2) + 
theme_bw() + 
labs(x=NULL, y="relative abundances")+
theme(axis.text.x = element_text(angle=45, vjust=1,hjust=1))
ggsave(file=paste(out_name, "_sAB_stackplot_by_group.pdf", sep=""))

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
ggsave(file=paste(out_name, "_sAB_stackplot_by_method_ordered.pdf", sep=""))
