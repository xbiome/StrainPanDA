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

#create log:202010403
#Objective: To generate the final figure in the for figures like fig2B
#Use output from cross all_group_MCC_paired_Fig2B.R on JSD or MCC (in fact, same as other matrix), generate boxplot for str groups. Each dataset has paired two boxplot
#generate one figure at a time, the order of dataset is controlled by the order of input matrix from the cross evaluation step.
#Note1: to recognize the group name, use the second set of - in the given dataset name of the cross evaluation step. 
#Note2: it is hard coded to replace XT and EST in the anno name, other strings will caused the column facet fail.
#developed from the previous constrain_type_plot-jinhui.R improved plot details by Jinhui

#usage example: Rscript /home/yxtan/StrainSIM/constrain_type_plot.R MCC_all_list_f=MCC_all_list_f 
#the output will be the final plots and matrixes
#library(vegan) library(ggplot2) library(phyloseq) library(reshape2)  library(stringr) library(ggpubr) are requried
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

#check whether the file in para is exist
if (!exists("MCC_all_list_f")) {
    stop(paste("\nRscript all_group_MCC_paired_Fig2B.R MCC_all_list_f=",MCC_all_list_f," \nWarning: Usage: The MCC file para is not given, please check the input para. \n\n",sep=""))
}

if (!exists("dist_type")) {
    stop(paste("\nRscript all_group_MCC_paired_Fig2B.R dist_type=",dist_type," \nWarning: Usage: The dist_type para is not given, please check the input para. \n\n",sep=""))
}

#check whether the file in is exist
if (!file.exists(MCC_all_list_f)) {
    stop(paste("\nRscript all_group_MCC_paired_Fig2B.R MCC_all_list_f=",MCC_all_list_f," \nWarning: Usage: The MCC file list file is not exist, please check the path. \n\n",sep=""))
}

auto_groups <- function(data_in){
    library(reshape2)
    library(stringr)
    data_ml <- melt(data_in)
    
    data_ml$treatment=sapply(str_split(as.character(data_ml$variable),"\\.|-"),"[",2)
    
    data_ml$group=str_replace_all(as.character(data_ml$variable),".StrainEst|.StrainPanDA|.PStrain","")

    return(data_ml)
}

boxplot_meta= function(data_in,out_name,dist_type,color_set){
    library(ggplot2)
    library(ggpubr)

    data_ml <- auto_groups(data_in)
    
    p <- ggboxplot(data_ml,x="variable",y="value", color="treatment", add = "jitter",palette = color_set)+ 
    facet_grid(rows=vars(data_ml$nstrs),cols=vars(data_ml$group),scales = "free")+
    xlab("") +
    ylab(dist_type) +
    labs(title=dist_type,colour ="Groups")+
        theme( panel.grid.minor = element_blank(),
               panel.background = element_blank(),panel.border = element_blank(),
               axis.line = element_line(colour = "black",linetype="solid",size = 1),
               panel.grid.major=element_line(colour=NA),axis.text.x=element_blank(), #axis.text.y=element_text(size = 10, face= "bold",family = "",angle = 90),
               legend.text = element_text(size = 15),legend.title = element_text(size=15),axis.title=element_text(size=15))
    
     ggsave(file=paste(out_name,"_boxplot.pdf",sep = "")) 

    p <- ggboxplot(data_ml,x="variable",y="value", color="treatment", add = "jitter",palette = color_set)+ 
    coord_cartesian(ylim = c(0, 1)) + 
    facet_grid(rows=vars(data_ml$nstrs),cols=vars(data_ml$group),scales = "free")+
    xlab("") +
    ylab(dist_type) +
    labs(title=dist_type,colour ="Groups")+
        theme( panel.grid.minor = element_line(colour=NA),
               panel.background = element_blank(),panel.border = element_blank(),
               axis.line = element_line(colour = "black",linetype="solid",size = 1),
               panel.grid.major=element_line(colour=NA),axis.text.x=element_blank(), #axis.text.y=element_text(size = 10, face= "bold",family = "",angle = 90),
               legend.text = element_text(size = 15),legend.title = element_text(size=15),axis.title=element_text(size=15))
        
     ggsave(file=paste(out_name,"_y01_boxplot.pdf",sep = "")) 

     ylim_min <- ceiling(min(data_ml$value)*100/5-1)*5/100
     ylim_max <- ceiling(max(data_ml$value)*100/5+1)*5/100
     p <- ggboxplot(data_ml,x="variable",y="value", color="treatment", add = "jitter",palette = color_set)+ 
    coord_cartesian(ylim = c(ylim_min, ylim_max)) + 
    facet_grid(rows=vars(data_ml$nstrs),cols=vars(data_ml$group),scales = "free")+
    xlab("") +
    ylab(dist_type) +
    labs(title=dist_type,colour ="Groups")+
        theme( panel.grid.minor = element_line(colour=NA),
               panel.background = element_blank(),panel.border = element_blank(),
               axis.line = element_line(colour = "black",linetype="solid",size = 1),
               panel.grid.major=element_line(colour=NA),axis.text.x=element_blank(), #axis.text.y=element_text(size = 10, face= "bold",family = "",angle = 90),
               legend.text = element_text(size = 15),legend.title = element_text(size=15),axis.title=element_text(size=15))
        
     ggsave(file=paste(out_name,"_yauto_boxplot.pdf",sep = ""))
}
violinplot_meta= function(data_in,out_name,dist_type){
    library(ggplot2)
    library(ggpubr)

    data_ml <- auto_groups(data_in)    
        
    p <- ggviolin(data_ml, x="variable", y="value", fill = "treatment", add = "boxplot", add.params = list(fill="white"))+
    facet_grid(rows=vars(data_ml$nstrs),cols=vars(data_ml$group),scales = "free")+
    xlab("") +
    ylab(dist_type) +
    labs(title=dist_type,colour = "Groups")+
        theme( panel.grid.minor = element_blank(),
               panel.background = element_blank(),panel.border = element_blank(),
               axis.line = element_line(colour = "black",linetype="solid",size = 1),
               panel.grid.major=element_line(colour=NA),axis.text.x=element_blank(),
               legend.text = element_text(size = 15),legend.title = element_text(size=15),axis.title=element_text(size=15))
    
    
    ggsave(file=paste(out_name,"_violin.pdf",sep = "")) 
    
}

#Step 1: read in the stat all matrixs and do cross method comparison
plot_para= function(dis_all_list_f,dist_type,mypal){
    dis_all_list <- read.table(dis_all_list_f,as.is=TRUE,stringsAsFactors = FALSE)
    for(i in 1:dim(dis_all_list)[1]){
        stat_all_tmp<-read.table(dis_all_list[i,1],header=TRUE,as.is=TRUE,stringsAsFactors = FALSE)
        assign(dis_all_list[i,2],stat_all_tmp)
        if(i==1){
            #sample_names<-colnames(get(dis_all_list[i,2]))
            JSD_melt <- get(dis_all_list[i,2])
            JSD_melt <- cbind(JSD_melt, rep(dis_all_list[i,2],dim(JSD_melt)[1]))
            colnames(JSD_melt)[dim(JSD_melt)[2]] <- "nstrs"
        }else {
            #sample_names<-union(sample_names,colnames(get(dis_all_list[i,2])))
            JSD_melt_tmp <- get(dis_all_list[i,2])
            JSD_melt_tmp <- cbind(JSD_melt_tmp, rep(dis_all_list[i,2],dim(JSD_melt_tmp)[1]))
            colnames(JSD_melt_tmp)[dim(JSD_melt_tmp)[2]] <- "nstrs"
            JSD_melt <- rbind(JSD_melt,JSD_melt_tmp)
        }
    }
    write.table(JSD_melt,file=paste("All_groups","_strain_",dist_type,"_matr.txt",sep = ""), sep='\t',quote=FALSE,row.names=FALSE, fileEncoding="UTF-8")

    boxplot_meta(JSD_melt,paste("All_groups","_",dist_type,sep = ""),dist_type,mypal)
    violinplot_meta(JSD_melt,paste("All_groups","_",dist_type,sep = ""),dist_type)
    
}


library(ggsci)
mypal =pal_npg("nrc", alpha =0.7)(4)

#MCC_all_list_f <- "MCC_all_file_list.txt"
plot_para(MCC_all_list_f,dist_type,mypal)



