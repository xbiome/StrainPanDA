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

#create log:20210406
#Objective: To generate the final figure 2A stackplot in the constrain style
#Use output from cross SIM_evaluate_across_methods.R on relative abundance by the abund_metr, generate grided stack plot for str groups.
#Updated from constrain_type_plot-jinhui.R


#usage example: Rscript /home/yxtan/StrainSIM/constrain_type_plot.R all_abun_list_f=all_abun_list_f sample_order_f=sample_order_f
#the output will be the final plots and matrixes
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

#check whether the file in para is exist
if (!exists("all_abun_list_f")) {
    stop(paste("\nRscript SIM_evaluate_across_methods.R all_abun_list_f=",all_abun_list_f," \nWarning: Usage: The relative abundance file para is not given, please check the input para. \n\n",sep=""))
}

if (!exists("sample_order_f")) {
    stop(paste("\nRscript SIM_evaluate_across_methods.R sample_order_f=",sample_order_f," \nWarning: Usage: The sample order file para is not given, please check the input para. \n\n",sep=""))
}

#check whether the file in is exist
if (!file.exists(all_abun_list_f)) {
    stop(paste("\nRscript SIM_evaluate_across_methods.R all_abun_list_f=",all_abun_list_f," \nWarning: Usage: The relative abundance file list file is not exist, please check the path. \n\n",sep=""))
}

if (!file.exists(sample_order_f)) {
    stop(paste("\nRscript SIM_evaluate_across_methods.R sample_order_f=",sample_order_f," \nWarning: Usage: The sample order list file is not exist, please check the path. \n\n",sep=""))
}

###### change
merge_nb_not_bench= function(Ab_all_ml_list){
    library(reshape2)
    library(dplyr)
    Ab_all_ml_list_dc = dcast(Ab_all_ml_list,Rank + variable ~ Row.names)
    Ab_all_ml_list_dc_bench <- Ab_all_ml_list_dc[Ab_all_ml_list_dc$variable=="bench",]
    
    rownames(Ab_all_ml_list_dc_bench)<-Ab_all_ml_list_dc_bench$Rank
    
    Ab_all_ml_list_dc_bench_an <- Ab_all_ml_list_dc_bench[,3:dim(Ab_all_ml_list_dc_bench)[2]]
    
    bench_names <- names(which(rowSums(Ab_all_ml_list_dc_bench_an)>0))
    Ab_all_ml_list_dc_select <- Ab_all_ml_list_dc[Ab_all_ml_list_dc$Rank %in% bench_names,]
    
    #### change
    ### not select data as other FPs, and summary again
    other_FPs <- Ab_all_ml_list_dc[Ab_all_ml_list_dc$Rank %in% setdiff(unique(Ab_all_ml_list_dc$Rank),bench_names),]
        
    
    # other_FPs$variable=rownames(other_FPs)
    other_FPs$Rank1=other_FPs$Rank
    other_FPs$Rank="neighbors_not_bench" 
    
    Ab_all_list_select <- full_join(Ab_all_ml_list_dc_select,other_FPs)
    
    Ab_all_ml_list_select <- melt(Ab_all_list_select)
    colnames(Ab_all_ml_list_select)[4] <- "sample" 
    sample_ID <- as.character(Ab_all_ml_list_select$sample)
    
    for (i in 1:length(sample_ID)){
        sample_ID[i] <- strsplit(sample_ID[i],"_",fixed=TRUE)[[1]][1]
    }
    
    Ab_all_ml_list_stat<-cbind(Ab_all_ml_list_select,sample_ID)
    colnames(Ab_all_ml_list_stat)[6] <- "V2"
    return(Ab_all_ml_list_stat)
}


### change 2021.04.09
scale_str=function(Ab_all_ml_list_stat_ordered,var1,var2){
    df=Ab_all_ml_list_stat_ordered
    for (i in unique(df[,var1])) {
        for (j in unique(df[,var2])) {
          df[df[,var1]==i & df[,var2]==j,"scale_value"]=df[df[,var1]==i & df[,var2]==j,"value"]/sum(df[df[,var1]==i & df[,var2]==j,"value"])
        }
    }
    return(df)
}

others_plot=function(Ab_all_ml_list,out_name, sample_order){
    library(ggplot2)
    library(stringr)
    library(ggpubr)
    Ab_all_ml_list_stat <- merge_nb_not_bench(Ab_all_ml_list)
    
    ## for 2strains
    Ab_all_ml_list_stat_ordered <- Ab_all_ml_list_stat
    
    ### change
    Ab_all_ml_list_stat_ordered$variable=str_replace_all(Ab_all_ml_list_stat_ordered$variable,"bench","Ground Truth")
    
    
    Ab_all_ml_list_stat_ordered$variable=factor(Ab_all_ml_list_stat_ordered$variable,
                                                levels = sample_order,ordered = TRUE)
    
    
    Ab_all_ml_list_stat_ordered$Rank=factor(Ab_all_ml_list_stat_ordered$Rank,
                                            levels=unique(Ab_all_ml_list_stat_ordered$Rank),
                                            labels=c(paste("Strain",seq(1,length(unique(Ab_all_ml_list_stat_ordered$Rank))-1),sep = ""),"Extras"))
    
    
    Ab_all_ml_list_stat_ordered=na.omit(Ab_all_ml_list_stat_ordered) 
    
    Ab_all_ml_list_stat_ordered=scale_str(Ab_all_ml_list_stat_ordered,
                                          var1 = "variable",var2="V2") %>% filter(scale_value>0)    
    
    
    ggbarplot(Ab_all_ml_list_stat_ordered,x="V2",y="scale_value",fill = "Rank1")+
        facet_wrap(~variable,nrow=1)+ 
        # guides(fill="none")+
        labs(x=NULL, y="Relative Abundance",fill="Strain")+
        theme( axis.text.x=element_blank(),
               legend.text = element_text(size = 15),legend.title = element_text(size=15),axis.title=element_text(size=15))
    
    

    ggsave(file=paste(out_name, "_abun_stackplot_by_method_ordered_others.pdf", sep=""),width=14,height=8)
}



stack_plot_row= function(Ab_all_ml_list,out_name, sample_order){
     library(ggplot2)
    library(stringr)
    library(ggpubr)
    Ab_all_ml_list_stat <- merge_nb_not_bench(Ab_all_ml_list)
    
    ## for 2strains
    Ab_all_ml_list_stat_ordered <- Ab_all_ml_list_stat
        
    ### change
    Ab_all_ml_list_stat_ordered$variable=str_replace_all(Ab_all_ml_list_stat_ordered$variable,"bench","Ground Truth")
    
    
    Ab_all_ml_list_stat_ordered$variable=factor(Ab_all_ml_list_stat_ordered$variable,
             levels = sample_order,ordered = TRUE)
    

    Ab_all_ml_list_stat_ordered$Rank=factor(Ab_all_ml_list_stat_ordered$Rank,
                                            levels=unique(Ab_all_ml_list_stat_ordered$Rank),
                                            labels=c(paste("Strain",seq(1,length(unique(Ab_all_ml_list_stat_ordered$Rank))-1),sep = ""),"Extras"))
    


    ### based on one strain for order, 4str select the strain4, and 2str select the strain1
    if(length(unique(Ab_all_ml_list_stat_ordered$Rank))==5){
      h=Ab_all_ml_list_stat_ordered[Ab_all_ml_list_stat_ordered$Rank==
                                      levels(Ab_all_ml_list_stat_ordered$Rank)[4],] %>%
        arrange(.,desc(value))
      ## change factor levels 
      Ab_all_ml_list_stat_ordered$Rank=factor(Ab_all_ml_list_stat_ordered$Rank,
                                               levels = c("Strain4","Strain3","Strain2","Strain1","Extras") )
    }else{
      h=Ab_all_ml_list_stat_ordered[Ab_all_ml_list_stat_ordered$Rank==
                                      levels(Ab_all_ml_list_stat_ordered$Rank)[1],] %>%
        arrange(.,desc(value))
    }

    
    order_first= unique(h$V2)[1]
    
    Ab_all_ml_list_stat_ordered$new_order=factor(Ab_all_ml_list_stat_ordered$V2,
                                                 levels=as.character(unique(h$V2)),
                                                 labels = seq(1,length(unique(h$V2))))
    
    Ab_all_ml_list_stat_ordered=scale_str(Ab_all_ml_list_stat_ordered,
                                          var1 = "variable",var2="new_order")
    
    # sum the data of Extras for ggpurb, 
    Ab_all_ml_list_stat_ordered1 = aggregate(scale_value~Rank+variable+sample+V2+new_order,
                                             data=Ab_all_ml_list_stat_ordered,sum)
    
 
    
    ggbarplot(Ab_all_ml_list_stat_ordered1,x="new_order",y="scale_value",fill = "Rank")+
        facet_wrap(~variable,nrow=1)+
      # scale_fill_manual(values =  c("#00AFBB", "#E7B800", "#FC4E07","green","purple"))+
        labs(x=NULL, y="Relative Abundance",fill="Strain")+
            theme( axis.text.x=element_blank(),
                   legend.text = element_text(size = 15),legend.title = element_text(size=15),axis.title=element_text(size=15))
      
                    
    ggsave(file=paste(out_name, "_abun_stackplot_by_method_ordered.pdf", sep=""),width=14,height=8)
}


######################################################################################
#Constrain format: 
#forece order
library(stringr)
sample_order<-unlist(read.table(sample_order_f,as.is=TRUE,stringsAsFactors = FALSE))
sample_order=str_replace_all(sample_order,"bench","Ground Truth")

library(ggplot2)
#all_abun_list_f = "3g_all_abun_list_f.txt"
all_abun_list <- read.table(all_abun_list_f,as.is=TRUE,stringsAsFactors = FALSE)
for(i in 1:dim(all_abun_list)[1]){
    Ab_all_tmp<-read.table(all_abun_list[i,1],header=TRUE,as.is=TRUE,stringsAsFactors = FALSE)
    if(length(which(Ab_all_tmp$variable=="shuffle"))>0){
        Ab_all_tmp <- Ab_all_tmp[-which(Ab_all_tmp$variable=="shuffle"),]
    }
    assign(all_abun_list[i,2],Ab_all_tmp)
    if(i==1){
        stack_plot_row(Ab_all_ml_list=get(all_abun_list[i,2]),out_name=all_abun_list[i,2], sample_order)
        others_plot(Ab_all_ml_list=get(all_abun_list[i,2]),out_name=all_abun_list[i,2], sample_order)
    }else {
        stack_plot_row(get(all_abun_list[i,2]),all_abun_list[i,2], sample_order)
        others_plot(Ab_all_ml_list=get(all_abun_list[i,2]),out_name=all_abun_list[i,2], sample_order)
    }
}

#plot all groups in one figure 
all_abun_list <- read.table(all_abun_list_f,as.is=TRUE,stringsAsFactors = FALSE)
for(i in 1:dim(all_abun_list)[1]){
    Ab_all_tmp<-read.table(all_abun_list[i,1],header=TRUE,as.is=TRUE,stringsAsFactors = FALSE)
    assign(all_abun_list[i,2],Ab_all_tmp)
    
    
    RA_melt_tmp <- merge_nb_not_bench(get(all_abun_list[i,2]))%>% subset(.,select=- Rank1)
    
    RA_melt_tmp <- cbind(RA_melt_tmp, rep(all_abun_list[i,2],dim(RA_melt_tmp)[2]))
    colnames(RA_melt_tmp)[6] <- "nstrs"
    
    RA_melt_tmp$Rank=factor(RA_melt_tmp$Rank,
                           levels=unique(RA_melt_tmp$Rank),
                           labels=c(paste("Strain",seq(1,length(unique(RA_melt_tmp$Rank))-1),sep = ""),"Extras"))


    if(RA_melt_tmp$nstrs=="4str"){
      h=RA_melt_tmp[RA_melt_tmp$Rank==levels(RA_melt_tmp$Rank)[1],] %>%
        arrange(.,desc(value))
    }else{
      
      h=RA_melt_tmp[RA_melt_tmp$Rank==levels(RA_melt_tmp$Rank)[1],] %>%
        arrange(.,desc(value))
    }

    #order by whole (all groups)
    #RA_melt_tmp$new_order=factor(RA_melt_tmp$V2, levels=as.character(unique(h$V2)), labels = seq(1,length(unique(h$V2))))

    #order by bench, in order to make the sample order different between str goups, need to change the sample lables
    RA_melt_tmp$new_order=factor(RA_melt_tmp$V2, levels=as.character(unlist(subset(h,variable=="bench",select=V2))), labels = seq(1,length(unique(h$V2))))
    
    
    RA_melt_tmp=scale_str(RA_melt_tmp,var1 = "variable",var2="new_order")
    
    # sum the data of others for ggpurb, 
    RA_melt_tmp = aggregate(scale_value~Rank+variable+sample+V2+new_order+nstrs,
                                             data=RA_melt_tmp,sum)

    if(i==1){
     RA_melt=RA_melt_tmp
    }else {

    RA_melt <- rbind(RA_melt,RA_melt_tmp)
    }
}



if(length(which(RA_melt$variable=="shuffle"))>0){
    RA_melt <- RA_melt[-which(RA_melt$variable=="shuffle"),]
}


RA_melt$variable=str_replace_all(RA_melt$variable,"bench","Ground Truth")

RA_melt$variable <- factor(RA_melt$variable,levels = sample_order,ordered = TRUE)

RA_melt$Rank=factor(RA_melt$Rank,levels = c(paste("Strain",seq(1,length(levels(RA_melt$Rank))-1),sep = ""),"Extras"))

ggbarplot(RA_melt,x="new_order",y="scale_value",fill="Rank")+
    facet_grid(nstrs~variable) +
    labs(x=NULL, y="Relative Abundance",fill="Strain")+
    theme( axis.text.x=element_blank(),
           legend.text = element_text(size = 15),legend.title = element_text(size=15),axis.title=element_text(size=15))

ggsave(file=paste("All_groups", "_abun_stackplot_by_method_ordered.pdf", sep=""),width=14,height=8)

