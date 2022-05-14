#functions
###########simplex plot section###############
simplex_plot_timepoint_highlight <- function(str_comp,str_comp_meta_filter,timepoint_T,color_code_T,out_header,sub_header_simplex_plot){
    T_runs<-intersect(rownames(str_comp),rownames(str_comp_meta_filter)[str_comp_meta_filter$meta_timepoint==timepoint_T])
    T_runs_minus<-intersect(rownames(str_comp),rownames(str_comp_meta_filter)[str_comp_meta_filter$meta_timepoint!=timepoint_T])
    p<-ggtern(data=str_comp_meta_filter,aes(x=strain1,y=strain2,z=strain3)) +
    geom_mask() + # 可将超出边界的点正常显示出来
    geom_point(data = str_comp_meta_filter[T_runs_minus,], size=2, color = "grey", alpha = 0.6, show.legend = FALSE) +
    geom_point(data = str_comp_meta_filter[T_runs,], aes( color = meta_timepoint), size = 2,show.legend = TRUE) +
    scale_colour_manual(values = color_code_T) +
    theme_bw( ) +    #绘图主题设置，可以通过 ?theme_bw 查看其他主题
    theme_legend_position("topright") +   #图例位置参数
    guides(size = "none") +
    theme(axis.title=element_text(size=20,face = "bold"), axis.text = element_blank(), axis.ticks = element_blank()) + #axis.title change the font type of axis names
    labs(title = paste(timepoint_T,":",length(T_runs),sep=""))    #标题
    ggsave(p, file=paste(out_header, sub_header_simplex_plot,"/simplex_plot_",timepoint_T,".pdf", sep=""))
}

simplex_plot_timepoint_subset <- function(str_comp,str_comp_meta_filter,timepoint_T,color_code_T,out_header,sub_header_simplex_plot){
    T_runs<-intersect(rownames(str_comp),rownames(str_comp_meta_filter)[str_comp_meta_filter$meta_timepoint==timepoint_T])
    T_runs_minus<-intersect(rownames(str_comp),rownames(str_comp_meta_filter)[str_comp_meta_filter$meta_timepoint!=timepoint_T])
    p<-ggtern(data=str_comp_meta_filter,aes(x=strain1,y=strain2,z=strain3)) +
    geom_mask() + # 可将超出边界的点正常显示出来
    #geom_point(data = str_comp_meta_filter[T_runs_minus,], size=2, color = "grey", alpha = 0.6, show.legend = FALSE) +
    geom_point(data = str_comp_meta_filter[T_runs,], aes( color = meta_timepoint), size = 2,show.legend = TRUE) +
    scale_colour_manual(values = color_code_T) +
    theme_bw( ) +    #绘图主题设置，可以通过 ?theme_bw 查看其他主题
    theme_legend_position("topright") +   #图例位置参数
    guides(size = "none") +
    theme(axis.title=element_text(size=20,face = "bold"), axis.text = element_blank(), axis.ticks = element_blank()) + #axis.title change the font type of axis names
    labs(title = paste(timepoint_T,":",length(T_runs),sep=""))    #标题
    ggsave(p, file=paste(out_header, sub_header_simplex_plot,"/simplex_plot_",timepoint_T,"_only.pdf", sep=""))
}

##########dist_to_M section###############
#same as pair_dist_gen, but calculate aitchison dist
pair_ait_dist_gen <- function(sample_ID_matrix,meta_data,str_comp_psuedo,group_M_n,group_N_n,extreme_psuedo){
    library(coda.base)
    MB_dist_list <- rep(100,length(unique(meta_data$meta_person)))
    names(MB_dist_list) <- rownames(sample_ID_matrix)
    for (i in 1:length(MB_dist_list)) {
        if(sample_ID_matrix[i,group_M_n]!="0" && sample_ID_matrix[i,group_N_n]!="0"){
            if(sum(str_comp_psuedo[sample_ID_matrix[i,group_M_n],])>0.5 && sum(str_comp_psuedo[sample_ID_matrix[i,group_N_n],])>0.5){
                pair_sam <- rbind(str_comp_psuedo[sample_ID_matrix[i,group_M_n],],str_comp_psuedo[sample_ID_matrix[i,group_N_n],])
                pair_sam[pair_sam==0]<-extreme_psuedo #add pseudo for log
                MB_dist_list[i] <- dist(pair_sam, method = "aitchison")
            }
        }
    }

    MB_dist_list<- MB_dist_list[MB_dist_list!=100]
    return(MB_dist_list)
    detach("package:coda.base", unload=TRUE)
}

#the following function is good for pairbase distance for ERR project with MB, M4M and M12M groups. For subset such as vaginal or CS or mother not 0,just subset the three input data(which is just a list of distance for each group.
#note: the data must be cleaned of extreme value such as 100 (used to mark empty samples)
#it is used only in dist_to_M_timepoints section
pair_bray_time <- function(MB_dist_data,M4M_dist_data,M12M_dist_data,out_csv,out_pdf_head,color_distM){
    pair_dist_pvalue <- matrix(1,nrow=4,ncol=4)
    colnames(pair_dist_pvalue) <- c("M_B","M_4M","M_12M","Null")
    rownames(pair_dist_pvalue) <- c("M_B","M_4M","M_12M","Null")
    pair_dist_list <- list(MB=MB_dist_data, M4M=M4M_dist_data, M12M=M12M_dist_data, Null=0)
    for(i in 1:dim(pair_dist_pvalue)[1]){
        for(j in 1:dim(pair_dist_pvalue)[2]){
            if(length(pair_dist_list[[i]])>1 && length(pair_dist_list[[j]])>1){
                if(var(pair_dist_list[[i]])>0 || var(pair_dist_list[[j]])>0){
                    pair_dist_pvalue[i,j] <- t.test(pair_dist_list[[i]],pair_dist_list[[j]])$p.value
                }
            }else{
                if(length(pair_dist_list[[i]])>1 && j==4 && i!=j){
                    pair_dist_pvalue[i,j] <- t.test(pair_dist_list[[i]])$p.value
                }

            }
        }
        if(i==1){
            pair_dist_data <- as.data.frame(cbind(pair_dist_list[[i]],rep(colnames(pair_dist_pvalue)[i],length(pair_dist_list[[i]]))),stringsAsFactors = FALSE)
        }else{
            pair_dist_data_tmp <- as.data.frame(cbind(pair_dist_list[[i]],rep(colnames(pair_dist_pvalue)[i],length(pair_dist_list[[i]]))),stringsAsFactors = FALSE)
            pair_dist_data <- rbind(pair_dist_data,pair_dist_data_tmp)
        }
    }
    pair_dist_pvalue[4,] <- pair_dist_pvalue[,4]
    write.csv(pair_dist_pvalue,file=out_csv)
    pair_dist_pvalue_line <- cbind(pair_dist_pvalue[1,2],pair_dist_pvalue[1,3],pair_dist_pvalue[2,3],pair_dist_pvalue[1,4],pair_dist_pvalue[2,4],pair_dist_pvalue[3,4])
    colnames(pair_dist_pvalue_line)<-c("M_BvsM_4M","M_BvsM_12M","M_4MvsM_12M","M_BvsNull","M_4MvsNull","M_12MvsNull")
    write.csv(pair_dist_pvalue_line,file=paste(out_csv,"_line_value",sep=""),row.names = FALSE)
    colnames(pair_dist_data)<-c("value","V2")
    pair_dist_data[,1]<-as.numeric(pair_dist_data[,1])
    
    #remove the Null group
    pair_dist_data_3g <- pair_dist_data[pair_dist_data[,2]!="Null",]
    pair_dist_data_3g$V2 <- factor(pair_dist_data_3g$V2, levels=c("M_B","M_4M","M_12M"), ordered=TRUE) #mark the order at the beginning

    #boxplot 改成用ggpubr的，加violin plot
    plot_meta_pair_bray_time(pair_dist_data_3g,out_pdf_head,color_distM)
    #加violine plot也是用ggpubr的
    return(pair_dist_pvalue)
}

plot_meta_pair_bray_time= function(data_in,out_name,color_distM){
    library(reshape2)
    library(ggplot2)
    data_ml <- melt(data_in)
    summary_data <- table(data_ml$V2)
    my_comparisons=list(c("M_B", "M_4M"), c("M_4M", "M_12M"), c("M_B", "M_12M"))
    p <- ggboxplot(data_ml,x="V2",y="value", color="V2", palette = color_distM, add = "jitter")+ #增加了jitter点，点shape由dose映射
    #stat_compare_means(comparisons = my_comparisons)+
    stat_compare_means(comparisons = my_comparisons,method = "t.test",label = "p.signif")
    ggsave(p, file=paste(out_name,"_boxplot.pdf",sep = "")) 
    p <- ggviolin(data_ml, x="V2", y="value", fill = "V2", palette = color_distM, add = "boxplot", add.params = list(fill="white"))+
    stat_compare_means(comparisons = my_comparisons,method = "t.test", label = "p.signif")
    ggsave(p, file=paste(out_name,"_violin.pdf",sep = "")) 

}


###################paired_dist on each strain for breast stop section##########
#get pairwise difference for group_M_n time point minus group_N_n on straini
pair_diff_gen <- function(sample_ID_matrix,meta_data,str_comp_psuedo,group_M_n,group_N_n,str_i){
    MB_diff_list <- rep(100,length(unique(meta_data$meta_person)))
    for (i in 1:length(MB_diff_list)) {
        #filter the pair with at least one sample had no result from XT
        if(sample_ID_matrix[i,group_M_n]!="0" && sample_ID_matrix[i,group_N_n]!="0"){
            MB_diff_list[i] <- str_comp_psuedo[sample_ID_matrix[i,group_M_n],str_i]-str_comp_psuedo[sample_ID_matrix[i,group_N_n],str_i]
        }
    }

    names(MB_diff_list) <- rownames(sample_ID_matrix)
    MB_diff_list<- MB_diff_list[MB_diff_list<100]
    return(MB_diff_list)
}

#do stat from the pair_dist matrix
pair_dist_stat <- function(MB_dist_list,p_G1,p_G2,out_pdf_head,color_BS){
    match_samples <- names(MB_dist_list)[MB_dist_list!=100]
    MB_dist_G1 <- MB_dist_list[intersect(p_G1,match_samples)]
    MB_dist_G2 <- MB_dist_list[intersect(p_G2,match_samples)]
    if(length(MB_dist_G1)>1 && length(MB_dist_G2)>1){
        MB_dist_stat <- t.test(MB_dist_G1,MB_dist_G2)$p.value
    }else{
        MB_dist_stat <- "Not enough value" #not enough 'y' observations的报错，在大多数菌种都会出现，因此改成直接用pvalue
    }
    MB_dist_data <- as.data.frame(cbind(c(MB_dist_G1,MB_dist_G2),c(rep("Stop_breast",length(MB_dist_G1)),rep("Cont_breast",length(MB_dist_G2)))))
    colnames(MB_dist_data)<-c("value","V2")
    MB_dist_data[,1]<-c(MB_dist_G1,MB_dist_G2)
    MB_dist_data$V2 <- factor(MB_dist_data$V2, levels=c("Stop_breast","Cont_breast"), ordered=TRUE) #mark the order at the beginning

    plot_meta_pair(MB_dist_data,out_pdf_head,color_BS)
    return(MB_dist_stat)
}

both_pair_dist_stat <- function(B4M_diff_list,p_BFF,p_others,M4M12_diff_list,p_BBF,p_BBM,out_pdf_head,color_BS){
    match_samples1 <- names(B4M_diff_list)[B4M_diff_list!=100]
    MB_dist_G1_1 <- B4M_diff_list[intersect(p_BFF,match_samples1)]
    MB_dist_G2_1 <- B4M_diff_list[intersect(p_others,match_samples1)]
    match_samples2 <- names(M4M12_diff_list)[M4M12_diff_list!=100]
    MB_dist_G1_2 <- M4M12_diff_list[intersect(p_BBF,match_samples2)]
    MB_dist_G2_2 <- M4M12_diff_list[intersect(p_BBM,match_samples2)]
    
    MB_dist_G1 <- c(MB_dist_G1_1,MB_dist_G1_2)
    MB_dist_G2 <- c(MB_dist_G2_1,MB_dist_G2_2)

    if(length(MB_dist_G1)>1 && length(MB_dist_G2)>1){
        MB_dist_stat <- t.test(MB_dist_G1,MB_dist_G2)$p.value
    }else{
        MB_dist_stat <- "Not enough value" #not enough 'y' observations的报错，在大多数菌种都会出现，因此改成直接用pvalue
    }
    MB_dist_data <- as.data.frame(cbind(c(MB_dist_G1,MB_dist_G2),c(rep("Stop_breast",length(MB_dist_G1)),rep("Cont_breast",length(MB_dist_G2)))))
    colnames(MB_dist_data)<-c("value","V2")
    MB_dist_data[,1]<-c(MB_dist_G1,MB_dist_G2)
    MB_dist_data$V2 <- factor(MB_dist_data$V2, levels=c("Stop_breast","Cont_breast"), ordered=TRUE) #mark the order at the beginning

    plot_meta_pair(MB_dist_data,out_pdf_head,color_BS)
    return(MB_dist_stat)
}

plot_meta_pair= function(data_in,out_name,color_BS){
    library(reshape2)
    library(ggplot2)
    data_ml <- melt(data_in)
    summary_data <- table(data_ml$V2)
    my_comparisons=list(c("Stop_breast", "Cont_breast"))
    p <- ggboxplot(data_ml,x="V2",y="value", color="V2", palette = color_BS, add = "jitter")+ #增加了jitter点，点shape由dose映射
    #stat_compare_means(comparisons = my_comparisons)+
    stat_compare_means(comparisons = my_comparisons,method = "t.test", label = "p.signif")
    ggsave(p, file=paste(out_name,"_boxplot.pdf",sep = "")) 
    p <- ggboxplot(data_ml,x="V2",y="value", color="V2", palette = color_BS, add = "jitter")+ #增加了jitter点，点shape由dose映射
    coord_cartesian(ylim = c(-1, 1)) + 
    stat_compare_means(comparisons = my_comparisons,method = "t.test", label = "p.signif") +
    labs(title = paste(table(data_ml$V2),collapse=";"))
    ggsave(p, file=paste(out_name,"_yfull_boxplot.pdf",sep = "")) 
    p <- ggviolin(data_ml, x="V2", y="value", fill = "V2", palette = color_BS, add = "boxplot", add.params = list(fill="white"))+
    stat_compare_means(comparisons = my_comparisons,method = "t.test", label = "p.signif")
    ggsave(p, file=paste(out_name,"_violin.pdf",sep = "")) 

}
 #this if for N*M groups (which means the rownames and colnames are different, so will go through all the grids)
check_significant_general <- function(P_summary){
    group_nr <- dim(P_summary)[1]
    group_nc <- dim(P_summary)[2]  
    for(i in 1:group_nr){
        for(j in 1:group_nc){
            if(P_summary[i,j]<0.05){
                if(exists("sig_list")){
                    sig_list <- rbind(sig_list,c(rownames(P_summary)[i],colnames(P_summary)[j],P_summary[i,j]))
                }else{
                    sig_list <- c(rownames(P_summary)[i],colnames(P_summary)[j],P_summary[i,j])
                }
            }
        }    
    }
    if(exists("sig_list")){
        if(sum(dim(sig_list))>0){
            colnames(sig_list)<-c("G1","G2","p-value")
        }
    }else{
        sig_list <- ""
    }
    return(sig_list)
    rm(sig_list)
}



################### M-B dist boxplot and stat  ##########
#do three different distance all together. p_vaginal and p_CS can be replaced by other people group
three_dis_gen <- function(sample_ID_matrix,meta_data,str_comp_psuedo,sub_header_VC_dist, group1_n,group2_n,groupname,comparename,p_vaginal,p_CS, extreme_psuedo){
    MB_dist_list<- pair_dist_gen(sample_ID_matrix,meta_data,str_comp_psuedo,group1_n,group2_n)
    out_pdf_head=paste(out_header, sub_header_VC_dist,"/",groupname,"_braydist_",comparename, sep="")
    MB_dist_stat <-pair_dist_stat_VC(MB_dist_list,p_vaginal,p_CS,out_pdf_head,color_VC)
    MB_ait_dist_list<- pair_ait_dist_gen(sample_ID_matrix,meta_data,str_comp_psuedo,group1_n,group2_n,extreme_psuedo)
    out_pdf_head=paste(out_header, sub_header_VC_dist,"/",groupname,"_aitdist_",comparename, sep="")
    MB_dist_stat_ait <-pair_dist_stat_VC(MB_ait_dist_list,p_vaginal,p_CS,out_pdf_head,color_VC)
    MB_euc_dist_list<- pair_euc_dist_gen(sample_ID_matrix,meta_data,str_comp_psuedo,group1_n,group2_n)
    out_pdf_head=paste(out_header, sub_header_VC_dist,"/",groupname,"_eucdist_",comparename, sep="")
    MB_dist_stat_euc <-pair_dist_stat_VC(MB_euc_dist_list,p_vaginal,p_CS,out_pdf_head,color_VC)
    list_pvalue <- c(MB_dist_stat, MB_dist_stat_ait, MB_dist_stat_euc)
    names(list_pvalue) <- c("bray","aitchison","euclidean")
    return(list_pvalue)
}

#generate pairwise dist for 2 groups using three inputs: sample_matrix(sample_ID_matrix),meta_data of person, string composition(str_comp_psuedo). group_M_n and group_N_n are the two column number of selected group in the sample_matrix
pair_dist_gen <- function(sample_ID_matrix,meta_data,str_comp_psuedo,group_M_n,group_N_n){
    MB_dist_list <- rep(100,length(unique(meta_data$meta_person)))
    for (i in 1:length(MB_dist_list)) {
        if(sample_ID_matrix[i,group_M_n]!="0" && sample_ID_matrix[i,group_N_n]!="0"){
            if(sum(str_comp_psuedo[sample_ID_matrix[i,group_M_n],])>0.5 && sum(str_comp_psuedo[sample_ID_matrix[i,group_N_n],])>0.5){ #there could not have empty sample in a pair
                pair_sam <- rbind(str_comp_psuedo[sample_ID_matrix[i,group_M_n],],str_comp_psuedo[sample_ID_matrix[i,group_N_n],])
                MB_dist_list[i] <- vegdist(pair_sam, method = "bray")
            }
        }
    }

    names(MB_dist_list) <- rownames(sample_ID_matrix)
    MB_dist_list<- MB_dist_list[MB_dist_list<100]
    return(MB_dist_list)
}

#same as pair_dist_gen, but calculate euclidean dist
pair_euc_dist_gen <- function(sample_ID_matrix,meta_data,str_comp_psuedo,group_M_n,group_N_n){
    MB_dist_list <- rep(100,length(unique(meta_data$meta_person)))
    for (i in 1:length(MB_dist_list)) {
        if(sample_ID_matrix[i,group_M_n]!="0" && sample_ID_matrix[i,group_N_n]!="0"){
            if(sum(str_comp_psuedo[sample_ID_matrix[i,group_M_n],])>0.5 && sum(str_comp_psuedo[sample_ID_matrix[i,group_N_n],])>0.5){
                pair_sam <- rbind(str_comp_psuedo[sample_ID_matrix[i,group_M_n],],str_comp_psuedo[sample_ID_matrix[i,group_N_n],])
                MB_dist_list[i] <- vegdist(pair_sam, method = "euclidean")
            }
        }
    }

    names(MB_dist_list) <- rownames(sample_ID_matrix)
    MB_dist_list<- MB_dist_list[MB_dist_list<100]
    return(MB_dist_list)
}

pair_dist_stat_VC <- function(MB_dist_list,p_G1,p_G2,out_pdf_head,color_BS){
    match_samples <- names(MB_dist_list)[MB_dist_list!=100]
    MB_dist_G1 <- MB_dist_list[intersect(p_G1,match_samples)]
    MB_dist_G2 <- MB_dist_list[intersect(p_G2,match_samples)]
    if(length(MB_dist_G1)>1 && length(MB_dist_G2)>1){
        MB_dist_stat <- t.test(MB_dist_G1,MB_dist_G2)$p.value
    }else{
        MB_dist_stat <- "Not enough value" #not enough 'y' observations的报错，在大多数菌种都会出现，因此改成直接用pvalue
    }
    MB_dist_data <- as.data.frame(cbind(c(MB_dist_G1,MB_dist_G2),c(rep("Vaginal",length(MB_dist_G1)),rep("C-section",length(MB_dist_G2)))))
    colnames(MB_dist_data)<-c("value","V2")
    MB_dist_data[,1]<-c(MB_dist_G1,MB_dist_G2)
    MB_dist_data$V2 <- factor(MB_dist_data$V2, levels=c("Vaginal","C-section"), ordered=TRUE) #mark the order at the beginning

    plot_meta_pair_VC(MB_dist_data,out_pdf_head,color_BS)
    return(MB_dist_stat)
}


plot_meta_pair_VC= function(data_in,out_name,color_BS){
    library(reshape2)
    library(ggplot2)
    data_ml <- melt(data_in)
    summary_data <- table(data_ml$V2)
    p <- ggboxplot(data_ml,x="V2",y="value", color="V2", palette = color_BS, add = "jitter")+ #增加了jitter点，点shape由dose映射
    #stat_compare_means(comparisons = my_comparisons)+
    stat_compare_means(method = "t.test", label = "p.signif")
    ggsave(p, file=paste(out_name,"_boxplot.pdf",sep = "")) 
    p <- ggviolin(data_ml, x="V2", y="value", fill = "V2", palette = color_BS, add = "boxplot", add.params = list(fill="white"))+
    stat_compare_means(method = "t.test", label = "p.signif")
    ggsave(p, file=paste(out_name,"_violin.pdf",sep = "")) 

}







