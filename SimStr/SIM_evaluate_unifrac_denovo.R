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


#usage example: Rscript /home/yxtan/StrainSIM/SIM_evaluate_unifrac_denovo.R bench_tb=/mnt/data2/LD_lab/yxtan/StrainAnal_tools/MockData/Xtrain_bench/Ecoli98/bench_data/bench_Ecoli98_WGS_all_str_prof.csv method_anno_tb=/mnt/data2/LD_lab/yxtan/StrainAnal_tools/MockData/Xtrain_bench/Ecoli98/bench_data/WGS/XT_all/xstrain_out/Escherichia-coli-202006_xstrainr_out/Ecoli98_WGS_all_str_anno_prof.csv denovo_ref=/mnt/data2/LD_lab/yxtan/StrainAnal_tools/MockData/Xtrain_bench/Ecoli98/Eval/Xtrain_str_list.txt tree_nwk=/mnt/data2/LD_lab/yxtan/StrainAnal_tools/MockData/Xtrain_bench/Ecoli98/Eval/P_2020_07_09_204331864390/parsnp.tree name_out=Ecoli98_WGS_all_Xtrain
#The converted output from a de novo methed will be compared with the bench with the integrated tree.
#All methods that generate de novo strains will use this solution. (From which the reference strains were not extractable)
#Inputs: 1. converted bench table;2. converted/annotated abun matrix; 3. integrated tree. 4. reference list.
#the output will be the distance and stat values of each sample and plots
#library(vegan) library(ggplot2) library(reshape2) library(ape) , library (rbiom) are requried

#update log:
#1. compared to oldv1, the bench_neighbor function was corrected, fingding the neighbor of predicted from the ref_list and SIM, rather than finding the neighbor of SIM from the predicted (this will overestimate the coreectness). 
#2. Add the sum_dup_row function in, to merge strains with same SIM neighbor (only consider the case when the neighbor is from SIM, all others will just keep the original)
#delete the old useless commented codes.
#update log 20200928
#add JSD by phyloseq
#add MCC
#add relative abundance(without ranked)

#update log 20201017
#add the FP group for ranked data by add rank N+1. This will make the across method plot have the rank N+1 group as the representative of the sum of FPs

#update log 20201207
#modify the neighbor function, use given range of SIM first, if not within the range, find the strs' closest neighbor instead (not necessary to be SIM str now). The default range now is 0.01 as an inputable parameter.

#update log 20210505
#updatede the bench_neighbor function to avoid multi output for each ID, because multi output made the neighbor relationship incorrect.

#update log 20210524
#add the JSD function for sorted relative abundance
#add the str_match parameter to control the calculation of matched relative abundance of strains
#pack functions(GCF_GCA_convert,annotation_check)
#add low_abun_filter for filtering abundance values smaller than it

#update log 20210603
#replace the first column test by checking whether the first unit is factor instead of using the name of first column(which will cause error for some methods' input)

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

#check whether the paras are exist
if (!exists("bench_tb")) {
    stop(paste("\nRscript SIM_evaluate_unifrac_denovo.R bench_tb=",bench_tb," \nWarning: Usage: simulation bench para is not given, please check the input parah. \n\n",sep=""))
}

if (!exists("method_anno_tb")) {
    stop(paste("\nRscript SIM_evaluate_unifrac_denovo.R method_anno_tb=",method_anno_tb," \nWarning: Usage: Xtrain annotated output file para is not given, please check the input para. \n\n",sep=""))
}

if (!exists("denovo_ref")) {
    stop(paste("\nRscript SIM_evaluate_unifrac_denovo.R denovo_ref=",denovo_ref," \nWarning: Usage: The ref database file in XXXX format of the de novo method para is not given, please check the input para. \n\n",sep=""))
}

if (!exists("tree_nwk")) {
    stop(paste("\nRscript SIM_evaluate_unifrac_denovo.R tree_nwk=",tree_nwk," \nWarning: Usage: The tree file of all strains para is not given, please check the input para. \n\n",sep=""))
}

if (!exists("dist_range")) {
    dist_range=0.05
    print(paste("\nRscript SIM_evaluate_unifrac_denovo.R dist_range=",dist_range," \nWarning: The dist_range parameter was not input, use its default value 0.05. \n\n",sep=""))
}

if (!exists("low_abun_filter")) {
    low_abun_filter=0.01
    print(paste("\nRscript SIM_evaluate_unifrac_denovo.R low_abun_filter=",low_abun_filter," \nWarning: The low_abun_filter parameter was not input, use its default value 0.01. \n\n",sep=""))
}

if (!exists("str_match")) {
    str_match=TRUE
    print(paste("\nRscript SIM_evaluate_unifrac_denovo.R str_match=",str_match," \nWarning: The str_match parameter was not input, use its default value is TRUE to calculate matched relative abundance for strains. \n\n",sep=""))
}


if (!exists("name_out")) {
    name_out="out"
    print(paste("\nRscript SIM_evaluate_unifrac_denovo.R name_out=",name_out," \nWarning: The name_out parameter was not input, use its default value out. \n\n",sep=""))
}


#check whether the file in is exist
if (!file.exists(bench_tb)) {
    stop(paste("\nRscript SIM_evaluate_unifrac_denovo.R bench_tb=",bench_tb," \nWarning: Usage: simulation bench file is not exist, please check the path. \n\n",sep=""))
}

if (!file.exists(method_anno_tb)) {
    stop(paste("\nRscript SIM_evaluate_unifrac_denovo.R method_anno_tb=",method_anno_tb," \nWarning: Usage: Xtrain annotated output file is not exist, please check the path. \n\n",sep=""))
}

if (!file.exists(denovo_ref)) {
    stop(paste("\nRscript SIM_evaluate_unifrac_denovo.R denovo_ref=",denovo_ref," \nWarning: Usage: The ref database file in XXXX format of the de novo method is not exist, please check the path. \n\n",sep=""))
}

if (!file.exists(tree_nwk)) {
    stop(paste("\nRscript SIM_evaluate_unifrac_denovo.R tree_nwk=",tree_nwk," \nWarning: Usage: The tree file of all strains is not exist, please check the path. \n\n",sep=""))
}

#functions
#This function is used to find the pairwised relation between predicted and SIM,each SIM match to 1 predicted ID.
bench_neighbor = function(denovo_ref_ID,SIM_ref_ID,dist_matrix,ref_list_ID,dist_range){
    ID_inter <- intersect(denovo_ref_ID,SIM_ref_ID)
    ID_diff <-  setdiff(denovo_ref_ID,SIM_ref_ID)
    ID_union_SIM <- union(denovo_ref_ID,SIM_ref_ID)
    
    #check whehter the key strains(denovoID\SIMID) are in the tree (in case the strains were not covered in the download list or droped by parsnp)
    tree_ID <- rownames(dist_matrix)
    denovo_diff_ID <- setdiff(denovo_ref_ID,tree_ID)
    SIM_diff_ID <- setdiff(SIM_ref_ID,tree_ID)
    ref_diff_ID <- setdiff(ref_list_ID,tree_ID)
    ref_inter_ID <- intersect(ref_list_ID,tree_ID)
    if (length(SIM_diff_ID)>0){
        stop(paste("\nIn SIM_evaluate_unifrac_denovo.R SIM_ID=(",paste(SIM_diff_ID,collapse=";"),") were missed in the tree, exit. \n Please check the version of the tree file and its dependent files. \n\n",sep=""))
    }
    if (length(denovo_diff_ID)>0){
        stop(paste("\nIn SIM_evaluate_unifrac_denovo.R predicted_ID=(",paste(denovo_diff_ID,collapse=";"),") were missed in the tree, exit. \n Please check the version of the tree file and its dependent files. \n\n",sep=""))
    }
    if (length(ref_diff_ID)>0){
        print(paste("\nIn SIM_evaluate_unifrac_denovo.R ref_ID=(",paste(ref_diff_ID,collapse=";"),") were missed in the tree, but won't affect the evaluation too much. \n If you want to make the evaluation perfect, please check the version of the tree file and its dependent files to be updated. \n\n",sep=""))
    }

    ID_union <- union(ID_union_SIM,ref_inter_ID)
    dist_m<-dist_matrix[ID_diff,ID_union]
    if (length(ID_inter)==0){
        if(length(ID_diff)==1){
            #check the distance first, if it is closed to one of the SIM—ref (within range), use it directly. Or, use the neighbor
            if(min(dist_m[SIM_ref_ID]) < dist_range){
                neighbor_l=names(dist_m[which(dist_m==min(dist_m[SIM_ref_ID]))])
                if(length(neighbor_l)>1){
                    if(length(intersect(neighbor_l,SIM_ref_ID))>0){
                        neighbor_l <- intersect(neighbor_l,SIM_ref_ID)[1]
                    }else{neighbor_l <-neighbor_l[1]}
                }
            } else{
                neighbor_l=names(sort(dist_m)[2])                
            }
            names(neighbor_l) <- ID_diff
        }else{
            if(min(dist_m[ID_diff[1],SIM_ref_ID]) < dist_range){
                neighbor_l=names(which(dist_m[ID_diff[1],]==min(dist_m[ID_diff[1],SIM_ref_ID])))
                if(length(neighbor_l)>1){
                    if(length(intersect(neighbor_l,SIM_ref_ID))>0){
                        neighbor_l <- intersect(neighbor_l,SIM_ref_ID)[1]
                    }else{neighbor_l <-neighbor_l[1]}
                }
            } else{
                neighbor_l=names(sort(dist_m[ID_diff[1],])[2])                
            }
            
            for(i in 2:length(ID_diff)){
                if(min(dist_m[ID_diff[i],SIM_ref_ID]) < dist_range){
                    neighbor_i=names(which(dist_m[ID_diff[i],]==min(dist_m[ID_diff[i],SIM_ref_ID])))
                    if(length(neighbor_i)>1){
                        if(length(intersect(neighbor_i,SIM_ref_ID))>0){
                            neighbor_i <- intersect(neighbor_i,SIM_ref_ID)[1]
                        }else{neighbor_i <-neighbor_i[1]}
                    }
                } else{
                    neighbor_i=names(sort(dist_m[ID_diff[i],])[2])                
                }

                neighbor_l <- c(neighbor_l,neighbor_i)
            }
            names(neighbor_l) <- ID_diff
        }
    } else {
        if(length(ID_diff)==1){
            neighbor_l <- ID_inter
            if(min(dist_m[SIM_ref_ID]) < dist_range){
                neighbor_i=names(dist_m[which(dist_m==min(dist_m[SIM_ref_ID]))])
                if(length(neighbor_i)>1){
                    if(length(intersect(neighbor_i,SIM_ref_ID))>0){
                        neighbor_i <- intersect(neighbor_i,SIM_ref_ID)[1]
                    }else{neighbor_i <-neighbor_i[1]}
                }
            } else{
                neighbor_i=names(sort(dist_m)[2])                
            }
            neighbor_l <- c(neighbor_l,neighbor_i)
            names(neighbor_l) <- c(ID_inter,ID_diff)
        }else{
            neighbor_l <- ID_inter
            for(i in ID_diff){
                if(min(dist_m[i,SIM_ref_ID]) < dist_range){
                    neighbor_i=names(which(dist_m[i,]==min(dist_m[i,SIM_ref_ID])))
                    if(length(neighbor_i)>1){
                        if(length(intersect(neighbor_i,SIM_ref_ID))>0){
                            neighbor_i <- intersect(neighbor_i,SIM_ref_ID)[1]
                        }else{neighbor_i <-neighbor_i[1]}
                    }
                } else{
                    neighbor_i=names(sort(dist_m[i,])[2])                
                }
                neighbor_l <- c(neighbor_l,neighbor_i)
            }
            names(neighbor_l) <- c(ID_inter,ID_diff)
        }
    }
    return(neighbor_l)
}
#filter empty rows, in order to filter out clustered ref strains for methods like EST
predict_filter = function(md_in){
    if(dim(md_in)[2]==1){
        md_tmp <- matrix(0, nrow=sum(md_in>0),ncol=1)
        rownames(md_tmp)<-names(md_in[md_in>0,1])
        md_tmp[,1]<-md_in[md_in>0,1]
        colnames(md_tmp)<-colnames(md_in)
    } else{
        md_tmp <- md_in[rowSums(md_in)>0,]
    }
    return(md_tmp)
}

#这个function，根据denovo_ref_ID确定是否有重复，其中names是对应md_in的行名，内容是可能出现重复的值。因此，使用时，对重复值，要找到其对应的行名，然后再进行合并。如果没有重复，或者全都是重复，就另外简单处理
sum_dup_row = function(md_in,denovo_ref_ID){
    dup_rows <- table(denovo_ref_ID)
    if (max(dup_rows)==1){
        md_tmp =md_in
        rownames(md_tmp) <- denovo_ref_ID
    }else{
        if (length(dup_rows)>1){
            md_tmp <- as.data.frame(matrix(1,nrow=length(dup_rows),ncol=dim(md_in)[2]))
            colnames(md_tmp) <- colnames(md_in)
            rownames(md_tmp) <- names(dup_rows)
            for (ni in 1:length(dup_rows)){
                if (dup_rows[ni]==1){
                    md_tmp[ni,]<-md_in[names(denovo_ref_ID[denovo_ref_ID==names(dup_rows[ni])]),]
                }else{
                    md_tmp[ni,]<-colSums(md_in[names(denovo_ref_ID[denovo_ref_ID==names(dup_rows[ni])]),])
                }
            }
        }else{
            md_tmp <- t(colSums(md_in))
            rownames(md_tmp) <- names(dup_rows)
        }
    }    
    return(md_tmp)
} 

#convert GCF into GCA format
GCF_GCA_convert = function(ID_in){
    if (sum(substr(ID_in,1,3)=="GCF")>0){
        ID_in <- chartr('GCF','GCA',ID_in)
    }
    for (i in 1:length(ID_in)){
        ID_in[i]<-strsplit(ID_in[i],split='.',fixed=TRUE)[[1]][1]
    }
    return(ID_in)
}

annotation_check = function(ID_in, file_name){
    if (sum(substr(ID_in,1,3)!="GCF")>0){
        if (sum(substr(ID_in,1,3)!="GCA")>0){
            stop(paste("\nThe ",file_name, " is not annotated in GCF or GCA format \n Please convert it into the right annotation, Exit! \n\n",sep=""))
        }
    }
}


#Step1: get tables
#SIM_bench
bench_in <- read.csv(bench_tb)
if (colnames(bench_in)[1]=="X"){
    rownames(bench_in) <- bench_in[,1]
    bench_in <- bench_in[,-1]
}
SIM_ref_ID <- rownames(bench_in)


#anno_matrix,all files other than CSV, will be consider tab deminated.
if (substr(method_anno_tb,nchar(method_anno_tb)-2,nchar(method_anno_tb))=="csv"){
        md_in <- read.csv(method_anno_tb)
        if (class(md_in[1,1])=="factor"){
            rownames(md_in) <- md_in[,1]
            if (dim(md_in)[2]==2){
                md_tmp <- matrix(0, nrow=dim(md_in)[1],ncol=1)
                md_tmp[,1] <- md_in[,2]
                rownames(md_tmp)<-md_in[,1]
                colnames(md_tmp)<-colnames(md_in)[2]
                md_in <- md_tmp
                
            } else {md_in <- md_in[,-1]}
        }
}else{
        md_in <- read.table(method_anno_tb, header=TRUE)
        rownames(md_in) <- md_in[,1]
        if (dim(md_in)[2]==2){
                md_tmp <- matrix(0, nrow=dim(md_in)[1],ncol=1)
                md_tmp[,1] <- md_in[,2]
                rownames(md_tmp)<-md_in[,1]
                colnames(md_tmp)<-colnames(md_in)[2]
                md_in <- md_tmp
                
        } else {md_in <- md_in[,-1]}
}
md_in <- predict_filter(md_in)

md_in_sample <- colnames(md_in)
for (i in 1:length(md_in_sample)){
    md_in_sample[i]<-strsplit(md_in_sample[i],split='.sorted.bam',fixed=TRUE)[[1]][1]
}
colnames(md_in)<-md_in_sample

##the following part is required for str_match=TRUE only
if (str_match==TRUE){
    ##work on the bench
    #check whether is in GCA or GCF format.
    annotation_check(SIM_ref_ID, "simulation table" )

    #because GCA is generally cover more strains than GCF
    #convert GCF into GCA format
    SIM_ref_ID <- GCF_GCA_convert(SIM_ref_ID)
    rownames(bench_in) <- SIM_ref_ID


    #work on the output from strain prediction method
    denovo_ref_ID<-rownames(md_in)
    annotation_check(denovo_ref_ID, "prediction output" )
    #convert GCF into GCA format
    denovo_ref_ID <- GCF_GCA_convert(denovo_ref_ID)
    rownames(md_in) <- denovo_ref_ID

    #ref_list
    ref_list <- read.table(denovo_ref, header=FALSE)
    ref_list_ID <- ref_list[,1]
    #check whether is in GCA or GCF format.
    annotation_check(ref_list_ID, "ref list table" )
    #convert GCF into GCA format
    ref_list_ID <- GCF_GCA_convert(ref_list_ID)


    #tree file
    library(ape)
    trefile <- ape::read.tree(tree_nwk)
    #drop the ref tip to avoid duplication
    ref_tip_ID<- trefile$tip.label[grepl(".ref",trefile$tip.label)]
    trefile = drop.tip(trefile,ref_tip_ID)
    #unify the strain names format
    dist_matrix_ID <- trefile$tip.label
    annotation_check(dist_matrix_ID, "tree file" )
    #convert GCF into GCA format
    dist_matrix_ID <- GCF_GCA_convert(dist_matrix_ID)
    trefile$tip.label <- dist_matrix_ID
    #get distance matrix
    dist_matrix <- cophenetic.phylo(trefile)
    dist_matrix_ID <- rownames(dist_matrix)

    #The orighinal assumption is that all the de novo outputs were converted to the neighbors on the bench already. Howver, this is not always true, especially when the bench str is not in the referecen database of the method.
    #As a result, need to get the neighbors of bench again and replace the original annotation.
    #get the neighbors for bench strs
    #######denovo_ref
    neighbor_str <- bench_neighbor(denovo_ref_ID,SIM_ref_ID,dist_matrix,ref_list_ID,dist_range)
    write.table(neighbor_str, file=(paste(name_out,"_all_neighbor.csv",sep="")),sep=",", quote=F, row.names= TRUE, col.names= FALSE, fileEncoding="UTF-8" )
    #rename the ID if their neighbor is SIM str; If not keep the original one.
    names(denovo_ref_ID) <- denovo_ref_ID
    denovo_ref_ID[names(neighbor_str[neighbor_str %in% SIM_ref_ID])]<-neighbor_str[neighbor_str %in% SIM_ref_ID]

    #merge strains with the same SIM neighbor
    m_out_anno<- sum_dup_row(md_in,denovo_ref_ID)
    m_out_anno[is.na(m_out_anno)]<-0
    write.table(m_out_anno, file=(paste(name_out,"_bench_neighbor_prof.csv",sep="")),sep=",", quote=F, row.names= TRUE, col.names= TRUE, fileEncoding="UTF-8" )

    #filter the low relative abudance strains, default is 0.01.
    m_out_anno[m_out_anno<low_abun_filter]=0

    #Step2: calculate distances
    otu_in<-union(rownames(m_out_anno),rownames(bench_in))
    dis_matrix_out <- matrix(1,nrow=5,ncol=length(colnames(bench_in))) #in case some samples were missed by xtrain, then the distances must all be 1.
    rownames(dis_matrix_out) <- c("UnW_D","Weight_D","Bray_D","JC_D","JSD")
    colnames(dis_matrix_out) <- colnames(bench_in)

    #generate tax matrix
    tax_matrix <- matrix(otu_in,nrow=length(otu_in),ncol=7)
    colnames(tax_matrix) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    rownames(tax_matrix) <- otu_in

    #calculate TP\FP\F-measure by anno result
    Postive_matrix <- matrix(0,nrow=20,ncol=length(colnames(bench_in))) #in case some samples were missed by xtrain, then the distances must all be 1.
    rownames(Postive_matrix) <- c("TP","FP","TN","FN","Sensitivity","Specificity","Precision","Negative-Predicted","F-measure","TP-adj","FP-adj","TN-adj","FN-adj","TPR-adj","TNR--adj","PPV-adj","NPV-adj","FM-adj","MCC","MCC-adj")
    colnames(Postive_matrix) <- colnames(bench_in)

    rown_bench<-length(rownames(bench_in))
    rank_abun_matrix <- matrix(0,nrow=(rown_bench+1)*2,ncol=length(colnames(bench_in)))
    colnames(rank_abun_matrix) <- colnames(bench_in)
    rownames(rank_abun_matrix) <- c(paste("bench_rank",1:(rown_bench+1),sep=""),paste("predict_rank",(1:(rown_bench+1)),sep=""))
    abun_matrix <- matrix(0,nrow=length(otu_in)*2,ncol=length(colnames(bench_in)))
    colnames(abun_matrix) <- colnames(bench_in)
    rownames(abun_matrix) <- c(paste("bench_",otu_in,sep=""),paste("predict_",otu_in,sep=""))


    #generate otu matrixs
    inter_samples <- intersect(colnames(m_out_anno),colnames(bench_in))
    library (rbiom) 
    library(vegan)
    library(phyloseq)
    for (sample_id in inter_samples){
    	otu_matrix <- matrix(0, nrow = length(otu_in), ncol = 2)
    	colnames(otu_matrix) <- c("Pred_out","Bench")
    	rownames(otu_matrix) <- otu_in
    	otu_matrix[rownames(m_out_anno),1]<-m_out_anno[rownames(m_out_anno),sample_id]
    	otu_matrix[rownames(bench_in),2]<-bench_in[rownames(bench_in),sample_id]
    	otu_matrix_binary <- (otu_matrix>0)*1
    	otu_matrix_binary_sum <- apply(otu_matrix_binary,1,sum)
    	otu_matrix_binary_minus <- otu_matrix_binary[,1]-otu_matrix_binary[,2]
    	otu_matrix_min <- apply(otu_matrix,1,min)
    	otu_matrix_minus <- otu_matrix[,1]-otu_matrix[,2]
    	Postive_matrix[1,sample_id] <- sum(otu_matrix_binary_sum==2)
    	Postive_matrix[2,sample_id] <- sum(otu_matrix_binary_minus ==1)
    	Postive_matrix[3,sample_id] <- length(dist_matrix_ID)-sum(otu_matrix_binary_sum==2)-sum(otu_matrix_binary_minus ==-1)-sum(otu_matrix_binary_minus ==1)
    	Postive_matrix[4,sample_id] <- sum(otu_matrix_binary_minus ==-1)
    	Postive_matrix[5,sample_id] <- Postive_matrix[1,sample_id]/(Postive_matrix[1,sample_id]+Postive_matrix[4,sample_id])
    	Postive_matrix[6,sample_id] <- Postive_matrix[3,sample_id]/(Postive_matrix[3,sample_id]+Postive_matrix[2,sample_id])
    	Postive_matrix[7,sample_id] <- Postive_matrix[1,sample_id]/(Postive_matrix[1,sample_id]+Postive_matrix[2,sample_id])
    	Postive_matrix[8,sample_id] <- Postive_matrix[3,sample_id]/(Postive_matrix[3,sample_id]+Postive_matrix[4,sample_id])
    	Postive_matrix[9,sample_id] <- 2*Postive_matrix[1,sample_id]/(2*Postive_matrix[1,sample_id]+Postive_matrix[2,sample_id]+Postive_matrix[4,sample_id])
    	Postive_matrix[10,sample_id] <- sum(otu_matrix_min)
    	Postive_matrix[11,sample_id] <- sum(otu_matrix_minus[which( otu_matrix_minus > 0)])
    	Postive_matrix[12,sample_id] <- 0 # There is no way to define true negative here, since it is value is 0。
    	Postive_matrix[13,sample_id] <- -sum(otu_matrix_minus[which(otu_matrix_minus < 0)])
    	Postive_matrix[14,sample_id] <- Postive_matrix[10,sample_id]/(Postive_matrix[10,sample_id]+Postive_matrix[13,sample_id])
    	Postive_matrix[15,sample_id] <- Postive_matrix[12,sample_id]/(Postive_matrix[12,sample_id]+Postive_matrix[11,sample_id])
    	Postive_matrix[16,sample_id] <- Postive_matrix[10,sample_id]/(Postive_matrix[10,sample_id]+Postive_matrix[11,sample_id])
    	Postive_matrix[17,sample_id] <- Postive_matrix[12,sample_id]/(Postive_matrix[12,sample_id]+Postive_matrix[13,sample_id])
    	Postive_matrix[18,sample_id] <- 2*Postive_matrix[10,sample_id]/(2*Postive_matrix[10,sample_id]+Postive_matrix[11,sample_id]+Postive_matrix[13,sample_id])
    	Postive_matrix[19,sample_id] <- (Postive_matrix[1,sample_id]*Postive_matrix[3,sample_id]-Postive_matrix[2,sample_id]*Postive_matrix[4,sample_id])/sqrt((Postive_matrix[1,sample_id]+Postive_matrix[2,sample_id])*(Postive_matrix[1,sample_id]+Postive_matrix[4,sample_id])*(Postive_matrix[2,sample_id]+Postive_matrix[3,sample_id])*(Postive_matrix[4,sample_id]+Postive_matrix[3,sample_id]))
        Postive_matrix[20,sample_id] <- (Postive_matrix[10,sample_id]*Postive_matrix[12,sample_id]-Postive_matrix[11,sample_id]*Postive_matrix[13,sample_id])/sqrt((Postive_matrix[10,sample_id]+Postive_matrix[11,sample_id])*(Postive_matrix[10,sample_id]+Postive_matrix[13,sample_id])*(Postive_matrix[11,sample_id]+Postive_matrix[12,sample_id])*(Postive_matrix[13,sample_id]+Postive_matrix[12,sample_id]))
        
    	
    	#calculate regular unifrac from annotaion data
        dis_matrix_out[1,sample_id] <- rbiom::unifrac(otu_matrix, weighted=FALSE, tree=trefile)
        dis_matrix_out[2,sample_id] <- rbiom::unifrac(otu_matrix, weighted=TRUE, tree=trefile)
        dis_matrix_out[3,sample_id] <- vegdist(t(otu_matrix), method = "bray")
        dis_matrix_out[4,sample_id] <- vegdist(t(otu_matrix), method = "jaccard",binary=TRUE) #if binary=FALSE，this distance is equal to 2B/(1+B) [B is bray distance], which is meaningless.
        phylo_test <- phyloseq(otu_table(t(otu_matrix), taxa_are_rows = F))
        dis_matrix_out[5,sample_id] <- distance(phylo_test, "jsd")
        abun_matrix[,sample_id] <- c(otu_matrix[,"Bench"], otu_matrix[,"Pred_out"])
        #Step3:calculate the rank abundance
        Pred_rank <- sort(otu_matrix[,"Pred_out"],decreasing = TRUE)
        bench_rank <- sort(otu_matrix[,"Bench"],decreasing = TRUE)
        rank_abun_matrix[,sample_id] <- c(bench_rank[1:rown_bench], 0,Pred_rank[1:rown_bench],1-sum(Pred_rank[1:rown_bench]))
    }
    rank_abun_diff_matrix <- rank_abun_matrix[1:rown_bench,]-rank_abun_matrix[(rown_bench+2):((rown_bench+1)*2-1),]
    rownames(rank_abun_diff_matrix) <- paste("Rank",1:rown_bench,sep="")

    abun_diff_matrix <- abun_matrix[1:length(otu_in),]-abun_matrix[(length(otu_in)+1):(length(otu_in)*2),]
    rownames(abun_diff_matrix) <- otu_in


    #get sample with no result
    none_emp <- intersect(colnames(m_out_anno),colnames(bench_in))
    em_rows <- setdiff(colnames(bench_in),colnames(m_out_anno))

    write.table(dis_matrix_out, file=(paste(name_out,"_dis_all.csv",sep="")),sep=",", quote=F, row.names= TRUE, col.names= TRUE, fileEncoding="UTF-8" )
    write.table(dis_matrix_out[,none_emp], file=(paste(name_out,"_dis_detected.csv",sep="")),sep=",", quote=F, row.names= TRUE, col.names= TRUE, fileEncoding="UTF-8" )
    write.table(em_rows, file=(paste(name_out,"_sample_missed.csv",sep="")),sep=",", quote=F, row.names= FALSE, col.names= FALSE, fileEncoding="UTF-8" )
    write.table(Postive_matrix, file=(paste(name_out,"_stat_all.csv",sep="")),sep=",", quote=F, row.names= TRUE, col.names= TRUE, fileEncoding="UTF-8" )
    write.table(Postive_matrix[,none_emp], file=(paste(name_out,"_stat_detected.csv",sep="")),sep=",", quote=F, row.names= TRUE, col.names= TRUE, fileEncoding="UTF-8" )

    write.table(rank_abun_matrix, file=(paste(name_out,"_rank_abun_all.csv",sep="")),sep=",", quote=F, row.names= TRUE, col.names= TRUE, fileEncoding="UTF-8" )
    write.table(rank_abun_matrix[,none_emp], file=(paste(name_out,"_rank_abun_detected.csv",sep="")),sep=",", quote=F, row.names= TRUE, col.names= TRUE, fileEncoding="UTF-8" )
    write.table(rank_abun_diff_matrix, file=(paste(name_out,"_rank_abun_diff_all.csv",sep="")),sep=",", quote=F, row.names= TRUE, col.names= TRUE, fileEncoding="UTF-8" )
    write.table(rank_abun_diff_matrix[,none_emp], file=(paste(name_out,"_rank_abun_diff_detected.csv",sep="")),sep=",", quote=F, row.names= TRUE, col.names= TRUE, fileEncoding="UTF-8" )

    write.table(abun_matrix, file=(paste(name_out,"_abun_all.csv",sep="")),sep=",", quote=F, row.names= TRUE, col.names= TRUE, fileEncoding="UTF-8" )
    write.table(abun_matrix[,none_emp], file=(paste(name_out,"_abun_detected.csv",sep="")),sep=",", quote=F, row.names= TRUE, col.names= TRUE, fileEncoding="UTF-8" )
    write.table(abun_diff_matrix, file=(paste(name_out,"_abun_diff_all.csv",sep="")),sep=",", quote=F, row.names= TRUE, col.names= TRUE, fileEncoding="UTF-8" )
    write.table(abun_diff_matrix[,none_emp], file=(paste(name_out,"_abun_diff_detected.csv",sep="")),sep=",", quote=F, row.names= TRUE, col.names= TRUE, fileEncoding="UTF-8" )


    #barplot for stat values
    library(ggplot2)
    library(reshape2)
    for (stat_name in rownames(Postive_matrix)){
    	data_in <- melt(Postive_matrix[stat_name,])
    	colnames(data_in)<-stat_name
    	data_in$samples <- rownames(data_in)
    	p <- ggplot(data_in,aes(x=samples,y=data_in[ ,stat_name])) + geom_col( position="dodge")+coord_flip()+ylab(stat_name)
    	ggsave(file=paste(name_out,"_",stat_name,"_barplot",".pdf",sep = ""),height=dim(data_in[1])*0.3,units ="cm",limitsize = FALSE)
    }

    #line plot for distance
    data_dis <- melt(dis_matrix_out)
    colnames(data_dis)<-c("Dist_type","Samples","Dist")
    data_dis$Samples <- factor(data_dis$Samples, levels=sort(rownames(data_in)), ordered=TRUE) #use data_in for leven and height is because the data_dis$Samples has duplications
    p <- ggplot(data_dis, aes(x=Samples, y=Dist, colour=Dist_type,group=Dist_type)) + geom_line()+coord_flip()
    ggsave(file=paste(name_out,"_distance_plot",".pdf",sep = ""),height=dim(data_in[1])*0.3,units ="cm",limitsize = FALSE)

    #lineplot for rank abundance different.
    data_rank_abun_diff <- melt(rank_abun_diff_matrix)
    colnames(data_rank_abun_diff)<-c("Rank","Samples","Dist")
    data_rank_abun_diff$Samples <- factor(data_rank_abun_diff$Samples, levels=sort(rownames(data_in)), ordered=TRUE)
    p <- ggplot(data_rank_abun_diff, aes(x=Samples, y=Dist, colour=Rank,group=Rank)) + geom_line()+coord_flip()
    ggsave(file=paste(name_out,"_rank_abund_diff_plot",".pdf",sep = ""),height=dim(data_in[1])*0.3,units ="cm",limitsize = FALSE)

}

##the JSD calculation on sorted relative abundance should work for any kind of strain prediction method with relative abundance reported
#filter the extreme low abundance
md_in[md_in<low_abun_filter]=0
md_in <- md_in[rowSums(md_in)>0,]

#max strain number selection
max_str_num <- max(dim(bench_in)[1],dim(md_in)[1])
bench_exp <- matrix(0,nrow=max_str_num,ncol=length(colnames(bench_in)))
md_in_exp <- matrix(0,nrow=max_str_num,ncol=length(colnames(bench_in)))
colnames(bench_exp)<-colnames(bench_in)
colnames(md_in_exp)<-colnames(bench_in)
rownames(bench_exp)<-paste("bench_GCA_str",1:max_str_num,sep="")
rownames(md_in_exp)<-paste("predict_GCA_str",1:max_str_num,sep="")


#sort by column
sort_by_column <- function(inter_samples,matrix_in,matrix_out){
    n_row <-dim(matrix_in)[1]
    for (sample_id in inter_samples){
        matrix_out[1:n_row,sample_id] <- sort(matrix_in[,sample_id],decreasing = TRUE)
    }
    return(matrix_out)
}
bench_exp <- sort_by_column(colnames(bench_in),bench_in,bench_exp)
inter_samples <- intersect(colnames(md_in),colnames(bench_in))
md_in_exp <- sort_by_column(inter_samples,md_in,md_in_exp)

#record the sorted abundance for generating stackplot
sort_abun_matrix <- rbind(bench_exp, md_in_exp)
write.table(sort_abun_matrix, file=(paste(name_out,"_sorted_RA_all.csv",sep="")),sep=",", quote=F, row.names= TRUE, col.names= TRUE, fileEncoding="UTF-8" )

#calculate sorted JSD for each sample 
library(phyloseq)
sJSD_matrix_out <- matrix(1,nrow=1,ncol=length(colnames(bench_in))) #in case some samples were missed by xtrain, then the distances must all be 1.
rownames(sJSD_matrix_out) <- "sorted_JSD"
colnames(sJSD_matrix_out) <- colnames(bench_in)
for (sample_id in inter_samples){
    otu_matrix <- cbind(bench_exp[,sample_id],md_in_exp[,sample_id])
    phylo_test <- phyloseq(otu_table(t(otu_matrix), taxa_are_rows = F))
    sJSD_matrix_out[1,sample_id] <- distance(phylo_test, "jsd")
}

#get sample with no result
none_emp <- intersect(colnames(md_in),colnames(bench_in))
em_rows <- setdiff(colnames(bench_in),colnames(md_in))

write.table(sJSD_matrix_out, file=(paste(name_out,"_sorted_JSD_all.csv",sep="")),sep=",", quote=F, row.names= TRUE, col.names= TRUE, fileEncoding="UTF-8" )
write.table(as.data.frame(sJSD_matrix_out)[,none_emp], file=(paste(name_out,"_sorted_JSD_detected.csv",sep="")),sep=",", quote=F, row.names= TRUE, col.names= TRUE, fileEncoding="UTF-8" )
write.table(em_rows, file=(paste(name_out,"_sample_missed.csv",sep="")),sep=",", quote=F, row.names= FALSE, col.names= FALSE, fileEncoding="UTF-8" )

#line plot 
library(ggplot2)
library(reshape2)
data_dis <- melt(sJSD_matrix_out)
colnames(data_dis)<-c("Dist_type","Samples","Dist")
p <- ggplot(data_dis, aes(x=Samples, y=Dist, colour=Dist_type,group=Dist_type)) + geom_line()+coord_flip()
ggsave(file=paste(name_out,"_sorted_JSD_plot",".pdf",sep = ""),height=dim(data_dis[1])*0.3,units ="cm",limitsize = FALSE)
