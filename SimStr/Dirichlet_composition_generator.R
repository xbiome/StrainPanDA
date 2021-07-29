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


#create log:202000403
#Objective: To generate the composition matrix following dirichlet distribution based on the number of strains, number of samples and the sum level of each sample
#Note1: to recognize the group name, use the second set of - in the given dataset name of the cross evaluation step. 
#Note2: it is hard coded to replace XT and EST in the anno name, other strings will caused the column facet fail.
#developed from the previous constrain_type_plot-jinhui.R improved plot details by Jinhui

#usage example: Rscript Dirichlet_composition_generator.R n_str=number_of_strains n_sample=number_of_samples sum_level=the_sum_of_depth_each_sample 
#the output will be the composition matrix

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


library(gtools) 

#check whether the file in is exist
if (!exists("n_str")) {
	stop("\nRscript Extract_candidate_strs.R n_str=n_str \nWarning: Usage: n_str para is not provied, please check the para \n\n")
}

if (!exists("n_sample")) {
  stop("\nRscript Extract_candidate_strs.R n_sample=n_sample \nWarning: Usage: n_sample para is not provied, please check the para \n\n")
}

if (!exists("sum_level")) {
  stop("\nRscript Extract_candidate_strs.R sum_level=sum_level \nWarning: Usage: sum_level para is not provied, please check the para \n\n")
}


#n_str=4
#n_sample=20
#sum_level=20
dist_matrix <- matrix(0,nrow=n_str,ncol=n_sample)
for (i in 1:n_sample){
tmp_dis <- rdirichlet(1,rep(1,n_str))
while(min(tmp_dis)<(1/sum_level)){
tmp_dis <- rdirichlet(1,rep(1,n_str))
}
dist_matrix[,i] <-tmp_dis
}
out_matrix <- dist_matrix
for (j in 1:n_str){
if(j!=n_str){
out_matrix[j,]<-round(dist_matrix[j,]*sum_level)
}else{
out_matrix[j,]<-sum_level-colSums(out_matrix[1:(j-1),])
}
}

write.table(out_matrix, file=(paste("dist_table_",n_str,"strs_",n_sample,"sample.csv",sep="")),sep=",", quote=F, row.names= FALSE, col.names= FALSE, fileEncoding="UTF-8" )