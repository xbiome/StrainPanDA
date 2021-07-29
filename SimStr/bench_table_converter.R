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

#usage example: Rscript bench_table_converter.R sim_table=benchmark_comp.csv add_name=_seqerr_r rep_add_time=3 name_out=bench_Ecoli98_seqerr_all
#the output will be a str profile table from the bench table, replicate names and data tye will be added and expanded

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
if (!exists("sim_table")) {
	stop("\nRscript bench_table_converter.R sim_table=benchmark_comp_WGS.csv \nWarning: Usage: sim_table para is not given. \n\n")
}

if (!file.exists(sim_table)) {
    stop(paste("\nRscript bench_table_converter.R sim_table=",sim_table," \nWarning: Usage: simulation bench file is not exist, please check the path. \n\n",sep=""))
}

#sim_table <- "benchmark_comp.csv"
#add_name <- "_seqerr_r"
#rep_add_time <- 3
#name_out <- "bench_Ecoli98_seqerr_all"

#get bench data
SIM_in <- read.csv(sim_table)
SIM_matrix <- SIM_in[-1,4:dim(SIM_in)[2]]
rownames(SIM_matrix) <- SIM_in[-1,2]

#add the data type and replicate into sample names
rep_add=0
col_n <-colnames(SIM_matrix)
if (rep_add_time>1){
    col_n_m <- paste(col_n,add_name,rep_add,sep="")
    col_n_m_l <- col_n_m
    rep_add=rep_add+1
    for (i in 2:rep_add_time){
        col_n_m <- paste(col_n,add_name,rep_add,sep="")
        col_n_m_l <-c(col_n_m_l,col_n_m)
        rep_add=rep_add+1
    }
} else {
    col_n_m <- paste(col_n,add_name,"0",sep="")
    col_n_m_l <- col_n_m
}

#expand the data dimension
if (rep_add_time>1){
    SIM_matrix_f <- SIM_matrix
    for (i in 2:rep_add_time){
        SIM_matrix_f <- cbind(SIM_matrix_f,SIM_matrix)
    }
    colnames(SIM_matrix_f) <- col_n_m_l
} else {
    SIM_matrix_f <- SIM_matrix
    colnames(SIM_matrix_f) <- col_n_m_l
}

SIM_matrix_out <- SIM_matrix_f
#turn into relative abund
c_s <- apply(SIM_matrix_f,2,sum)
for(i in names(c_s)){
  SIM_matrix_out[,i] <- SIM_matrix_f[,i]/c_s[i]
}

write.table(SIM_matrix_out, file=(paste(name_out,"_str_prof.csv",sep="")),sep=",", quote=F, row.names= TRUE, col.names= TRUE, fileEncoding="UTF-8" )

#然后再写一个比较两个这种matrix的脚本（BC距离），如果方法输出比模拟输出多，额外的行也一样算进去距离里。



