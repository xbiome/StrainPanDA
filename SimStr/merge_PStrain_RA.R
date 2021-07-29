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


#must be in the PS_out folder and have the sample_list.txt file in the parent folder

cDIR=getwd()
sample_list="../sample_list.txt"
samples <- unlist(read.table(sample_list, stringsAsFactors =FALSE))
out_name <- "PStrian_merged_RA.csv"
species_list =""
for (i in 1:length(samples)) {
  sample_name <- samples[i]
  result_path <- paste(cDIR,"/",sample_name,"/result/strain_RA.txt",sep = "")
  if (!file.exists(result_path)) {
    print(paste("\nThe following result (",result_path,") was not generated \n\n",sep=""))
  } else{
    print(result_path)
    tmp_matr <- read.table(result_path, stringsAsFactors =FALSE)
    colnames(tmp_matr)[4]<-sample_name
    for (species in unique(tmp_matr[,1])){
      if(!exists(paste("result_matr",species,sep="_"))){
        assign(paste("result_matr",species,sep="_"),tmp_matr[which(tmp_matr[,1]==species),c(3,4)])
        if(species_list ==""){
          species_list = paste("result_matr",species,sep="_")
        }else{
          species_list = c(species_list,paste("result_matr",species,sep="_"))
        }
      }else{
        assign(paste("result_matr",species,sep="_"), merge(get(paste("result_matr",species,sep="_")),tmp_matr[which(tmp_matr[,1]==species),c(3,4)],by=1,all = TRUE))
      }
    }  
  }
}


for (species_matrix in species_list){
  result_matr=get(species_matrix)
  result_matr[is.na(result_matr)]=0
  write.csv(result_matr, file=paste(species_matrix,"_",out_name,sep=""), quote=F, row.names= FALSE )
}