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

cDIR=getwd()
f_path=paste(cDIR,"/abund_prof_merge",sep = "")
ind_pattern <- "profile*"
out_name=paste(f_path,"/Est_merged_matrix.tsv",sep = "")
Est_merge <- function(f_path, out_name,ind_pattern){
  file_list <- list.files(path=f_path,pattern=ind_pattern)
  Count_Est=0
  for (Est_gene_table in file_list) {
    print(Est_gene_table)
    Est_gene <-read.table(paste(f_path,Est_gene_table,sep = "/"),header = T,sep="\t",comment.char = "")
    #Sample_ID <- colnames(Est_gene)[2]
    #Est_gene$Sample <- rep(Sample_ID,length(rownames(Est_gene)))
    if(Count_Est==0){
      Est_gene_merge <- Est_gene
    } else{
      Est_gene_merge <- merge(Est_gene_merge, Est_gene, by="OTU")
    }
    Count_Est=+1
  }
  write.table(Est_gene_merge, file=out_name,sep="\t", quote=F, row.names= FALSE )
  return(Est_gene_merge)
}

Est_gene_merge <- Est_merge(f_path, out_name,ind_pattern)
