#!/usr/bin/env Rscript


library(strainpandar)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(dplyr)

## replace with getopt
args <- commandArgs(trailingOnly=TRUE)

counts.file <- args[1]
pangenome.path <- args[2]
output <- args[3]
ncpu <- as.numeric(args[4])
rank <- as.numeric(args[5])
ko.profile <- FALSE

pangenome.file <- list.files(pangenome.path, "*pangenome.csv", full.names = TRUE)

if(length(pangenome.file) == 0){
  stop("Cannot locate pangenome file...")
}

profile <- read.csv(counts.file, row.names=1)

profile.preprocessed <- preprocess(profile, pangenome.file = pangenome.file)	

if(ncol(profile.preprocessed$data)<5){
  stop("Less than 5 samples left after preprocessing...")
}

res <- strain.decompose(profile.preprocessed, ncpu=ncpu,rank=rank)

## strain-sample plot
write.table(res$S, file = paste0(output, ".strain_sample.csv"), sep=",",
            quote=F, row.names = TRUE, col.names = TRUE)

p <- melt(res$S) %>%
  ggplot(aes(x=Var2, y=value, fill=Var1)) +
  geom_bar(stat="identity") +
  labs(x=NULL, y="Relative abundances") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        legend.title = element_blank())
ggsave(p, filename = paste0(output, ".strain_sample.pdf"), height = 10, width = 10)

## gene family-strain plot
refs <- profile.preprocessed$reference

p.filtered <- res$P[, colSums(res$P)!=0]
heatmap <- merge( p.filtered, refs, by=0, all.y=TRUE)
heatmap <- data.frame(heatmap, row.names = 1)
heatmap[is.na(heatmap)] <- 0

write.table(heatmap, file = paste0(output, ".genefamily_strain.csv"), sep=",",
            quote=F, row.names = TRUE, col.names = TRUE)

# prokka <- 
#   read.table(list.files(pangenome.path, "*prokka.gff", full.names=TRUE, recursive = T),
#              quote = "$$",
#              fill=TRUE,
#              sep="\t")[,c(1,9)] %>% 
#     filter(V9!="") %>% 
#     mutate(V9=gsub( ".*product=", "", V9)) %>% 
#     select(V1, prokka=V9)
# 
# panphlan.anno <- read.table(list.files(pangenome.path, "*annotations.csv", full.names=TRUE),
#            sep="\t", head=T, quote="$$")[, -c(2,3)]
# panphlan.anno$map <- gsub("\\.[0-9]+", "", panphlan.anno$Protein_id)
# mapping <- read.table(list.files(pangenome.path, "panphlan.annotation.tsv", full.names = TRUE), 
#                       sep="\t", head=F)
# panphlan.anno <- 
#   merge(panphlan.anno, mapping, by.x="map", by.y="V2", all.x=TRUE) %>% 
#   select(-V1, -V3, -map) %>% 
#   select(everything(), KO=V4) %>% 
#   mutate(KO=ifelse(KO=="", NA, as.character(KO))) 
# 
# p.filtered[rowSums(p.filtered)==1 | rowSums(p.filtered)==nrow(p.filtered)-1, ] %>%
#   merge(panphlan.anno, by.x=0, by.y=1, all.x=TRUE) %>% 
#   select(gene_family_id=Row.names, everything()) %>% head
# 
#   write.table(paste0(output, ".strain_specific_family_annotation.tsv"),
#               row.names = F,
#               sep="\t", quote=F)


if( nrow(heatmap) < 60000 ){
    pheatmap((heatmap>0)*1,filename = paste0(output, ".genefamily_strain.pdf"), show_rownames = FALSE, fontsize = 5)
}

if(ko.profile) {
    ko.map <- list.files(pangenome.path, "*genefamily_knumber.txt", full.names = TRUE)
    if(length(ko.map)==0){
      stop("No gene family-KO mapping file found...")
    }
    ko <- read.table(ko.map)

    for(i in 1:ncol(p.filtered)){
      tmp <- select(data.frame(p.filtered), i) %>%
      tibble::rownames_to_column("id")
      tmp[tmp[,2]!=0, ] %>% select(1) %>% merge(ko, by=1) %>%
      write.table(paste0(output, ".strain", i, ".ko.txt"), quote=F, row.names = F, col.names = F)
    }
}else{
    write.table(NULL, paste0(output, ".strain.dummy.ko.txt"), quote=F, row.names = F, col.names = F)
}

saveRDS(res, paste0(output, ".rds"))


#find neighbors
#rm extra strings in the colnames
name_list <- colnames(res$S)

for(i in 1:length(name_list)){
  n_list <- strsplit(name_list[i], split = "_", fixed = TRUE)[[1]]
  name_list[i] <- paste0(n_list[1:(length(n_list)-1)],collapse = "_")
}

m_out <- res$S
colnames(m_out)<-name_list


#library(tidyverse)
## True genome profie
pangenome <- read.table(pangenome.file, head = F) %>% dplyr::select(1,3) %>% distinct() %>%  mutate(presence=1) %>%
  dcast(V1~V3) %>%
  data.frame(row.names = 1)
pangenome[is.na(pangenome)] <- 0
#head(pangenome)
merged <- merge(pangenome, res$P, all=TRUE, by=0)
merged[is.na(merged)] <- 0
merged_m <- merged[,-1]
rownames(merged_m) <- merged[,1]
#generate tree
p_prof <- t(merged_m)
#out_dist <- dist(p_prof,method="euclidean")
out_dist <- vegan::vegdist(p_prof, method = "jaccard")

#get neighbor
strs_l<-colnames(res$P)
dist_m<-as.matrix(out_dist)
write.table(dist_m, file=(paste(output,"_all_dis.csv",sep="")),sep=",", quote=F, row.names= TRUE, col.names= TRUE, fileEncoding="UTF-8" )

neighbor_l <- names(which(dist_m[,1]==min(dist_m[-1,1])))
for(i in 2:dim(dist_m)[2]){
    neighbor_i <- names(which(dist_m[,i]==min(dist_m[-i,i])))[1]
    neighbor_l <- c(neighbor_l,neighbor_i)
}
names(neighbor_l) <- colnames(dist_m)
write.table(neighbor_l, file=(paste(output,"_all_neighbor.csv",sep="")),sep=",", quote=F, row.names= TRUE, col.names= FALSE, fileEncoding="UTF-8" )
#这里要输出两个表格，一个是所有的，一个是只有str_n的，作为注释的
dist_nostr <- dist_m
dist_nostr <- dist_m[setdiff(rownames(dist_m),strs_l),]
neighbor_str <- names(which(dist_nostr[,strs_l[1]]==min(dist_nostr[,strs_l[1]])))
if (length(strs_l)>1){
    for(j in 2:length(strs_l)){
        neighbor_j <- names(which(dist_nostr[,strs_l[j]]==min(dist_nostr[,strs_l[j]])))[1]
        neighbor_str <- c(neighbor_str,neighbor_j)
    }
}
names(neighbor_str) <- strs_l
write.table(neighbor_str, file=(paste(output,"_str_neighbor.csv",sep="")),sep=",", quote=F, row.names= TRUE, col.names= FALSE, fileEncoding="UTF-8" )

#correct the names of strains with predicted neighbors
rn_m_dup <- neighbor_str[duplicated(neighbor_str)]
#check duplicated rn_m_out, merge the rows with same strain names.
if(sum(duplicated(neighbor_str))==0){
    m_out_anno <- m_out
} else {
    if(sum(!duplicated(neighbor_str))>1){
        m_out_anno <- m_out[names(neighbor_str[!duplicated(neighbor_str)]),]
        for(dup_n in 1:length(rn_m_dup)){
            dup_i <- names(which(neighbor_str[!duplicated(neighbor_str)]==rn_m_dup[dup_n]))
            print(names(rn_m_dup[dup_n]))
            m_out_anno[dup_i,] <- m_out_anno[dup_i,]+m_out[names(rn_m_dup[dup_n]),]
        }
    } else{
        m_out_anno <- matrix(1,nrow=1,ncol=dim(m_out)[2])
        colnames(m_out_anno) <- colnames(m_out)
        rownames(m_out_anno) <- "strain1"
    }
}
rownames(m_out_anno) <- names(neighbor_str[!duplicated(neighbor_str)])
write.table(m_out_anno, file=(paste(output,"_str_merged_prof.csv",sep="")),sep=",", quote=F, row.names= TRUE, col.names= TRUE, fileEncoding="UTF-8" )
rownames(m_out_anno) <- neighbor_str[!duplicated(neighbor_str)]
write.table(m_out_anno, file=(paste(output,"_str_anno_prof.csv",sep="")),sep=",", quote=F, row.names= TRUE, col.names= TRUE, fileEncoding="UTF-8" )

p <- melt(m_out_anno) %>%
  ggplot(aes(x=Var2, y=value, fill=Var1)) +
  geom_bar(stat="identity") +
  labs(x=NULL, y="Relative abundances") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        legend.title = element_blank())
ggsave(p, filename = paste0(output, ".anno_strain_sample.pdf"), height = 10, width = 10)
