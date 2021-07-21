library(tidyverse)

data <- read_csv('data-raw/Faecalibacterium-prausnitzii-201911.counts.csv') %>%
  data.frame(row.names=1) %>%
  select(2:10, 1)

pangenome <- read_tsv('data-raw/panphlan_Faecalibacterium-prausnitzii-201911_pangenome.csv', col_names = F)

ans <- data.frame(matrix(0, 5,10))
colnames(ans) <- colnames(data)
ans[,1] <- c(0,0,0,0,1)
ans[,2] <- c(0,0.4,0.3,0.2,0.1)
ans[,3] <- c(0.1,0.1,0.1,0.1,0.6)
ans[,4] <- c(0.1,0.2,0.3,0.4,0)
ans[,5] <- c(0.2,0.2,0.2,0.2,0.2)
ans[,6] <- c(0.2,0.2,0,0.4,0.2)
ans[,7] <- c(1,0,0,0,0)
ans[,8] <- c(0,1,0,0,0)
ans[,9] <- c(0,0,1,0,0)
ans[,10] <- c(0,0,0,1,0)
rownames(ans) <- c('GCF_002549755', 'GCF_002549975', 'GCF_003287415', 'GCF_000162015', 'GCF_003433905')
fprausnitzii <- list(data=data.frame(data), ans=ans, pangenome=data.frame(pangenome))


usethis::use_data(fprausnitzii, overwrite = TRUE)
