##' @title findNMFRank
##'
##' @import NMF
##' @import dplyr
##' @importFrom vegan vegdist
##' @param obj strainpandar object (output of `preprocess`)
##' @param jc.threshold jaccard distance threshold
##' @param abun.threshold mean abundance threshold
##' @param gene.family.min gene family number used for gene presence detection
##' @param presence.thre gene presence threshold
##' @param ncpu number of processes running in parallel (default: 1)
##' @param max.strain the maximum number of strains to try
##' @param seed used for nmf
##' @param ... other parameters for nmf
##' @export
##' @author Xbiome
findNMFRank <- function(obj, gene.family.min=3000, presence.thre=0.5,
                        jc.threshold=0.1, abun.threshold=0.1,
                        ncpu=1, max.strain=12, seed=123456, ...){
  est <- obj
  if(class(obj)!="NMF.rank"){
    ## ensure number of samples > max.strain
    max.strain <- min(max.strain, ncol(obj$data))
    est <- NMF::nmf(obj$data, 1:max.strain, nrun=10, seed=seed, method = "snmf/r",
               .pbackend='par', .options=paste0('p', ncpu), ...)
  }

  tmp <- foreach(r=2:(length(est$fit)-1), .combine = "rbind") %do%{
    s <- scoef(est$fit[[r]])
    p_est <- basis(est$fit[[r]])
    ## 1. jaccard distance > 0.1 for all strains
    ## crude estimation of gene presence-absence matrix
    p <- apply(p_est, 2, function(x) gene.presence(x, gene.family.min, adjust = presence.thre)$binary)
    jc.dist <- c(vegdist(t(p), method="jaccard"))
    ## 2. mean abundance
    abun.min <- min(apply(s, 1, function(x) mean(x)))
    ## 3. prevalence
    prev.min <- min(apply(s, 1, function(x) sum(x>0.)/length(x) ))
    ## 4. # of gene families must not be to large or small
    n.gene <- colSums(p)
    data.frame(rank=r, jc.dist=min(jc.dist), abun.min=abun.min, prev.min=prev.min, gene.min=min(n.gene), gene.max=max(n.gene))
  }

  ## todo: handle case when there is only 1 strain
  #tmp
  r <- (filter(tmp, jc.dist>jc.threshold &
                 abun.min>abun.threshold &
                 gene.min>gene.family.min*0.5) %>%
          pull(rank) %>% rev() )[1]

  list(solutions=est, rank=r)
}
