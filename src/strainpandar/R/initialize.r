##' @title strainpanda.init
##'
##' @import dplyr
##' @import foreach
##' @importFrom NNLM nnlm
##' @importFrom vegan vegdist
##' @param obj strainpandar object generated using preprocess
##' @param abundance.thre threshold for relative abundance (default: 0.05)
##' @param fil.thre threshold for filtering abnormal gene families (default: 2)
##' @param cutree.h the height at which the clustering tree is cut (defualt: 0.1)
##' @param ncpu number of processes running in parallel (default: 1)
##' @author Xbiome
##' @description select most possible genomes from the reference database
##' @export
strainpanda.init <- function(obj, ncpu=1, abundance.thre=0.05, fil.thre=2, cutree.h=0.1){

  D <- obj$data

  ## cluster reference to remove too close ones
  clustering <- hclust(vegan::vegdist(t(obj$reference), method = "jaccard"))
  sel <- cutree(clustering, h=cutree.h) %>% data.frame(c=.) %>%
    tibble::rownames_to_column("genome") %>%
    group_by(c) %>% mutate(n=1:n()) %>%
    filter(n==1) %>% pull(1)
  P <- obj$reference[, sel]
  res <- solveStrainSampleRegression(D, P, ncpu, fil.thre=fil.thre, check.x=F)
  rownames(res) <- colnames(P)
  colnames(res) <- colnames(D)
  ## second round
  sel <- names(which(rowSums(res[,] >= abundance.thre)>0))
  P <- obj$reference[, sel]
  res <- solveStrainSampleRegression(D, P, ncpu, fil.thre=fil.thre, check.x=F)
  rownames(res) <- colnames(P)
  colnames(res) <- colnames(D)
  res
}
