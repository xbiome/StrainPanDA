##' @title solveStrainSampleRegression
##'
##' @import foreach
##' @import dplyr
##' @importFrom NNLM nnlm
##' @param D a normalized count matrix (output of `preprocess`)
##' @param P current gene-strain matrix
##' @param ncpu number of processes running in parallel (default: 1)
##' @param fil.thre threshold for filtering abnormal gene families (default: 2)
##' @param ... other parameters for nnlm
##' @export
##' @author Xbiome
solveStrainSampleRegression <- function(D, P, ncpu=1, fil.thre=2, ...){
  N <- ncol(P)  # number of strains
  S <- ncol(D)  # number of samples

  idx <- colSums(P) != 0
  tmp <- foreach(i=1:S, .combine = cbind) %do% {
    ## some virtual strains might not be present
    coefs <- rep(0, N)
    grouping <- apply(P, 1, function(x) paste0(x, collapse = ""))
    reg.dat <- data.frame(Y=(D[,i]), P,  g=grouping) %>%
      group_by(g) %>%
      filter(n()>30) %>%
      filter( abs(Y-median(Y))/(mad(Y)+1e-15 ) < fil.thre) %>%
      ungroup() %>%
      select(-g)

    coefs[idx] <- coef(NNLM::nnlm(y=as.matrix(reg.dat[,1]),
                    x=as.matrix(reg.dat[,-1]),
                    n.threads = ncpu,
                    alpha=c(0,0,1),
                    ...))[,1]

    if(all(coefs == 0)){
      warning(paste("Bad status for sample no.", i,"\n"))
    }else{
      coefs <- coefs/sum(coefs)
      coefs[coefs < 0.01] <- 0
    }
    coefs
  }
  tmp
}

