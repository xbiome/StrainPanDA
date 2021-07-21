##' @title strain.decompose
##'
##' @import NMF
##' @importFrom MASS ginv
##' @importFrom vegan vegdist
##' @param obj strainpandar object (output of `preprocess`)
##' @param rank the rank of NMF (number of strains)
##' @param presence.thre in iterative process, for a gene to be considered as present, this threshold is used as `adjsut` for the function `gene.presence` (default: 0.5)
##' @param ncpu number of processes running in parallel (default: 1)
##' @param seed random seed for NMF
##' @param debug print debugging information (default: FALSE)
##' @author Xbiome
##' @description Decovolute strain information based on pangenome coverage information
##' @export
strain.decompose <- function(obj, rank=NULL,
                             ## parameters for determining the number of strains
                             jc.threshold=0.1, abun.threshold=0.1, max.strain=12,
                             presence.thre=0.5,
                             ## run controls
                             ncpu=1, seed=123456, debug=FALSE){
  library(NMF)
  if (ncol(obj$data) < max.strain) {
    message(paste0("Setting the maximum number of strains to ", ncol(obj$data), " (number of samples)..."))
    max.strain <- ncol(obj$data)
  }
  if (is.null(rank)){
    nmf.rank <- findNMFRank(obj, jc.threshold=jc.threshold, abun.threshold=abun.threshold,
                        ncpu=ncpu, max.strain=max.strain, seed=seed)
    rank <- nmf.rank$rank
  
    message(paste0("\nUsing NMF with rank ", rank, "..."))

    S <- scoef(nmf.rank$solutions$fit[[rank]])
    P_est <- basis(nmf.rank$solutions$fit[[rank]])
  }
  else {
    nmf.rank <- NMF::nmf(obj$data, rank, nrun=10, seed=seed, method = "snmf/r",
                    .pbackend='par', .options=paste0('p', ncpu))
    S <- scoef(nmf.rank)
    P_est <- basis(nmf.rank)
  }
  P <- apply(P_est, 2, function(x) gene.presence(x, obj$min.genefamily, presence.thre)$binary)

  rownames(S) <- colnames(P) <- colnames(P_est) <- paste0('strain',1:rank)
  colnames(S)  <- colnames(obj$data)
  rownames(P) <- rownames(P_est) <- rownames(obj$data)

  list(S=S, P=P, P_est=P_est, nmf.rank=nmf.rank)
}
