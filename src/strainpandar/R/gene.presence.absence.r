
##' @title gene.presence
##'
##' @importFrom pracma findpeaks
##' @param p one column from the predicted P matrix
##' @param adjust adjust the threshold between start of the presence peak and the presence peak, must be [0,1]. Smaller value means more relaxed. [Default: 0]
##' @param n.gene.expected expected number of genes
##' @author Xbiome
##' @description Determine whether a gene is present or not using a peak finding approach
##' @export
gene.presence <- function(p, n.gene.expected, adjust=0){
  p.origin <- p
  #p[p < (hist(p.origin))$breaks[2]] <- NA
  if(all(p<=0)){
    ## this strain is not present here
    return(list(binary=p, scores=p))
  }
  bw <- tryCatch(d <- density(p, na.rm=TRUE),
           error = function(c) NA
  )
  if(length(bw)==1){
    ### some times the default "nrd0" method doesn't works
    d <- density(p, na.rm=TRUE, bw="ucv")
  }

  ## find all peaks
  peaks <- pracma::findpeaks(d$y)
  ## cumulate to the right of each peak
  cutoffs <- d$x[sort(peaks[,2])]
  ## number of genes for each cutoff
  n.genes <- vapply(cutoffs, function(x) sum(p>x), 1)
  ## closest number of genes to the expected gene family number
  cutoff <- cutoffs[which.min(abs(n.genes-n.gene.expected))]
  cutoff.adjusted <- sum(c(cutoff, 0) * c(1-adjust, adjust))

  ## TODO: add assertion if all peaks are far from the min number of genes?
  scores <- vapply(p-cutoff.adjusted, function(x) min(x, 0), 0.1)
  scores <- (scores-min(scores))/(max(scores)-min(scores))

  ## TODO: check which one works better (cutoff or cutoff.adjusted)

  list(binary=(p.origin > cutoff.adjusted)*1, scores=scores)
}
