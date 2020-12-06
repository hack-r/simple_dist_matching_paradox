mydist <- function(a, b){
  return(sqrt(sum((a-b)^2)))
}

euclidean.distfun <- function(mydat){
  mytreat <- mydat[mydat$treat==1,]
  mycontrol <- mydat[mydat$treat==0,]
  distmat <- matrix(NA, nrow=nrow(mytreat), ncol=nrow(mycontrol))
  rownames(distmat) <- rownames(mytreat)
  colnames(distmat) <- rownames(mycontrol)
  for(i in 1:nrow(mytreat)){
    a <- mytreat[i, c("a")]
    for(j in 1:nrow(mycontrol)){
      b <- mycontrol[j, c("a")]
      distmat[i, j] <- mydist(a, b)
    }
  }
  return(distmat)
}
