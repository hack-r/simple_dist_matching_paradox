
###################################################################
## Functions for analysis below
## function from C:\Users\Richard Nielsen\Desktop\Papers\pscore paradox\aid shocks analysis\spacegraph function w replacement.R

mdistfun <- function(formula, data){
  ## I think this function has a dependence on "tvar"

  ## a mahalanobis function
  ## myMH function from http://www.stat.lsa.umich.edu/~bbh/optmatch/doc/mahalanobisMatching.pdf

  ## the second stopifnot condition was in ben hansen's code
  ## but the code seems to work fine for even one variable
  ## and I have samples with just one covariate.
  myMH <- function(Tnms, Cnms, inv.cov, data) {
   stopifnot(!is.null(dimnames(inv.cov)[[1]]),# dim(inv.cov)[1] > 1,
   all.equal(dimnames(inv.cov)[[1]], dimnames(inv.cov)[[2]]),
   all(dimnames(inv.cov)[[1]] %in% names(data)))
   covars <- dimnames(inv.cov)[[1]]
   xdiffs <- as.matrix(data[Tnms, covars])
   xdiffs <- xdiffs - as.matrix(data[Cnms, covars])
   rowSums((xdiffs %*% inv.cov) * xdiffs)
  }
    ## pull out the covariates
  X1 <- data
  X2 <- data.frame(model.frame(formula, data))
  tvar <- as.character(formula)[2]
    ## just the covariates
  X3 <- subset(data.frame(model.matrix(formula,data)),select=-1)
  mvrs <- colnames(X3)
    ## covariates and treatment variable
  X <- cbind(X2[,1],X3)
  names(X) <- c(tvar,mvrs)

  ## from http://tolstoy.newcastle.edu.au/R/help/05/05/3996.html
  matrix.rank <- function(A, eps=sqrt(.Machine$double.eps)){
    sv. <- abs(svd(A)$d)
    sum((sv./max(sv.))>eps)
  }
  #if(dim(X)[2]!=matrix.rank(X)){
  #  #print("MDM failed: The matrix is not of full rank.")
  #}
  #stopifnot(dim(X)[2]==matrix.rank(X), silent=T)

    ## if it can't invert, this is where it will fail
  try(icv <- solve(cov(subset(X,select=mvrs))))
  stopifnot(exists("icv"))
  trtnms <- row.names(X)[X[[tvar]]==1]
  ctlnms <- row.names(X)[X[[tvar]]==0]
    ## use the data X3 that has just the covariates

    ## run the main internal function
  mydist <- outer(trtnms, ctlnms, FUN = myMH, inv.cov = icv, data = X3)
  stopifnot(exists("mydist"))
  dimnames(mydist) <- list(trtnms, ctlnms)
  return(mydist)
}


psm <- function(formula, data, ratio=1, caliper=0, linear.pscore=F){
    ## estimate the propensity score
  if(nrow(na.omit(data))!=nrow(data)){
    print("Missing values in propensity score model.  Observations with missing data ignored.")
    data <- na.omit(data)
  }
  mymod <- (glm(formula, data, family=(binomial(link="logit"))))
  
  if(mymod$converged==F){break}
  ps <- mymod$fitted.values
  if(linear.pscore==T){
    ps <- log(ps/(1-ps))
  }
    ## set the number of matches, 1-to-k
  k <- ratio

  tvar <-  as.character(as.formula(formula))[2]
  Tnms <- row.names(data)[data[[tvar]]==1]
  Cnms <- row.names(data)[data[[tvar]]==0]

    ## make some holders
  match.matrix <- matrix(NA, length(Tnms), k)
  rownames(match.matrix) <- Tnms
  dist.matrix <- match.matrix

  t0 <- Sys.time()
  ixnay <- c()
    ## loop over the match ratio
  for(j in 1:k){
      ## loop over the treated units
    for(i in 1:length(Tnms)){
      #if(i%%100==0){
      #  print(paste("Match",j,"of",k,", T",i,"of",length(Tnms),":",round(Sys.time()-t0,2)))
      #}
        ## define a vector of the Cnms that are still available
      #elig <- Cnms[(Cnms%in%ixnay)==F]
      elig <- Cnms
      #if(i%%500==0){print(paste("treated unit",i,"of",length(Tnms)))}
        ## make the distance vector for treated t and all the eligable
        ## controls
      x <-abs(ps[Tnms[i]] - ps[elig])
        ## capture the colnames of the ones that have
        ## minimum distances and take the first
      min.x <- min(x)
      m <- (elig[(x==min.x)==T])[1]
        ## add the things to the output matrices
      match.matrix[i,j] <- m
      dist.matrix[i,j] <- min.x
        ## add the new match to ixnay
      ixnay <- c(ixnay, m)
        ## stop the iterations if there are no more control units
      #if(length(ixnay)==length(Cnms)){break}
    }
    ## stop the iterations if there are no more control units
    if(length(ixnay)==length(Cnms)){break}
  }


  ## make the matched data
    ## first, select all the treated obs
  mh.data <- data[rownames(match.matrix)[is.na(match.matrix[,1])==F],]
    ## then add the control obs
  for(i in 1:ncol(match.matrix)){
    mh.data <- rbind(mh.data,data[na.omit(match.matrix[,i]),])
  }


  ## add distances -- mean of all distances for a matched group
  dd <- apply(dist.matrix, MAR=1, mean, na.rm=T)
  distance <- matrix(NA,nrow(mh.data),1)
  rownames(distance) <- rownames(mh.data)
  distance[names(na.omit(dist.matrix[,1])),] <- dd[names(na.omit(dist.matrix[,1]))]
  for(i in 1:ncol(dist.matrix)){
    names(dd) <- match.matrix[,i]
    distance[na.omit(match.matrix[,i]),] <- dd[na.omit(match.matrix[,i])]
  }
  mh.data$distance <- distance

    ## add weights to the data
  weights <- matrix(NA,nrow(mh.data),1)
  rownames(weights) <- rownames(mh.data)
  weights[names(na.omit(match.matrix[,1])),]<-1
    ## calculate the weights
  countf <- function(x){sum(is.na(x))}
  wts <- 1/(k-apply(match.matrix, MAR=1, countf))
    ## fill in the weights
    ## this has me all mixed up, but I think it works
  for(i in 1:ncol(dist.matrix)){
    names(wts) <- match.matrix[,i]
    weights[na.omit(match.matrix[,i]),] <-wts[na.omit(match.matrix[,i])]
  }
  mh.data$weights <- weights

  ## Begin calipering
  match.obj <- list("match.matrix"=match.matrix, "dist.matrix"=dist.matrix,
              "match.dat"=mh.data)
  rc <- randomCaliper(match.obj, tvar, rawCP=rawCP, caliper=caliper)

  ## Need to make calipered match matrix, data, etc.
  keepobs <- rownames(rc$caliperdat)
  match.matrix <- as.matrix(match.matrix[which(rownames(match.matrix) %in% keepobs),])
  dist.matrix <- as.matrix(dist.matrix[which(rownames(dist.matrix) %in% keepobs),])
  
 
    ## RETURN
  out <- list("match.matrix"=match.matrix, "dist.matrix"=dist.matrix,
              "match.dat"=rc$caliperdat,call=match.call(),
              N=rc$N)
  class(out) <- "psm.match"
  return(out)
}


## This function takes a "match.obj" for either myMH or myPS and
## applies calipers.

  ## 1-to-1 matching ONLY at the moment.
randomCaliper <- function(match.obj, treatment=treatment,rawCP=rawCP, caliper=0){
  mdat <- match.obj$match.dat
     ## pull out the solution
  solution <- match.obj
     ## make the dataset
  t.units <- rownames(solution$match.matrix)
  c.units1 <- solution$match.matrix[,1]
    ## this makes sure the names in the match.matrix get matched to
    ## the data rownames correctly
  matchfind <- function(x,vec){
    if(is.na(x)){return(NA)}
    if(!is.na(x)){return(which(vec==x))}
  }
  t.rows <- apply(as.matrix(t.units),MAR=1,FUN=matchfind,vec=rownames(mdat))
  if(is.list(t.rows)){
    t.rows <- unlist(t.rows)
  }
  c.rows <- apply(as.matrix(c.units1),MAR=1,FUN=matchfind,vec=rownames(mdat))
    ## Then, make sure to only include treated units that have at
    ## least one match
  matched <- which(is.na(solution$match.matrix[,1])==F)
    ## then combine them to make the treated data matrix
  tdat <- mdat[t.rows[matched],]
    ## the cdat1 will be the same, because it has to be 1-to-1
  cdat1 <- mdat[c.rows[matched],]
  psdiff <- tdat$distance  ## diff from pscore

    ## order the datasets so that they are in order of best
    ## to worst matches
  tdat <- tdat[order(psdiff),]
  cdat1 <- cdat1[order(psdiff),]

    ## calculate the ATT, L1s, etc

  #if(exists(psrunsR[i])){
    if(caliper > (nrow(tdat) - 1)){
      caliper <- (nrow(tdat) - 1)
      print("Caliper exceeds number of matches.  Automatically reset to the maximum possible")
    }
    caliperdat <-  na.omit(rbind(tdat[1:(nrow(tdat)-caliper),],
                     cdat1[1:(nrow(cdat1)-caliper),]))
      ## calculate the ATT, L1, etc
      ## this now calls "my.imbalance" which takes out the chi-sq
      ## calculation that was leading to bugs.
    #L1 <- (imbalance(caliperdat[[treatment]],caliperdat,drop=c(treatment,"distance","weights"),breaks=rawCP))$L1$L1
    #L1 <- (my.imbalance(caliperdat[[treatment]],caliperdat,drop=c(treatment,"distance","weights"),breaks=rawCP))$L1$L1
    N <- table(caliperdat[[treatment]])
    matched.rows <- rownames(caliperdat)

  return(list(N=N,matched.rows=matched.rows, caliperdat=caliperdat))

}

## A mahalanobis discrepancy function
## adist is a matrix of mah distances produced by mdistfun()

mdisc <- function(caliperdat, adist=alldist, wt = NULL){
  if(is.null(wt)){wt <- rep(1,nrow(caliperdat))}
  matched.names <- rownames(caliperdat)[wt>0]
  mdis <- adist[rownames(adist)[rownames(adist) %in% matched.names],
          colnames(adist)[colnames(adist) %in% matched.names]]
  if(is.matrix(mdis)){
    try(minT <- apply(mdis,MAR=1,min), silent=T)
    ## If the vector is too big for apply
    if(exists("minT")==F){
      minT <- rep(NA,nrow(mdis))
      for(jj in 1:nrow(mdis)){minT[jj] <- min(mdis[jj,])}
    }
    ## 
    try(minC <- apply(mdis,MAR=2,min), silent=T)
    if(exists("minC")==F){
      minC <- rep(NA,ncol(mdis))
      for(jj in 1:ncol(mdis)){minC[jj] <- min(mdis[,jj])}
    }
  } else { 
    minT <- min(mdis)
    minC <- min(mdis)
  }    
  return(mean(c(minT,minC)))
}


## a second version that doesn't need the distances calculated outside
mdisc2 <- function(caliperdat, alldat,mvars,tvar, wt = NULL){
  ff <- as.formula(paste(tvar,"~",paste(mvars,collapse="+")))
  try(adist <- mdistfun(ff,alldat), silent=T)
  if(exists("adist")==F){
    adist <- mdistfun2(ff,alldat)
  }
  if(is.null(wt)){wt <- rep(1,nrow(caliperdat))}
  matched.names <- rownames(caliperdat)[wt>0]
  mdis <- adist[rownames(adist)[rownames(adist) %in% matched.names],
          colnames(adist)[colnames(adist) %in% matched.names]]
  if(is.matrix(mdis)){
    minT <- apply(mdis,MAR=1,min)
    minC <- apply(mdis,MAR=2,min) 
  } else { 
    minT <- min(mdis)
    minC <- min(mdis)
  }     
  return(mean(c(minT,minC)))
}
## END FUNCTIONS
####################################################################
