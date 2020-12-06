#a function to generate a dataset
makedat <- function(){
	# First generate a dataset as follows
	# 1) Generate 50 unbalanced control units
	# that are distributed uniform on X1 in [-8, -4] and X2 in [-8, 2]
	# 2) Generate 25 treated and control units that are "completely randomized"
	# where both covariates are uniform on X1 in [-2, 2] and X2 in [-8, -4]
	# 3) Generate 25 treated and control units that are "pair randomized"
	# where the control units are uniform on X1 in [-2, 2] and X2 in [-2, 2] and the
	# control units are jittered from the treated units
	# Second Generate Y1...Y20 where Yi is an ith order polynomial function of X1, X2
	# plus a constant treatment effect theta=5
	

	X1 <- X2 <- treat <- rep(NA, 150)
	dat <- data.frame(X1, X2, treat)
	
	#unbalanced
	x1 <- runif(n=50, min=-8, max=-4)
	x2 <- runif(n=50, min=-8, max=2)
	mytreat <- rep(0, 50)
	
	dat[1:50, 1:3] <- cbind(x1, x2, mytreat)
	
	#fully randomized
	x1 <- runif(n=50, min=-2, max=2)
	x2 <- runif(n=50, min=-8, max=-4)
	mytreat <- c(rep(0, 25), rep(1, 25))

	dat[51:100, 1:3] <- cbind(x1, x2, mytreat)

	#pair randomized
	x1c <- runif(n=25, min=-2, max=2)
	x2c <- runif(n=25, min=-2, max=2)
	x1t <- jitter(x1c)
	x2t <- jitter(x2c)
	x1 <- c(x1c, x1t)
	x2 <- c(x2c, x2t)
	
	dat[101:150, 1:3] <- cbind(x1, x2, mytreat)
	
	theta <- 5
	syst <- theta * dat$treat
	
	for (i in 1:20){
		syst <- syst + (1/i^10)*dat$X1^i + (1/i^10)*dat$X2^i
		dat[,paste("Y", i, sep="")] <- rnorm(n=150, mean=syst, sd=5)
	}

	return(dat)
}

#takes a data set with propensity scores and creates a matrix that gives the propensity
#distance from every treated unit to every control unit
psdistfun <- function(mydat, myps){
	
	distmat <- matrix(NA, nrow=sum(mydat$treat), ncol=nrow(mydat)-sum(mydat$treat))
	rownames(distmat) <- rownames(mydat)[mydat$treat==1]
	colnames(distmat) <- rownames(mydat)[mydat$treat==0]	
	for(i in 1:nrow(distmat)){
		tu <- rownames(distmat)[i]
		for(j in 1:ncol(distmat)){
			cu <- colnames(distmat)[j]
			distmat[i,j] <- abs(myps[names(myps)==tu] - myps[names(myps)==cu])			
		}
	}
	return(distmat)
}

#Returns the euclidean distance between two vectors
mydist <- function(a, b){
	return(sqrt(sum((a-b)^2)))
}


edistfun <- function(mydat){
	mytreat <- mydat[mydat$treat==1,]
	mycontrol <- mydat[mydat$treat==0,]
	distmat <- matrix(NA, nrow=nrow(mytreat), ncol=nrow(mycontrol))
	rownames(distmat) <- rownames(mytreat)
	colnames(distmat) <- rownames(mycontrol)
	for(i in 1:nrow(mytreat)){
		a <- mytreat[i, c("a", "b")]
		for(j in 1:nrow(mycontrol)){
			b <- mycontrol[j, c("a", "b")]
			distmat[i, j] <- mydist(a, b)
		}
	}
	return(distmat)
}

# takes a distance matrix produces by psdist (or some other distance matrix function)
# and returns a greedy match. first finds the best match, then the best available 
# match and so on until all pairs have been found. If Nt != Nc, will find 
# pairs equal to min(Nt, Nc)

greedymatch <- function(distmat){
  controls <- colnames(distmat)
  treated  <- rownames(distmat)
  N        <- min(dim(distmat))
  maxiter  <- N-1
  
  bestpairs <- data.frame(matrix(NA, nrow=N, ncol=2))
  colnames(bestpairs) <- c("Treated", "Control")
  for(i in 1:maxiter){
    minElem <- which(distmat==min(distmat), arr.ind=T)
    mt <- minElem[,"row"]
    mc <- minElem[,"col"]
    bestpairs[i,1] <- treated[mt][1]  # added [1] to break ties
    bestpairs[i,2] <- controls[mc][1] # added [1] to break ties
    controls       <- controls[-mc[1]]# added [1] to break ties
    treated        <- treated[-mt[1]] # added [1] to break ties
    distmat        <- distmat[-mt[1], -mc[1]] # added [1] to break ties
  }
  bestpairs[N,1] <- treated
  bestpairs[N,2] <- controls
  return(bestpairs)
}

#returns the difference between the max and min of a vector
myrange <- function(x){
	return(abs(max(x) - min(x)))
}
