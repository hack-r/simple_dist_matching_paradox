
# Sim. #2: Actual PS Matching ---------------------------------------------
sessionInfo()
# R version 3.6.3 (2020-02-29)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 18.04.4 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1
# 
# locale:
#   [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8        LC_COLLATE=C.UTF-8    
# [5] LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8    LC_PAPER=C.UTF-8       LC_NAME=C             
# [9] LC_ADDRESS=C           LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# loaded via a namespace (and not attached):
#   [1] MLmetrics_1.1.1 compiler_3.6.3  tools_3.6.3     packrat_0.5.0 

setwd("~/PSM_research/simple_dist_matching_paradox/sim2")

#Load Functions
source("aux_scripts/sim_functions.R")
source("aux_scripts/matchingFunctions_20oct2014.R")

# King's Data + True PSM -----------------------------------------
require(Matching)
allMDM1 <- allPSM1 <- matrix(NA, nrow=n, ncol=85)

for(z in 1:100){
  # Announce iteration
  print(paste("Beginning to make Data Set #", z, sep=""))
  
  #set seed for replication
  set.seed(z)
  
  #generate data
  a       <- c(runif(n=100, min=0, max=5), runif(n=100, min=1, max=6))
  b       <- c(runif(n=100, min=0, max=5), runif(n=100, min=1, max=6))
  treat   <- c(rep(0, 100), rep(1, 100))
  plotlab <- c(rep("C", 100), rep("T", 100))
  plotcol <- c(rep("red", 100), rep("blue", 100))
  Y       <- rnorm(n=200, mean=a+b + 2*treat, sd=1)
  
  psformula <- as.formula("treat ~ a + b")
  
  #set up dataframe with all variables needed for models
  mydat     <- data.frame(Y, treat, a, b, plotlab, plotcol)
  mydat$a2  <- a^2
  mydat$b2  <- b^2
  mydat$ab  <- a*b
  mydat$a3  <- a^3
  mydat$b3  <- b^3
  mydat$a2b <- a*a*b
  mydat$ab2 <- a*b*b
  mydat$ps  <- fitted(glm(psformula, family=binomial(link="logit"), data=mydat))
  ps.model  <- glm(psformula, family=binomial(link="logit"), data=mydat)
  model.sig <- with(ps.model, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = FALSE))
  
  #make distance matrices
  md          <- mdistfun(psformula, mydat)
  ed          <- edistfun(mydat)
  myps        <- mydat$ps
  names(myps) <- rownames(mydat)
  pd          <- psdistfun(mydat, myps)
  
  #rank matches
  mrank <- greedymatch(md)
  erank <- greedymatch(ed)
  prank <- greedymatch(pd)
  
  #covariates for all models
  mycovars <- c("a", "b", "a2", "b2", "ab", "a3", "b3", "a2b", "ab2")
  
  #make a matrix of all binary numbers 1 <= n <= 512
  modmatrix <- matrix(NA, nrow=512, ncol=9)
  for(i in 1:512) modmatrix[i,] <- as.numeric(intToBits(i-1))[1:9]
  
  #holder matrices
  mestmat <- eestmat <- pestmat <- matrix(NA, nrow=512, ncol=85)
  
  #run the caliper
  for(i in 1:85){
    print(paste("caliper: ", i, sep=""))
    pairs  <- i+15
    mnames <- c(mrank[1:pairs,1], mrank[1:pairs,2])
    enames <- c(erank[1:pairs,1], erank[1:pairs,2])
    pnames <- c(prank[1:pairs,1], prank[1:pairs,2])
    
    #retain the best i+15 matches 
    mdat <- mydat[rownames(mydat) %in% mnames,]
    edat <- mydat[rownames(mydat) %in% enames,]
    pdat <- mydat[rownames(mydat) %in% pnames,]
    
    #the difference in mean estimator
    myformula    <- as.formula("Y ~ treat")
    mestmat[1,i] <- coef(lm(myformula, data=mdat))[2]
    eestmat[1,i] <- coef(lm(myformula, data=edat))[2]
    effect       <- Match(Y        = pdat$Y, 
                          Tr       = pdat$treat,
                          X        = pdat$ps,
                          estimand = "ATT",
                          M        = 1, 
                          ties     = F,
                          replace  = F,
                          exact    = F#,caliper  = 0.00001
    )
    pestmat[1,i] <- effect$est
    
    #run the other 511 models using all combinations of covariates
    for(j in 2:512){
      mcvrs        <- mycovars[modmatrix[j,]==1]
      myformula    <- as.formula(paste("Y ~ treat + ", paste(mcvrs, collapse="+")))
      mestmat[j,i] <- coef(lm(myformula, data=mdat))[2]
      eestmat[j,i] <- coef(lm(myformula, data=edat))[2]
      #pestmat[j,i] <- coef(lm(myformula, data=pdat))[2]
      # PSM Does not use regression to estimate the treatment effect!
      #   That's the whole point...
      pestmat[j,i] <- effect$est
    }
  }
  saveRDS(mestmat, paste("saved_objects/mdm_kingdata", z, ".rds", sep=""))
  saveRDS(eestmat, paste("saved_objects/edm_kingdata", z, ".rds", sep=""))
  saveRDS(pestmat, paste("saved_objects/psm_kingdata", z, ".rds", sep=""))
}

# Quick viz
plot(rev(apply(pestmat,2,var))) # use rev() because the loop above worked "backwards"
plot(rev(apply(mestmat,2,var)))
plot(rev(apply(pestmat,2,mean)))
plot(rev(apply(mestmat,2,mean)))

plot(rev(abs(1-apply(pestmat,2,mean))))
lines(rev(abs(1-apply(mestmat,2,mean))),col="blue")

# Max Estimate
plot(rev(apply(mestmat,2,max)),type="l")
lines(rev(apply(pestmat,2,max)),lty=2)

# Min Estimate
plot(rev(apply(mestmat,2,min)),type="l", ylim=c(0,3))
lines(rev(apply(pestmat,2,min)),lty=2)

# Mean Estimate
plot(rev(apply(mestmat,2,mean)),type="l", ylim=c(0,3))
lines(rev(apply(pestmat,2,mean)),lty=2)


# Data Visualizations -----------------------------------------------------
n <- 100 # number of sims

#### Between-model Variance (Mean Across 100 Data Sets)
allMDM <- allPSM <- matrix(NA, nrow=n, ncol=85)

for(j in 1:n){
  mymdm <- readRDS(paste("saved_objects/mdm_kingdata", j, ".rds", sep=""))
  mypsm <- readRDS(paste("saved_objects/psm_kingdata", j, ".rds", sep=""))
  
  allMDM[j,] <- rev(apply(mymdm, 2, var))
  allPSM[j,] <- rev(apply(mypsm, 2, var))
  
}

var.m <- apply(allMDM, 2, mean)
var.p <- apply(allPSM, 2, mean)

pdf("../figs/sim2_true_psm_between_model_variance_kingdata.pdf")
plot(var.m, type="l", lwd=2, col="black", ylab="Variance", 
     xlab="Number of Units Pruned", xaxt="n", ylim=c(0,0.05))
axis(side=1,
     at=c(0, 20, 40, 60, 80),
     labels=c(0, 40, 80, 120, 160))
lines(var.p, type="l", lwd=2, col=rgb(red=1, green=0, blue=0, alpha=.7))
text(x=25, y=0.02, labels="MDP-OLS")
text(x=20, y=0.003, labels="PSM")
dev.off()

#### bias plots (Estimates)
# Since PSM will be the same (max,mean,min) - 1 plot for all 3 
#   makes sense this time 
allMDM.max <- allMDM.min <- allMDM.mean <- allPSM <- matrix(NA, nrow=n, ncol=85)

for(j in 1:n){
  mymdm <- readRDS(paste("saved_objects/mdm_kingdata", j, ".rds", sep=""))
  mypsm <- readRDS(paste("saved_objects/psm_kingdata", j, ".rds", sep=""))
  
  allPSM[j,] <- rev(apply(mypsm, 2, mean))	

  allMDM.max[j,]  <- rev(apply(mymdm, 2, max))
  allMDM.mean[j,] <- rev(apply(mymdm, 2, mean))
  allMDM.min[j,]  <- rev(apply(mymdm, 2, min))
}

est.m.min  <- apply(allMDM.min, 2, mean)
est.m.mean <- apply(allMDM.mean, 2, mean)
est.m.max  <- apply(allMDM.max, 2, mean)
est.p      <- apply(allPSM, 2, mean)

pdf("../figs/sim2_true_psm_allest_kingdata.pdf")
plot(est.p, type="l", lwd=2, col="red", 
     ylab="PSM vs Max, Mean, and Min MDP-OLS Coef.s", 
     xlab="Number of Units Pruned", xaxt="n",ylim=c(1.7,4))
axis(side=1,
     at=c(0, 20, 40, 60, 80),
     labels=c(0, 40, 80, 120, 160))
lines(est.m.mean, type="l", lwd=2, col="black", lty=3)
lines(est.m.min, type="l", lwd=2, col="black", lty=2)
lines(est.m.max, type="l", lwd=2, col="black", lty=1)
text(x=20, y=2.1, labels="True effect = 2")
legend(x="topright",  c("PSM", "MDP-OLS Mean", "MDP-OLS Max", "MDP-OLS Min"), lty=c(1,3,1,2),
       col=c("red","black","black","black"),cex=0.5,inset = c(0, 0),y.intersp=0.25)
abline(h=2,lty=2)
dev.off()

### Max Estimate
allMDM <- allPSM <- matrix(NA, nrow=n, ncol=85)

for(j in 1:n){
  mymdm <- readRDS(paste("saved_objects/mdm_kingdata", j, ".rds", sep=""))
  mypsm <- readRDS(paste("saved_objects/psm_kingdata", j, ".rds", sep=""))

  allMDM[j,] <- rev(apply(mymdm, 2, max))
  allPSM[j,] <- rev(apply(mypsm, 2, max))
}

est.m <- apply(allMDM, 2, mean)
est.p <- apply(allPSM, 2, mean)

pdf("../figs/sim2_true_psm_max_est_kingdata.pdf")
plot(est.m, type="l", lwd=2, col="black",
     ylab="Maximum Coefficient across 512 Specifications",
     xlab="Number of Units Pruned", xaxt="n")
axis(side=1,
     at=c(0, 20, 40, 60, 80),
     labels=c(0, 40, 80, 120, 160))
lines(est.p, type="l", lwd=2, col=rgb(red=1, green=0, blue=0, alpha=.7))
text(x=42, y=2.3, labels="MDP-OLS")
text(x=30, y=2.1, labels="PSM")
text(x=10, y=2.1, labels="True effect = 2")
abline(h=2,lty=2)
dev.off()

## MSE - PSM (invariant) vs MDP-OLS Max Estimator
# Follows Rich Nielsen's approach
allMDM <- allPSM <- matrix(NA, nrow=n, ncol=85)

for(j in 1:n){
  mymdm <- readRDS(paste("saved_objects/mdm_kingdata", j, ".rds", sep=""))
  mypsm <- readRDS(paste("saved_objects/psm_kingdata", j, ".rds", sep=""))
  
  allMDM[j,] <- rev(apply(mymdm, 2, function(x){max(x)}))
  allPSM[j,] <- rev(apply(mypsm, 2, function(x){max(x)}))
}

mse.m <- apply(allMDM, 2, function(x){mean((x-2)^2)})
mse.p <- apply(allPSM, 2, function(x){mean((x-2)^2)})

pdf("../figs/sim2_true_psm_mse_kingdata.pdf")
plot(mse.m, type="l", lwd=2, col="black", 
     ylab="MSE: PSM vs. MDP-OLS Max Estimator", 
     xlab="Number of Units Pruned", xaxt="n") #,ylim=c(0,3)
axis(side=1,
     at=c(0, 20, 40, 60, 80),
     labels=c(0, 40, 80, 120, 160))
lines(mse.p, type="l", lwd=2, col=rgb(red=1, green=0, blue=0, alpha=.7))
text(x=24, y=1, labels="MDP-OLS")
text(x=20, y=0.3, labels="PSM")
dev.off()
