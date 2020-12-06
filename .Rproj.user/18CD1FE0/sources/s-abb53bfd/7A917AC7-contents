# Sim. #1: Selction Bias, Unequal Importance of Variables ---------------------

sessionInfo()
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

setwd("~/PSM_research/simple_dist_matching_paradox/sim1")

#Load Functions
source("aux_scripts/sim_functions.R")
source("aux_scripts/matchingFunctions_20oct2014.R")

# Run the simulations
n <- 100 # Number of sims
for(z in 1:n){
  print(paste("Beginning to make Data Set #", z, sep=""))
  
  #set seed for replication
  set.seed(z)
  
  #generate data
  a       <- c(runif(n=100, min=0, max=5), runif(n=100, min=1, max=6))
  b       <- c(runif(n=100, min=0, max=5), runif(n=100, min=1, max=6))
  
  mydat     <- data.frame(a, b)
  mydat$a2  <- a^2
  mydat$b2  <- b^2
  mydat$ab  <- a*b
  mydat$a3  <- a^3
  mydat$b3  <- b^3
  mydat$a2b <- a*a*b
  mydat$ab2 <- a*b*b
  
  # Simulate selection bias
  #   to avoid breaking greedymatch
  #   keep nC=nT, but let treat=1 
  #   be more prevalent among those
  #   with high values of a
  tmp  <- mydat[mydat$a > (mean(mydat$a)),]
  tmp1 <- mydat[mydat$a <= (mean(mydat$a)),]
  
  tmp$treat  <- sample(c(0,1,1), nrow(tmp), replace=T)
  ntreat     <- 100-sum(tmp$treat) #num. of remaining treatment cases to create
  tmp1$treat <- sample(c(rep(1,ntreat), rep(0,nrow(tmp1)-ntreat)))
  
  mydat <- rbind(tmp, tmp1)
  
  # The selection bias variable must be correlated with Y
  #   To simulate uneqal feature importance, leave b
  #     unrelated to Y
  Y <- rnorm(n = 200, mean = mydat$a + mydat$treat, sd=1) # add to mydat after PS logit
  
  # Simulate a human selecting a model via stepwise
  full.model <- glm("treat ~ a + b", data = mydat, family = binomial)
  summary(full.model)
  
  library(MASS)
  library(magrittr)
  step.model <- step(full.model, direction = "both", test = "F")
  summary(step.model)
  
  # add Y, plotting features to data set
  plotlab <- c(rep("C", 100), rep("T", 100))
  plotcol <- c(rep("red", 100), rep("blue", 100))
  
  mydat$plotlab <- plotlab
  mydat$plotcol <- plotcol
  mydat$Y       <- Y
  
  mydat$ps  <- fitted(glm(step.model$formula, family=binomial(link="logit"), data=mydat))
  
  #make distance matrices
  md   <- mdistfun(formula("treat~a+b"), mydat)
  ed   <- edistfun(mydat)
  myps <- mydat$ps
  names(myps) <- rownames(mydat)
  pd          <- psdistfun(mydat, myps)
  
  #rank matches
  mrank <- greedymatch(md) # I had to add [1] to break ties
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
    pestmat[1,i] <- coef(lm(myformula, data=pdat))[2]
    
    #run the other 511 models using all combinations of covariates
    for(j in 2:512){
      mcvrs        <- mycovars[modmatrix[j,]==1]
      myformula    <- as.formula(paste("Y ~ treat + ", paste(mcvrs, collapse="+")))
      mestmat[j,i] <- coef(lm(myformula, data=mdat))[2]
      eestmat[j,i] <- coef(lm(myformula, data=edat))[2]
      pestmat[j,i] <- coef(lm(myformula, data=pdat))[2]
    }
  }
  saveRDS(mestmat, paste("saved_objects/mdm", z, ".rds", sep=""))
  saveRDS(eestmat, paste("saved_objects/edm", z, ".rds", sep=""))
  saveRDS(pestmat, paste("saved_objects/psm", z, ".rds", sep=""))
}

# Quick viz
plot(rev(apply(pestmat,2,var))) # use rev() because the loop above worked "backwards"
plot(rev(apply(mestmat,2,var)))
plot(rev(apply(pestmat,2,mean)))
plot(rev(apply(mestmat,2,mean)))

plot(rev(abs(1-apply(pestmat,2,mean))))
plot(rev(abs(1-apply(mestmat,2,mean))))
plot(rev(abs(1-apply(pestmat,2,mean))))
lines(rev(abs(1-apply(mestmat,2,mean))),col="blue")

# Max Estimate
plot(rev(apply(mestmat,2,max)),type="l", ylim=c(0,3))
lines(rev(apply(pestmat,2,max)),lty=2)

# Min Estimate
plot(rev(apply(mestmat,2,min)),type="l", ylim=c(0,3))
lines(rev(apply(pestmat,2,min)),lty=2)

# Mean Estimate
plot(rev(apply(mestmat,2,mean)),type="l", ylim=c(0,3))
lines(rev(apply(pestmat,2,mean)),lty=2)


# Data Visualizations -----------------------------------------------------

n <- 100 ## the number of sims

#### Between-model Variance (Mean Across 100 Data Sets)
allMDM <- allPSM <- matrix(NA, nrow=n, ncol=85)

for(j in 1:n){
  mymdm <- readRDS(paste("saved_objects/mdm", j, ".rds", sep=""))
  mypsm <- readRDS(paste("saved_objects/psm", j, ".rds", sep=""))
  
  allMDM[j,] <- rev(apply(mymdm, 2, var))
  allPSM[j,] <- rev(apply(mypsm, 2, var))
  
}

var.m <- apply(allMDM, 2, mean)
var.p <- apply(allPSM, 2, mean)

pdf("../figs/sim1_between_model_variance_selection_bias_sw.pdf")
plot(var.m, type="l", lwd=2, col="black", ylab="Variance", 
     xlab="Number of Units Pruned", xaxt="n", ylim=c(0,0.025))
axis(side=1,
     at=c(0, 20, 40, 60, 80),
     labels=c(0, 40, 80, 120, 160))
lines(var.p, type="l", lwd=2, col=rgb(red=1, green=0, blue=0, alpha=.7))
text(x=25, y=0.0007, labels="MDP-OLS")
text(x=40, y=0.004, labels="PSP-OLS")
dev.off()


#### bias plots (Estimates)

## Max Estimate
allMDM <- allPSM <- matrix(NA, nrow=n, ncol=85)

for(j in 1:n){
  mymdm <- readRDS(paste("saved_objects/mdm", j, ".rds", sep=""))
  mypsm <- readRDS(paste("saved_objects/psm", j, ".rds", sep=""))
  
  allMDM[j,] <- rev(apply(mymdm, 2, max))
  allPSM[j,] <- rev(apply(mypsm, 2, max))	
}

est.m <- apply(allMDM, 2, mean)
est.p <- apply(allPSM, 2, mean)

pdf("../figs/sim1_max_est_selection_bias_sw.pdf")
plot(est.m, type="l", lwd=2, col="black", 
     ylab="Maximum Coefficient across 512 Specifications", 
     xlab="Number of Units Pruned", xaxt="n",
     ylim=c(0,3))
axis(side=1,
     at=c(0, 20, 40, 60, 80),
     labels=c(0, 40, 80, 120, 160))
lines(est.p, type="l", lwd=2, col=rgb(red=1, green=0, blue=0, alpha=.7))
text(x=70, y=0.7, labels="MDP-OLS")
text(x=70, y=1.3, labels="PSP-OLS")
text(x=10, y=1.05, labels="True effect = 1")
abline(h=1,lty=2)
dev.off()

## Mean Estimator
allMDM <- allPSM <- matrix(NA, nrow=n, ncol=85)

for(j in 1:n){
  mymdm <- readRDS(paste("saved_objects/mdm", j, ".rds", sep=""))
  mypsm <- readRDS(paste("saved_objects/psm", j, ".rds", sep=""))
  
  allMDM[j,] <- rev(apply(mymdm, 2, mean))
  allPSM[j,] <- rev(apply(mypsm, 2, mean))	
}

est.m <- apply(allMDM, 2, mean)
est.p <- apply(allPSM, 2, mean)

pdf("../figs/sim1_mean_est_selection_bias_sw.pdf")
plot(est.m, type="l", lwd=2, col="black", 
     ylab="Mean Coefficient across 512 Specifications", 
     xlab="Number of Units Pruned", xaxt="n",
     ylim=c(0.5,1.5))
axis(side=1,
     at=c(0, 20, 40, 60, 80),
     labels=c(0, 40, 80, 120, 160))
lines(est.p, type="l", lwd=2, col=rgb(red=1, green=0, blue=0, alpha=.7))
text(x=70, y=0.9, labels="MDP-OLS")
text(x=70, y=1.1, labels="PSP-OLS")
text(x=9, y=0.95, labels="True effect = 1")
abline(h=1,lty=2)
dev.off()

## MSE
allMDM <- allPSM <- matrix(NA, nrow=n, ncol=85)

for(j in 1:n){
  mymdm <- readRDS(paste("saved_objects/mdm", j, ".rds", sep=""))
  mypsm <- readRDS(paste("saved_objects/psm", j, ".rds", sep=""))
  
  allMDM[j,] <- rev(apply(mymdm, 2, function(x){mean((x-1)^2)}))
  allPSM[j,] <- rev(apply(mypsm, 2, function(x){mean((x-1)^2)}))
}

mse.m <- apply(allMDM, 2, mean)
mse.p <- apply(allPSM, 2, mean)

pdf("../figs/sim1_mse_selection_bias_sw.pdf")
plot(mse.m, type="l", lwd=2, col="black", 
     ylab="Mean Square Error across 512 Specifications", 
     xlab="Number of Units Pruned", xaxt="n",
     ylim=c(0,0.2))
axis(side=1,
     at=c(0, 20, 40, 60, 80),
     labels=c(0, 40, 80, 120, 160))
lines(mse.p, type="l", lwd=2, col=rgb(red=1, green=0, blue=0, alpha=.7))
text(x=65, y=0.04, labels="MDP-OLS")
text(x=70, y=0.08, labels="PSP-OLS")
dev.off()

