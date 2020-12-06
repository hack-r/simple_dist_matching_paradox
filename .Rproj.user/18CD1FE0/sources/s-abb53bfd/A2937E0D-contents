# Same as sim4, but now there are 2 subsegments with different treatment effects

# Sim. #4: ATT for PSM and MDM + New Data + Alternate Y--------------------

sessionInfo()
# R version 3.6.3 (2020-02-29)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 18.04.4 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1
# 
# Random number generation:
#   RNG:     Mersenne-Twister 
# Normal:  Inversion 
# Sample:  Rounding 
# 
# locale:
#   [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
# [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
# [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
# [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] magrittr_1.5   Matching_4.9-7 MASS_7.3-51.5 
# 
# loaded via a namespace (and not attached):
#   [1] compiler_3.6.3 tools_3.6.3    packrat_0.5.0 

setwd("~/PSM_research/simple_dist_matching_paradox/sim4")

#Load Functions
source("aux_scripts/sim_functions.R")
source("aux_scripts/matchingFunctions_20oct2014.R")
require(Matching)

# Run the simulations
n <- 100 # Number of sims
for(z in 1:n){
  print(paste("Beginning to make Data Set #", z, sep=""))
  
  # set seed for replication
  set.seed(z)
  
  # generate data
  a       <- rnorm(n=200, mean = 10, sd = 1)
  b       <- rnorm(n=200, mean = 10, sd = 10)
  c       <- rnorm(n=200, mean = 0, sd = 3)
  d       <- c*.48 + rnorm(n=200, mean = 0, sd = 10)
  e       <- (b*.5)^2 +  rnorm(n=200, mean = 0, sd = 10)
  
  mydat     <- data.frame(a, b, c, d, e)
  
  # Simulate selection bias
  #   to avoid breaking greedymatch
  #   keep nC=nT, but let treat=1 
  #   be more prevalent among those
  #   with high values of a
  tmp  <- mydat[mydat$a > (mean(mydat$a)),]
  tmp1 <- mydat[mydat$a <= (mean(mydat$a)),]
  
  tmp$treat  <- sample(c(0,1,1,1), nrow(tmp), replace=T)
  ntreat     <- 100-sum(tmp$treat) #num. of remaining treatment cases to create
  tmp1$treat <- sample(c(rep(1,ntreat), rep(0,nrow(tmp1)-ntreat)))
  
  mydat <- rbind(tmp, tmp1)
  
  # The selection bias variable must be correlated with Y
  #   To simulate uneqal feature importance, leave b
  #     unrelated to Y
  Y <- rnorm(n = 200, mean = ifelse(mydat$a > mean(mydat$a), 5, 10) + mydat$treat, sd=1) # add to mydat after PS logit
  
  # Simulate a human selecting a model
  full.model <- glm("treat ~ a + b + c + d + e", data = mydat, family = binomial)
  summary(full.model)
  
  form <- capture.output(cat("treat ~ 1",intersect(names(coef(summary(full.model))[,4][coef(summary(full.model))[,4]<0.05]),c("a", "b", "c", "d", "e")),sep = "+"))
  if(substr(form,nchar(form),nchar(form))=="+") form <- "treat~1"
  form <- formula(form)
  
  final.model <- glm(form, data = mydat, family = binomial)
  summary(final.model)
  
  # add Y, plotting features to data set
  plotlab <- c(rep("C", 100), rep("T", 100))
  plotcol <- c(rep("red", 100), rep("blue", 100))
  
  mydat$plotlab <- plotlab
  mydat$plotcol <- plotcol
  mydat$Y       <- Y
  
  mydat$ps  <- fitted(glm(final.model$formula, family=binomial(link="logit"), data=mydat))
  
  #make distance matrices
  md          <- mdistfun(formula("treat~a+b+c+d+e"), mydat) # no statistical model so you don't know to exclude b/c
  myps        <- mydat$ps
  names(myps) <- rownames(mydat)
  pd          <- psdistfun(mydat, myps)
  
  #rank matches
  mrank <- greedymatch(md)
  prank <- greedymatch(pd)
  
  #holder matrices
  mestmat <- pestmat <- matrix(NA, nrow=512, ncol=85)
  
  #run the caliper
  for(i in 1:85){
    print(paste("caliper: ", i, sep=""))
    pairs  <- i+15
    mnames <- c(mrank[1:pairs,1], mrank[1:pairs,2])
    pnames <- c(prank[1:pairs,1], prank[1:pairs,2])
    
    #retain the best i+15 matches 
    mdat <- mydat[rownames(mydat) %in% mnames,]
    pdat <- mydat[rownames(mydat) %in% pnames,]
    
    #ATT
    effect       <- Match(Y        = mdat$Y, 
                          Tr       = mdat$treat,
                          X        = rep(0,nrow(mdat)),
                          estimand = "ATT",
                          M        = 1, 
                          ties     = F,
                          replace  = F
    )
    mestmat[1,i] <- effect$est
    
    effect       <- Match(Y        = pdat$Y, 
                          Tr       = pdat$treat,
                          X        = pdat$ps,
                          estimand = "ATT",
                          M        = 1, 
                          ties     = F,
                          replace  = F
    )
    pestmat[1,i] <- effect$est
    
    #since the rest of the code assumes 512 rows we'll temporarily
    #patch it like this. this can come out later.
    for(j in 2:512){
      mestmat[j,i] <- mestmat[1,i]
      pestmat[j,i] <- pestmat[1,i]
    }
  }
  saveRDS(mestmat, paste("saved_objects/sim4_att_mdm_newdata_alty", z, ".rds", sep=""))
  saveRDS(pestmat, paste("saved_objects/sim4_att_psm_newdata_alty", z, ".rds", sep=""))
}

# Quick viz

# Mean Estimate
plot(rev(apply(mestmat,2,mean)),type="l", ylim=c(-2,3))
lines(rev(apply(pestmat,2,mean)),lty=1,col="red")
abline(1,0,lty=2)

# MAE
plot(rev(abs(1-apply(pestmat,2,mean))),col="red", ylim=c(0,3.5)) 
lines(rev(abs(1-apply(mestmat,2,mean))),col="black")
abline(0,0)



# Data Visualizations -----------------------------------------------------
n <- 100 # number of sims

#### bias plots (Estimates)
allMDM <- allPSM <- matrix(NA, nrow=n, ncol=85)

for(j in 1:n){
  mymdm <- readRDS(paste("saved_objects/sim4_att_mdm_newdata_alty", j, ".rds", sep=""))
  mypsm <- readRDS(paste("saved_objects/sim4_att_psm_newdata_alty", j, ".rds", sep=""))
  
  allPSM[j,] <- rev(apply(mypsm, 2, mean))	
  allMDM[j,] <- rev(apply(mymdm, 2, mean))
}

est.m <- apply(allMDM, 2, mean)
est.p <- apply(allPSM, 2, mean)

pdf("../figs/sim4_true_psm_mdm_att_newdata_alty.pdf")
plot(est.p, type="l", lwd=2, col="red", 
     ylab="PSM vs MDM ATT across Pruning Levels", 
     xlab="Number of Units Pruned", xaxt="n") #,ylim=c(0.9,2.3)
axis(side=1,
     at=c(0, 20, 40, 60, 80),
     labels=c(0, 40, 80, 120, 160))
lines(est.m, type="l", lwd=2, col="black", lty=1)
text(x=20, y=1.5, labels="MDM")
text(x=20, y=1.2, labels="PSM")
text(x=10, y=1.1, labels="True effect = 1")
abline(h=1,lty=2)
dev.off()

# MAE
for(j in 1:n){
  mymdm <- readRDS(paste("saved_objects/sim4_att_mdm_newdata_alty", j, ".rds", sep=""))
  mypsm <- readRDS(paste("saved_objects/sim4_att_psm_newdata_alty", j, ".rds", sep=""))
  
  allPSM[j,] <- rev(abs(1-apply(mypsm, 2, mean)))	
  allMDM[j,] <- rev(abs(1-apply(mymdm, 2, mean)))
}

est.m <- apply(allMDM, 2, mean)
est.p <- apply(allPSM, 2, mean)

pdf("../figs/sim4_true_psm_mdm_mae_newdata_alty.pdf")
plot(est.p, type="l", lwd=2, col="red", 
     ylab="PSM vs MDM MAE across Pruning Levels", 
     xlab="Number of Units Pruned", xaxt="n") #,ylim=c(0.9,2.3)
axis(side=1,
     at=c(0, 20, 40, 60, 80),
     labels=c(0, 40, 80, 120, 160))
lines(est.m, type="l", lwd=2, col="black", lty=1)
text(x=60, y=0.5, labels="PSM")
text(x=60, y=1, labels="MDM")
dev.off()
