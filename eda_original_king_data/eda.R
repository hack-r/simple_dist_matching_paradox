
# EDA of King and Nielsen's Data ------------------------------------------

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
#   [1] xgboost_1.2.0.1 Matrix_1.2-18   optmatch_0.9-13 survival_3.1-11 Matching_4.9-7 
# [6] MASS_7.3-51.5   MatchIt_3.0.2  
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.4.6      lattice_0.20-40   packrat_0.5.0     digest_0.6.25    
# [5] grid_3.6.3        xtable_1.8-4      magrittr_1.5      stringi_1.4.6    
# [9] svd_0.5           data.table_1.12.8 SparseM_1.78      RItools_0.1-17   
# [13] splines_3.6.3     tools_3.6.3       abind_1.4-5       compiler_3.6.3 

# Working Directory:
#     Set to the replication directory on your computer
setwd("~/PSM_research/simple_dist_matching_paradox/eda_original_king_data/")

#Load Functions
source("aux_scripts/sim_functions.R")
source("aux_scripts/matchingFunctions_20oct2014.R")


# MSE ---------------------------------------------------------------------
n <- 100 # Number of sims
allMDM1 <- allPSM1 <- matrix(NA, nrow=n, ncol=85)

for(j in 1:n){
  mymdm <- readRDS(paste("saved_objects/mdm", j, ".rds", sep=""))
  mypsm <- readRDS(paste("saved_objects/psm", j, ".rds", sep=""))
  
  allMDM1[j,] <- rev(apply(mymdm, 2, mean))
  allPSM1[j,] <- rev(apply(mypsm, 2, mean))
}

mse.m <- apply(allMDM1, 2, function(x){mean((x-2)^2)})
mse.p <- apply(allPSM1, 2, function(x){mean((x-2)^2)})

pdf("../figs/eda_king_nielsen_MSE.pdf")
plot(mse.m, type="n", lwd=2, col="black", ylab="MSE", 
     xlab="Number of Units Pruned", xaxt="n")
axis(side=1,
     at=c(0, 20, 40, 60, 80),
     labels=c(0, 40, 80, 120, 160))
lines(mse.p, type="l", lwd=2, col=rgb(red=1, green=0, blue=0, alpha=.7))
lines(mse.m, type="l", lwd=2, col="black")
legend("topleft",c("MDP-OLS","PSP-OLS"),col=c("black","red"),pch=c(1,2))
dev.off()
