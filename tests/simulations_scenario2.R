#install.packages("fda", repos="https://mirror.rcg.sfu.ca/mirror/CRAN/")
#install.packages("fds", repos="https://mirror.rcg.sfu.ca/mirror/CRAN/")
#install.packages("foreach", repos="https://mirror.rcg.sfu.ca/mirror/CRAN/")
#install.packages("doParallel", repos="https://mirror.rcg.sfu.ca/mirror/CRAN/")
#install.packages("doFuture", repos="https://mirror.rcg.sfu.ca/mirror/CRAN/")
#install.packages("future.batchtools", repos="https://mirror.rcg.sfu.ca/mirror/CRAN/")
#install.packages("doRNG", repos="https://mirror.rcg.sfu.ca/mirror/CRAN/")

#******************************* REQUIRED LIBRARY *****************************#
library(funHDDC)
library(fda)
library(fds)
library(doParallel)
library(doFuture)
library(foreach)
library(mclust)
library(future.batchtools)
library(doRNG)


registerDoFuture() # registers packages and functions that may be used
plan(multicore)
#plan(sequential)
library(funclustweight)
#library(funweightclust)
rand_seeds <- read.csv("random_seeds.csv")$x
models=c("AKJBKQKDK","AKJBQKDK", "AKBKQKDK", "ABKQKDK", "AKBQKDK", "ABQKDK")
modelsy=c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "VEE", "EVE", "EEV", "VVE", "VEV","EVV","VVV")


res_reg=readRDS("input1n.RData")
#ntrials should be 100 but it takes hours
ntrials = 1
#ntrials=100
all_thresh <- c(0.005, 0.01, 0.1, 0.2)
#faster with one threshold value
#all_thresh <- c(0.005)


fullclass <- foreach(i = 30:(ntrials+29), .combine=rbind, .init=NULL, .errorhandling="remove") %dorng% {
    current_class <- data.frame(matrix(ncol = 608,nrow=0))
 #   
    seed_fds <- rand_seeds[i+1]
    set.seed(seed_fds)
    
    # produce functional data

   model_reg <- funclustweight::genFromModelRegFD(res_reg, ncurves = c(300, 300),mu_div = 1.4)
  #plot the simulated data
   #funclustweight::plotRegFD(model_reg)
   
     # save the grouping data
    for(thresh in all_thresh) {
      ### GET RESULTS FROM FUNWEIGHT CLUST
      ##nb.rep should be 20 but it will take long
      res_funweight <- funclustweight(model_reg$xfd, model_reg$yfd,
                                      model = models, modely = modelsy,
                                      threshold = thresh, nb.rep = 1,
                                      itermax = 200, K = 2, init = "kmeans",
                                      verbose = F)
      if(!is.null(res_funweight$class)) {
        ccr <- sum(diag(table(model_reg$groupd, res_funweight$class)))/600
        ccr <- ifelse(ccr < 0.5, 1 - ccr, ccr)
        ari <- adjustedRandIndex(model_reg$groupd, res_funweight$class)
        current_class <- rbind(current_class, c(i, seed_fds, "funweightclust",
                                                thresh, ccr, ari,
                                                res_funweight$model,
                                                res_funweight$modely,
                                                res_funweight$class))
      }
      ### END OF GET RESULTS FROM FUNWEIGHT CLUST

      ### GET RESULTS FROM FUNHDDC LIST
      ##nb.rep should be 20 but it will take long
      res_funHDDC_list <- funHDDC(list(model_reg$xfd, model_reg$yfd),
                                      model = models,
                                      threshold = thresh, nb.rep = 1,
                                      itermax = 200, K = 2, init = "kmeans")
      if(!is.null(res_funHDDC_list$class)) {
      ccr <- sum(diag(table(model_reg$groupd, res_funHDDC_list$class)))/600
      ccr <- ifelse(ccr < 0.5, 1 - ccr, ccr)
      ari <- adjustedRandIndex(model_reg$groupd, res_funHDDC_list$class)
      current_class <- rbind(current_class, c(i, seed_fds, "funHDDC_list",
                                              thresh, ccr, ari,
                                              res_funHDDC_list$model,
                                              "funHDDC",
                                              res_funHDDC_list$class))
      }
      ### END OF GET RESULTS FROM FUNHDDC LIST

      
    }
    print(paste0("Finished Clustering - ", "Test: ", i))
    names(current_class) <- c("trial", "seed", "method", "threshold", "CCR", "ARI", "modelx", "modely", paste0("C", 1:600))
    write.csv(current_class, file = paste0(i, "_simScenario2.csv"))
    current_class
  } # end of for each clust_FUN
names(fullclass) <- c("trial", "seed", "method", "threshold", "CCR", "ARI", "modelx", "modely", paste0("C", 1:600))
write.csv(fullclass, file = "sim_fulltrial_Scenario2.csv")
