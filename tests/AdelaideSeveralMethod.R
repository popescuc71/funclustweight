library(flexmix)
library(fda.usc)
library(fds)
library(mclust)
library(stringr)
library(MASS)

library(funclustweight)
library(funHDDC)
library(fda)
library(mclust)
de_sun <- fds::sundaydemand

de_tue <- fds::tuesdaydemand

#****** setup labels for days ******#
# when checking dimensions using the $x nrow is timepoints
cls <- append(rep("sunday", ncol(de_sun$y)),
              rep("tuesday", ncol(de_tue$y)))

#****** combine and format data ******#
demand_raw <- cbind(de_sun$y, de_tue$y)
draw_x <- t(demand_raw[1:24,])
draw_y <- t(demand_raw[25:48,])

demand <- fitAdelaideFD()
plotAdelaideFD(demand)
demand_coef <- cbind(t(demand$fdx$coefs), t(demand$fdy$coefs))
#***********************USING funweightclust and funHDDC***************************#
models <- c("AKJBKQKDK","AKJBQKDK", "AKBKQKDK", "ABKQKDK", "AKBQKDK", "ABQKDK")
modelsy <- c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "VEE", "EVE",
             "EEV", "VVE", "VEV","EVV","VVV")
#at least nb.rep=20 is recommended, but will take longer; 
set.seed(12345)
res_reg <- funclustweight(demand$fdx, demand$fdy, K=2, model = models,
                          modely = modelsy, init = "kmeans", nb.rep = 1,
                          threshold = 0.1)
table(demand$groupd, res_reg$class)
adjustedRandIndex(demand$groupd, res_reg$class)

#***********************USING funHDDC***************************#
#at least nb.rep=20 is recommended, but will take longer;
#2D-funHDDC
set.seed(12345)
res=funHDDC(list(demand$fdx, demand$fdy), K=2, threshold=0.2, init="kmeans",nb.rep=20, model=models)
table(demand$groupd, res$class)
adjustedRandIndex(demand$groupd, res$class)
#at least nb.rep=20 is recommended, but will take longer;
set.seed(12345)
res=funHDDC(demand$fdfull, K=2, threshold=0.1, init="kmeans",nb.rep=20, model=models)
table(demand$groupd, res$class)
adjustedRandIndex(demand$groupd, res$class)
#******************************** USING OTHER METHODS *******************************#

flexmix_raw <- flexmix(~ draw_x, data = data.frame(draw_y), k = 2,
                       model = list(FLXMRglm(X25 ~ .),
                            FLXMRglm(X26 ~ .),
                            FLXMRglm(X27 ~ .),
                            FLXMRglm(X28 ~ .),
                            FLXMRglm(X29 ~ .),
                            FLXMRglm(X30 ~ .),
                            FLXMRglm(X31 ~ .),
                            FLXMRglm(X32 ~ .),
                            FLXMRglm(X33 ~ .),
                            FLXMRglm(X34 ~ .),
                            FLXMRglm(X35 ~ .),
                            FLXMRglm(X36 ~ .),
                            FLXMRglm(X37 ~ .),
                            FLXMRglm(X38 ~ .),
                            FLXMRglm(X39 ~ .),
                            FLXMRglm(X40 ~ .),
                            FLXMRglm(X41 ~ .),
                            FLXMRglm(X42 ~ .),
                            FLXMRglm(X43 ~ .),
                            FLXMRglm(X44 ~ .),
                            FLXMRglm(X45 ~ .),
                            FLXMRglm(X46 ~ .),
                            FLXMRglm(X47 ~ .),
                            FLXMRglm(X48 ~ .)))
table(cls, flexmix_raw@cluster)
adjustedRandIndex(cls,flexmix_raw@cluster)
# COEFFICIENTS
dcoef_x <- t(as.matrix(demand$fdx$coefs))
dcoef_y <- t(demand$fdy$coefs)
flexmix_coef <- flexmix(~ dcoef_x, data = data.frame(dcoef_y), k = 2,
                       model = list(FLXMRglm(bspl4.1 ~ .),
                                    FLXMRglm(bspl4.2 ~ .),
                                    FLXMRglm(bspl4.3 ~ .),
                                    FLXMRglm(bspl4.4 ~ .),
                                    FLXMRglm(bspl4.5 ~ .),
                                    FLXMRglm(bspl4.6 ~ .)))

table(cls, flexmix_coef@cluster)
adjustedRandIndex(cls,flexmix_coef@cluster)


#####

Mclust_raw <- Mclust(t(demand_raw), G=2)
table(cls, Mclust_raw$classification)
adjustedRandIndex(cls,Mclust_raw$classification)
Mclust_coef <- Mclust(demand_coef, G=2)
table(cls, Mclust_coef$classification)
adjustedRandIndex(cls,Mclust_coef$classification)

#####

kmeans_raw <- kmeans(t(demand_raw), centers = 2)
table(cls, kmeans_raw$cluster)
adjustedRandIndex(cls,kmeans_raw$cluster)
kmeans_coef <- kmeans(demand_coef, centers = 2)
table(cls, kmeans_coef$cluster)
adjustedRandIndex(cls,kmeans_coef$cluster)

