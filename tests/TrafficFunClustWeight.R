#******************************* LIBRARIES ************************************#
library(dplyr)
library(fda)
library(funclustweight)

#********************************** DATA **************************************#


fitTrafficFD <- function(splitx = 48, nbasis = 12) {
  # open the data
  weat_surv <- read.csv("survey_weather_joined_2017-2018.csv")

  # format into functional data
  # consider potentially turning the data into a smaller basis
  basis_x <- create.bspline.basis(rangeval = c(0, splitx/4), nbasis = nbasis)
  basis_y <- create.bspline.basis(rangeval = c((splitx+1)/4, 24), nbasis = nbasis)
  surv2 <- weat_surv %>% dplyr::select(contains("Bin.2"))
  surv2_fd_x <- smooth.basis(argvals = seq(0, splitx/4, length.out = splitx),
                           y = t(surv2[,1:splitx]),
                           fdParobj = basis_x)$fd
  surv2_fd_y <- smooth.basis(argvals = seq((splitx+1)/4, 24, length.out = ncol(surv2)-splitx),
                             y = t(surv2[,(splitx+1):96]),
                             fdParobj = basis_y)$fd

  surv4 <- weat_surv %>% dplyr::select(contains("Bin.4"))
  surv4_fd_x <- smooth.basis(argvals = seq(0, splitx/4, length.out = splitx),
                           y = t(surv4[,1:splitx]),
                           fdParobj = basis_x)$fd
  surv4_fd_y <- smooth.basis(argvals = seq((splitx+1)/4, 24, length.out = ncol(surv4)-splitx),
                             y = t(surv4[,(splitx+1):96]),
                             fdParobj = basis_y)$fd

  # sample 1000 curves from the data
  set.seed(54361)
  curv_samp <- sample(1:nrow(weat_surv), 1000, replace = F)
  full_fd <- list(b2x = surv2_fd_x[curv_samp,1:splitx],
                  b4x = surv4_fd_x[curv_samp,1:splitx],
                  b2y = surv2_fd_y[curv_samp,(splitx+1):96],
                  b4y = surv4_fd_y[curv_samp,(splitx+1):96],
                  fast = weat_surv$Speed.Limit[curv_samp])
  return(full_fd)
}

plotTrafficFD <- function(tfd, class = NULL, class_ind = F, pdf_name = "traffic_groups.pdf") {

  if(!is.null(class)) {
    if(!class_ind) {
      par(mfrow = c(2,2))
      plot(tfd$b2x, col = class)
      plot(tfd$b2y, col = class)
      plot(tfd$b4x, col = class)
      plot(tfd$b4y, col = class)
    } else {
      pdf(file = pdf_name)
      par(mfrow = c(2,4))
      for(i in 1:length(unique(class))){
        plot.fd(tfd$b2x[class == i], col = i)
        plot.fd(tfd$b2y[class == i], col = i)
        plot.fd(tfd$b4x[class == i], col = i)
        plot.fd(tfd$b4y[class == i], col = i)
      }
      dev.off()
    }
  } else {
    par(mfrow = c(2,2))
    plot(tfd$b2x)
    plot(tfd$b2y)
    plot(tfd$b4x)
    plot(tfd$b4y)
  }
}

# IMPORTANT #
# for less time processing call with reduced hyper-parameters


models <- c("AKJBKQKDK","AKJBQKDK", "AKBKQKDK", "ABKQKDK", "AKBQKDK", "ABQKDK")
modelsy <- c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "VEE", "EVE",
             "EEV", "VVE", "VEV","EVV","VVV")
#models <- c("AKJBKQKDK")
#modelsy <- c("EII")
# fitting settings
cutting_points <- c(24, 48, 72)
basis_nums <- c(4, 6, 12, 18)
# clustering settings
inits <- c("kmeans", "random")
thresholds <- c(0.005, 0.01, 0.1, 0.2)

for(cutx in cutting_points) for(nbasis in basis_nums) {
    traffic_fd <- fitTrafficFD(splitx = cutx, nbasis = nbasis)
    temp_res <- matrix(ncol = 8+ncol(traffic_fd$b2x$coefs), nrow = 0)
    for(init in inits) for(thresh in thresholds) {
      reg_res <- funclustweight(list(traffic_fd$b2x, traffic_fd$b4x),
                                list(traffic_fd$b2y, traffic_fd$b4y),
                                threshold = thresh,
                                init = init,
                                nb.rep = 50, mc.cores = 25,
                                model = models, modely = modelsy)
      temp_res <- rbind(temp_res,
                        c(cutx, nbasis, thresh, init,
                          reg_res$K, reg_res$model, reg_res$modely, reg_res$BIC,
                          reg_res$class))
      plotTrafficFD(traffic_fd, class = reg_res$class, class_ind = T,
                    pdf_name = paste0("traffic_", cutx, "-", nbasis, "_",
                                      init, "_", thresh, "_groups.pdf"))
    }
  names(temp_res) <- c("cutx", "nbasis", "threshold", "initialization", "K",
                       "modelx", "modely", "BIC",
                       paste0("C", 1:ncol(traffic_fd$b2x$coefs)))
  write.csv(temp_res, file = paste0("traffic_", cutx, "-", nbasis, "res.csv"))
}


