#******************************* REQUIRED LIBRARY *****************************#
library(fda)
library(fds)


fitAdelaideFD <- function() {
  #****** open required data ******#
  # data used in the paper
  de_sun <- fds::sundaydemand
  # data we chose as a different cluster demand exists for any weekday
  de_tue <- fds::tuesdaydemand

  #****** setup labels for days ******#
  # when checking dimensions using the $x nrow is timepoints
  cls <- append(rep("sunday", ncol(de_sun$y)),
                rep("tuesday", ncol(de_tue$y)))

  #****** combine and format data ******#
  demand <- cbind(de_sun$y, de_tue$y)

  # potentially adjust where the data cuts off
  # I set labels 1 to 24 for 24 hours, original data is in half-hour chunks
  # so the original data is 48 hours
  bbasis_x <- create.bspline.basis(rangeval = c(1, 12), nbasis = 6)
  # total splines is 24 so I use 18 and 6 to represent them
  bbasis_y <- create.bspline.basis(rangeval = c(13, 24), nbasis = 6)
  bbasis_fullx <- create.bspline.basis(rangeval = c(1, 24), nbasis = 6)
  # splits are made so no data is shared
  # rows represent recrodings per half-hours it is already transposed
  fd_demand_x <- smooth.basis(argvals = seq(1, 12, length.out = 24),
                              y = demand[1:24,], fdParobj = bbasis_x)$fd
  fd_demand_y <- smooth.basis(argvals = seq(13, 24, length.out = 24),
                              y = demand[25:48,], fdParobj = bbasis_y)$fd
  fd_demand_fullx <- smooth.basis(argvals = c(seq(1, 12, length.out = 24), seq(13, 24, length.out = 24)),
                              y = demand[1:48,], fdParobj = bbasis_fullx)$fd
  return( list(fdx = fd_demand_x,
               fdy = fd_demand_y,
               fdfull=fd_demand_fullx,
               groupd = cls) )
}

plotAdelaideFD <- function(fd_demand) {
  win.graph(width=9.5, height=6,pointsize=12)
  par(mfrow=c(1,2))
  # assuming 2 labels of saturday and other
  plot.fd(fd_demand$fdx, col = ifelse(fd_demand$groupd=="tuesday", "black", "red"))
  plot.fd(fd_demand$fdy, col = ifelse(fd_demand$groupd=="tuesday", "black", "red"))
  legend("topright", c("tuesday","sunday"), fill = c("black", "red"))

  win.graph(width=9.5, height=6,pointsize=12)
  par(mfrow=c(1,1))
  plot.fd(fd_demand$fdfull, col = ifelse(fd_demand$groupd=="tuesday", "black", "red"))
  win.graph(width=9.5, height=6,pointsize=12)
  par(mfrow=c(1,2))
  plot.fd(fd_demand$fdfull[fd_demand$groupd=="tuesday",], col = "red")
  plot.fd(fd_demand$fdfull[fd_demand$groupd=="sunday",], col = "black")
  win.graph(width=9.5, height=6,pointsize=12)
  par(mfrow=c(1,2))
  plot(fds::sundaydemand, main="Sunday",xaxt = "n")
  abline(v = 24)
  axis(1, at = c(0, 6, 12, 18, 24, 30, 36, 42, 48), labels=c("0", "3", "6", "9", "12", "15", "18", "21", "24"))
  plot(fds::tuesdaydemand, main="Tuesday",xaxt = "n")
  abline(v = 24)
  axis(1, at = c(0, 6, 12, 18, 24, 30, 36, 42, 48))
  }
