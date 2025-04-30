\name{genFromModelRegFD}
\alias{genFromModelRegFD}
\title{
genFromModelRegFD
}
\description{
Generate univariate functional data for X with coefficients distributed according to a finite mixture of normal distributions.
Also generate univariate functional data for Y with coefficients determined by X using simple regression coefficients. Parameters
for the model data are loaded from a file.
}
\arguments{
  \item{res_reg}{
An object of type funclustweightas the ones resulting from clustering data with two groups.
}
  \item{ncurves}{
A vector denoting number of curves for each of the two groups.
}
  \item{nsplines}{
The number of splines to fit to the simulated data.
}

  \item{mu_div}{
A number to divide mu value by to increase complexity of the problem.
}
}

\usage{
genModelFD(ncurves=c(200, 300), nsplines=12)
}
\details{
The number of clusters is fixed to \eqn{K=2}. 
}
\value{
  \item{xfd}{A functional data object representing the X part of simulated data.}
  \item{yfd}{A functional data object representing the Y part of simulated data.}
  \item{groupd}{Group classifications for each curve.}
}
\references{
- Cristina Anton, Iain Smith ?
}
\author{
Cristina Anton and Iain Smith
}

\examples{
res_reg=readRDS("input1n.RData")
model_reg <- genFromModelRegFD(res_reg, ncurves = c(300, 300),mu_div = 1)
plotRegFD(model_reg)
}
\keyword{datagen}
