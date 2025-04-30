\name{fitAdelaideFD}
\alias{fitAdelaideFD}
\title{
fitAdelaideFD
}
\description{
Extract Adelaide data from fds
}
\arguments{
  \item{nbasis}{
The number of basis functions to fit to the data.
}
}

\usage{
fitAdelaideFD()
}
\details{
Open Adelaide electricity demand data from fds for Sunday and Tuesday. Fit the data to 6 basis functions for the X and Y sides, splitting the data at noon mid-day.
}
\value{
  \item{fd}{A functional data object representing the fitted Adelaide Demand data.}
  \item{groupd}{Group classifications for each curve as a Sunday or Tuesday demand curve.}
}
\references{
- citation for paper labelled OG-FNN

- Cristina Anton, Iain Smith .
}
\author{
Cristina Anton and Iain Smith
}
\seealso{
\code{\link{plotAdelaideFD}}
}
\examples{
data1<-fitAdelaideFD()
plotAdelaideFD(data1)
}
\keyword{datagen}
