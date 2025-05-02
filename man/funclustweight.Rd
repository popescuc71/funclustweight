\name{funclustweight}
\alias{funclustweight}
\title{
funclustweight: Functional Regression Clustering.
}
\description{
funclustweight is ...
}
\usage{
funclustweight(datax, datay, K=1:10, model="AkjBkQkDk", modely = "VII", known=NULL,
           threshold=0.1, itermax=200, eps=1e-6, init='random',
           criterion="bic", d_select="Cattell", init.vector=NULL,
           show=TRUE, mini.nb=c(5, 10), min.individuals=2, mc.cores=1, nb.rep=2,
           keepAllRes=TRUE, kmeans.control = list(), d_max=100, d_range=2,
           verbose = TRUE)
}
\arguments{
  \item{data}{
In the univariate case: a functional data object produced by the fda package. In the multivariate case: a list of functional data objects. Represents the X part of the data used to predict datay.
}
  \item{datay}{
In the univariate case: a functional data object produced by the fda package. In the multivariate case: a list of functional data objects. Represents the Y part of the data predicted by data.
}
  \item{K}{
The number of clusters or list of clusters to try, for example K=2:10.
}
  \item{model}{
The chosen model among 'AkjBkQkDk', 'AkjBQkDk', 'AkBkQkDk', 'ABkQkDk', 'AkBQkDk', 'ABQkDk'. 'AkjBkQkDk' is the default. We can test multiple models at the same time with the command c(). For example c("AkjBkQkDk","AkjBQkDk").
}
  \item{modely}{
The chosen model among "EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "VEE", "EVE", "EEV", "VVE", "VEV","EVV", "VVV". 'VII' is the default. We can test multiple models at the same time with the command c(). For example c("VII","EII").
}
  \item{threshold}{
The threshold of the Cattell' scree-test used for selecting the group-specific intrinsic dimensions.
}
\item{known}{
A vector of known classifications that can be numeric or \code{NA}. It is optional for clustering. For classification, curves with unknown classification should be given the value \code{NA} within \code{known} (see the examples below). Must be the same length as the number of curves in the data set.
}
  \item{itermax}{
The maximum number of iterations.
}
  \item{eps}{
The threshold of the convergence criterion.
}
  \item{init}{
A character string. It is the way to initialize the EM algorithm. There are five ways of initialization: “kmeans” (default), “param”, “random”, “mini-em”, “vector”, or "tkmeans". See details for more information. It can also be directly initialized with a vector containing the prior classes of the observations.
}
  \item{criterion}{
The criterion used for model selection: bic (default) or icl.
}
  \item{d_select}{
“Cattell” (default), “BIC”, or "grid". This parameter selects which method to use to select the intrinsic dimensions of subgroups. "grid" will select d based on thecriterion value after running each combination of d1, d2, ..., dK for the groups. d used for each group is based on the values for d_range. "grid" will only work for a single value of K (not a list). See details for more information.
}
  \item{init.vector}{
A vector of integers or factors. It is a user-given initialization. It should be of the same length as of the data. Only used when init="vector".
}
  \item{show}{
Use show = FALSE to settle off the informations that may be printed.
}
  \item{mini.nb}{
A vector of integers of length two. This parameter is used in the “mini-em” initialization. The first integer sets how many times the algorithm is repeated; the second sets the maximum number of iterations the algorithm will do each time. For example, if init=“mini-em” and mini.nb=c(5,10), the algorithm wil be launched 5 times, doing each time 10 iterations; finally the algorithm will begin with the initialization that maximizes the log-likelihood.
}
  \item{min.individuals}{
This parameter is used to control for the minimum population of a class. If the population of a class becomes stricly inferior to 'min.individuals' then the algorithm stops and gives the message: 'pop<min.indiv.'. Here the meaning of "population of a class" is the sum of its posterior probabilities. The value of 'min.individuals' cannot be lower than 2.
}
  \item{mc.cores}{
Positive integer, default is 1. If mc.cores>1, then parallel computing is used, using mc.cores cores. Warning for Windows users only: the parallel computing can sometimes be slower than using one single core (due to how parLapply works).
}
  \item{nb.rep}{
A positive integer (default is 1 for kmeans initialization and 20 for random initialization). Each estimation (i.e. combination of (model, K, threshold)) is repeated nb.rep times and only the estimation with the highest log-likelihood is kept.
}
  \item{keepAllRes}{
Logical. Should the results of all runs be kept? If so, an argument all_results is created in the results. Default is TRUE.
}
  \item{kmeans.control}{
A list. The elements of this list should match the parameters of the \code{\link{kmeans}} initialization (see kmeans help for details). The parameters are “iter.max”, “nstart” and “algorithm”. "alpha" is an added parameter for the \code{\link{tkmeans}} initialization (see tkmeans help for details)
}
  \item{d_max}{
A positive integer. The maximum number of dimensions to be computed. Default is 100. It means that the instrinsic dimension of any cluster cannot be larger than d_max. It quickens a lot the algorithm for datasets with a large number of variables (e.g. thousands).
}
  \item{d_range}{
Vector of values to use for the intrinsic dimension for each group when d_select="grid".
}

  \item{verbose}{
Whether to print progress and approximate timing information as funclustweight executes. TRUE (default when running in serial) or FALSE (default when running parallel).
}
}
\details{
If we choose init="random", the algorithm is run 20 times with the same model options and the solution which maximises the log-likelihood is printed. This explains why sometimes with this initialization it runs a bit slower than with 'kmeans' initialization.

If the warning message: "In funclustweight(...) : All models diverged" is printed, it means that the algorithm found less classes that the chosen number (parameter K). Because the EM algorithm is used, it could be because of a bad initialization of the EM algorithm. So we have to restart the algorithm multiple times in order to check if with a new initialization of the EM algorithm the model converges, or if there is no solution with the chosen number K.


The different initializations are:

“mini-em”:
it is an initialization strategy for which the classes are randomly initialized and the EM algorithm is run for several iterations.  This action is repetead a few times (the default is 5 iterations and 10 times). At the end, the initialization chosen is the one which maximise the log-likelihood (see mini.nb for more information about its parameters).

“random”:
the classes are randomly given using a multinomial distribution

“kmeans”:
the classes are initialized using the kmeans function (with algorithm="Hartigan-Wong"; nstart=4; iter.max=50); note that the user can use his own arguments for kmeans using the dot-dot-dot argument

A prior class "vector":
It can also be directly initialized with a vector containing the prior classes of the observations. To do so use init="vector" and provide the vector in the argument init.vector.

Note that the BIC criterion used in this function is to be maximized and is defined as 2*LL-k*log(n) where LL is the log-likelihood, k is the number of parameters and n is the number of observations.


There are three methods for selecting the intrinsic dimension using d_select:

"Cattell":
Runs a Cattell's scree test to approximate the intrinsic dimension that yields the greatest improvement in clustering.

"BIC":
At each iteration we tests each value for each group's intrinsic dimension and sets the intrinsic dimension that yields the best BIC.

"grid":
Runs every combination of hyperparameters (eg. K=2, threshold = 0.05, model = ...) for every combination of intrinsic dimensions that can be set with the given d_range (with K = 2 and d_range = c(2, 10) it would set (2,2), (2, 10), (10, 2), and (10, 10)). Due to the sharp increase in test cases it is recommended that this mode is run in parallel if possible. Doing an intial short run to approximate the timing with verbose = TRUE is suggested as well.
}
\value{
  \item{d}{The number of dimensions for each cluster.}
  \item{a}{Values of parameter a for each cluster.}
  \item{b}{Values of parameter b for each cluster.}
  \item{model}{Model used for clustering X part of data.}
  \item{modely}{Model used for clustering regression coefficients of data.}
  \item{mu}{The mean of each cluster in the original space.}
  \item{prop}{The proportion of individuals in each cluster.}
  \item{loglik}{The maximum of log-likelihood.}
  \item{loglik_all}{The log-likelihood at each iteration.}
  \item{posterior}{The posterior probability for each individual to belong to each cluster.}
  \item{class}{The clustering partition.}
  \item{BIC}{The BIC value.}
  \item{ICL}{The ICL value.}
  \item{complexity}{the number of parameters that are estimated.}
  \item{all_results}{if multiple number of clusters or models are considered, results for each model are stored here}
  \item{gam}{The regression coefficients for each cluster.}
  \item{covy}{Covariance matrix of the clustered data estimated by the model.},
}
\references{
- C.Bouveyron and J.Jacques, Model-based Clustering of Time Series in Group-specific Functional Subspaces, Advances in Data Analysis and Classification, vol. 5 (4), pp. 281-300, 2011 <doi:10.1007/s11634-011-0095-6>

- Schmutz A, Jacques J, Bouveyron C, et al (2020) Clustering multivariate functional data in group-specific functional subspaces. Comput Stat 35:1101-1131

- Cristina Anton, Iain Smith ?

-Dang U J, Punzo A, et al (2017) Multivariate response and parsimony for Gaussian cluster-weighted models. Journal of Classification 34:4-34.

-Celeux G. and Govaert G. (1995). Gaussian Parsimonious Clustering Models. Pattern Recognition.
28(5), 781–793.

-Browne R. P. and McNicholas P. D. (2015). mixture: Mixture Models for Clustering and Classification. R package version 1.4.
}
\author{
Cristina Anton and Iain Smith
}
\seealso{
\code{\link{kmeans}}
}

\examples{
\donttest{
## ADELAIDE EXAMPLE
demand <- fitAdelaideFD()
models <- c("AKJBKQKDK","AKJBQKDK", "AKBKQKDK", "ABKQKDK", "AKBQKDK", "ABQKDK")
modelsy <- c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "VEE", "EVE",
             "EEV", "VVE", "VEV","EVV","VVV")
modelsys <- c("EII", "VII", "EEI", "VEI")
res_reg <- funclustweight(demand$fdx, demand$fdy, K=2, model = models,
                          modely = modelsys, init = "kmeans", nb.rep = 1,
                          threshold = 0.1)
set.seed(1027)
res_reg = funclustweight(demand$fdx, demand$fdy, K = 2, model = "AKJBKQKDK", modely = "VII", nb.rep = 1, verbose = T)
table(demand$groupd, res_reg$class)

## MODEL DATA
set.seed(10403)
model_demand <- genFromModelRegFD(res_reg, ncurves = c(300, 300), mu_div = 1)
res_reg = funclustweight(model_demand$xfd, model_demand$yfd, K = 2, model = models, modely = "VII", nb.rep = 5, verbose = T)

}
}
\keyword{cluster}
