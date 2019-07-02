\name{qif}
\alias{qif}
\alias{print.qif}
\alias{summary.qif}
\alias{print.summary.qif}
\title{
  Function to solve a Quadratic Inference Function Model
}
\description{
  Produces an object of class \code{"qif"} which is a Quadratic Inference Function fit of the balanced longitudinal data.
}
\usage{
qif(formula = formula(data), id = id, data = parent.frame(), b = NULL, tol = 1e-8, maxiter = 1000, family = gaussian, corstr = "independence", invfun = "finv")
}
\arguments{
\item{formula}{
  a formula expression as for other regression models, of the form
  \code{response ~ predictors}. See the documentation of \code{\link{lm}} and \code{\link{formula}} for details.
}
\item{id}{
  a vector which identifies the clusters.  The length of \code{id} should be
  the same as the number of observations.  Data are assumed to be sorted
  so that observations on a cluster are contiguous rows for all entities
  in the formula.
}
\item{data}{
  an optional data frame in which to interpret the variables occurring
  in the \code{formula}, along with the \code{id} variables.
}
\item{b}{
  an initial estimate for the parameters.
}
\item{tol}{
  the tolerance used in the fitting algorithm.
}
\item{maxiter}{
  the maximum number of iterations.
}
\item{family}{
  a \code{family} object: a list of functions and expressions for
  defining canonical link and variance functions.  Families supported
  in \code{qif} are \code{gaussian}, \code{binomial}, \code{poisson},
  and \code{gamma};
  see the \code{\link{glm}} and \code{\link{family}} documentation.
  Some links are not currently available: probit link for binomial family   and  log link for gamma family.
}
\item{corstr}{
  a character string specifying the correlation structure.
  The following are permitted:
     \code{"independence"},
     \code{"exchangeable"},
     \code{"AR-1"} and
     \code{"unstructured"}.
}
\item{invfun}{
  a character string specifying the matrix inverse function.
  The following are permitted:
     \code{"finv"} and
     \code{"ginv"}.
}
}
\value{
  A list containing:
  \item{title}{name of qif}
  \item{version}{the current version of qif}
  \item{model}{ analysis model for link function, variance function and correlation struture}
  \item{terms}{ analysis model for link function, variance function and correlation struture}
  \item{iteration}{the number of iterations}
  \item{coefficients}{beta esitmates value}
  \item{linear.perdictors}{linear predictor value}
  \item{fitted.value}{fitted value of y}
  \item{x}{the perdicted matrix}
  \item{y}{the response}
  \item{residuals}{y-mu}
  \item{pearson.resi}{pearson residuals}
  \item{scale}{the scale of fitted model}
  \item{family}{the type of distribution}
  \item{id}{model fitted value}
  \item{max.id}{max number of each steps}
  \item{xnames}{the values are X name of qif}
  \item{statistics}{The qif statistics}
  \item{Xnames}{the name X matrix in qif}
  \item{parameter}{parameter estimates and}
  \item{covariance}{Covariance of coefficients}
}
\note{
  This R package is created by transplanting a SAS macro QIF developed originally by Peter Song and Zhichang Jiang (2006).
  This is version 1.0 of this user documentation file, revised
  30/07/2009.
}
\author{
  Zhichang Jiang, Alberta Health Services, and Peter X.K. Song, University of Michigan.
}
\references{
  Qu A, Lindsay BG, Li B. Improving generalized estimating equations using quadratic inference functions. 
  \emph{Biometrika 2000}, \bold{87} 823-836. 

  Qu A, Song P X-K. Assessing robustness of generalised estimating equations and quadratic inference functions. 
  \emph{Biometrika 2004}, \bold{91} 447-459.

  Qu A, Lindsay BG. Building adaptive estimating equations when inverse of covariance estimation is difficult. 
  \emph{J. Roy. Statist. Soc. B 2003}, \bold{65}, 127-142.
}
\seealso{
  \code{\link{glm}}, \code{\link{lm}}, \code{\link{formula}}.
}
\examples{
## Marginal log-linear model for the epileptic seizures count data (Diggle et al., 2002, Analysis of Longitudinal Data, 2nd Ed., Oxford Press).
# fit <- qif(y ~ baseline + trt + logage + visit, id=id, data=epilepsy, family=poisson, corstr="AR-1")
# summary(fit)
}
\keyword{nonlinear}

