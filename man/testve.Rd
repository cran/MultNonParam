\name{testve}
\alias{testve}
\title{Diagnosis for multivariate stratified Kawaguchi - Koch - Wang method}
\description{
Diagnostic tool that verifies the normality of the estimates of the probabilities b with the Kawaguchi - Koch - Wang method. The diagnostic method is based on a Monte Carlo method.
}
\usage{
testve(n, m, k, nsamp = 100, delta = 0, beta = 0, disc = 0)
}
\arguments{
  \item{n}{number of observations in the first group.}
  \item{m}{number of observations in the second group.}
  \item{k}{number of stratum.}
  \item{nsamp}{The number of estimates that will be calculated. Must be enough to be sure that the results are interpretable. Default value is 10000.}
  \item{delta}{Offset that depends on group.}
  \item{beta}{Correlation between \code{x} and \code{y}.}
  \item{disc}{The Mann Whitney test is designed to handle continuous data, but this method applies to discretized data; \code{disc} adjusts the discreteness.}


}
\details{
This functions serves as a diagnosis to prove that the Kawaguchi - Koch - Wang method gives Gaussian estimates for b. It generates random data sets, to which the Mann Whitney test gets applied.
\code{y} is the generated response variable and \code{x} the generated covariable related to \code{y} through a regression model.
}
\value{
The function returns a normal QQ plot.
}
\references{
A. Kawaguchi, G. G. Koch and X. Wang (2012), \emph{Stratified Multivariate Mann-Whitney Estimators for the Comparison of Two Treatments with Randomization Based Covariance Adjustment}.\cr
J. E. Kolassa and Y. Seifu (2013), \emph{Nonparametric Multivariate Inference on Shift Parameters }, Academic Radiology vol.20 n7, 883-888.
}

\author{
John E. Kolassa
}

\seealso{\code{\link{probest}} is the function used to calculate the estimates with the Kawaguchi - Koch - Wang method.}

\examples{
#Dimensions of the data
n <- 10
m <- 15
k = 3
nsamp <- 100
testve(n,m,k,nsamp, beta=0.4)
}

