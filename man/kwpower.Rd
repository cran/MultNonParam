% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kwpower.R
\name{kwpower}
\alias{kwpower}
\title{Power for the Kruskal-Wallis test.}
\usage{
kwpower(
  nreps,
  shifts,
  distname = c("normal", "cauchy", "logistic"),
  level = 0.05,
  mc = 0,
  taylor = FALSE
)
}
\arguments{
\item{nreps}{The numbers in each group.}

\item{shifts}{The offsets for the various populations, under the alternative hypothesis.}

\item{distname}{The distribution of the underlying observations; normal, cauchy, and logistic are currently supported.}

\item{level}{The test level.}

\item{mc}{0 for asymptotic calculation, or positive for mc approximation.}

\item{taylor}{logical determining whether Taylor series approximation is used for probabilities.}
}
\value{
A list with components power, giving the power approximation, ncp, giving the noncentrality parameter, cv, giving the critical value, probs, giving the intermediate output from pairwiseprobability, and expect, the quantities summed before squaring in the noncentrality parameter.
}
\description{
\code{kwpower} approximates power for the Kruskal-Wallis test, 
using a chi-square approximation under the null, and a non-central chi-square approximation under the alternative.  The noncentrality parameter is calculated using alternative means and the null variance structure.
}
\details{
The standard noncentral chi-square power formula, or Monte Carlo, is used.
}
\examples{
#Calculate the power for the Kruskal Wallis test for normal observations,
#10 observations in each of three groups, with groups centered at 0, 1, 2.
#Level is 0.05 by default.
kwpower(rep(10,3),c(0,1,2),"normal")
}
