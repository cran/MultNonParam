% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mood.median.test.R
\name{mood.median.test}
\alias{mood.median.test}
\title{Mood's Median test, extended to odd sample sizes.}
\usage{
mood.median.test(x, y, exact = FALSE)
}
\arguments{
\item{x}{First data set.}

\item{y}{Second data set.}

\item{exact}{Indicator for whether the test should be done exactly or approximately.}
}
\value{
The two-sided p-value.
}
\description{
Test whether two samples come from the same distribution.  This version of Mood's median test is presented for pedagogical purposes only.  Many authors successfully argue that it is not very powerful.  The name "median test" is a misnomer, in that the null hypothesis is equality of distributions, and not just equality of median.  Exact calculations are not optimal for the odd sample size case.
}
\details{
The exact case reduces to Fisher's exact test.
}
