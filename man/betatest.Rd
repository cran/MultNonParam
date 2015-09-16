\name{betatest}
\alias{betatest}

\title{Permumation p-value for the test of association between two variables.}
\description{Calculate the p-value for the test of association between two variables using the permutation method.}

\usage{betatest(x, y)}

\arguments{
  \item{x, y}{ vectors for the two variables to test.}
}

\details{Based on a Fortran routine.}

\value{Returns a p-value.}

\author{John E. Kolassa}

\examples{
#Example using data from plant Qn1 from the CO2 data set.
betatest(CO2[CO2$Plant=="Qn1",4],CO2[CO2$Plant=="Qn1",5])
}