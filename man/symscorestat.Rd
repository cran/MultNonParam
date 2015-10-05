\name{symscorestat}
\alias{symscorestat}

\title{Generalization of Wilcoxon signed rank test}
\description{
\code{symscorestat} returns either exact or asymptotic p-values for score tests of the null hypothesis of univariate symmetry about 0.
}
\usage{symscorestat(y,scores,exact,sides)}

\arguments{
  \item{y}{Vector of data on which test will be run.}
  \item{scores}{Scores to be used for the test.  Defaults to integers 1:length(y).}
  \item{exact}{ Logical variable indicating whether the exact p-value should be calculate.  Default is false.}
  \item{sides}{ Integer; 1 for one sided test rejecting for large values of the statistic, and 2 for the two-sided test.  Defaults to 1.}
}

\details{
The statistic considered here is the sum of scores corresponding to those entries in y that are positive.
If exact=T, the function calls a Fortran code to cycle through all permutations.  If exact=F, the expectation of the statistic is calculated as half the sum of the scores, the variance is calculated as one quarter the sum of squares of scores about their mean, and the statistic is compared to its approximating normal distribution.

}
\value{
  \item{pv}{the p-value obtained with the permutation tests.}
  \item{tot}{the total number of rearrangements of the data considred in calculating the p-value.}
}

\references{ J.J. Higgins, (2004), \emph{Introduction to Modern Nonparametric Statistics}, Brooks/Cole, Cengage Learning.
}

\author{John E. Kolassa}

\examples{
symscorestat(y=c(1,-2,3,-4,5),exact=TRUE)
}
