\name{aov.P}
\alias{aov.P}

\title{One-way ANOVA using permutation tests}
\description{
\code{aov.P} uses permutation tests instead of classic theory tests to run a one-way ANOVA.
}
\usage{aov.P(dattab)}

\arguments{
  \item{dattab}{ The table on which the ANOVA has to be done.}
}

\details{
The function calls a Fortran code to perform the permutation tests and the ANOVA. 
The function has to be applied directly on a cross-table of two variables.

}
\value{
  \item{pv}{the p-value obtained with the permutation tests.}
  \item{tot}{}
}

\references{ J.J. Higgins, (2004), \emph{Introduction to Modern Nonparametric Statistics}, Brooks/Cole, Cengage Learning.
}

\author{John E. Kolassa}

\examples{
##Building an array for the example. The data is arbitrary.
dattab<-array(c(1.5,1.8,1.9,2.1,2.0,2.5,1.9,2.0,2.5
              ,2.8,2.7,2.6,1.4,1.6,2.1,1.8,2.3,2.4)
              ,dim=c(3,6)
              )
aov.P(dattab)
}