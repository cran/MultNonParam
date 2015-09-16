\name{higgins.fisher.kruskal.test}
\alias{higgins.fisher.kruskal.test}

\title{Fisher's Least Significant Difference (LSD) for Kruskal-Wallis test}
\description{Rank-based method for controlling experiment-wise error. Two hypothesis have to be respected: normality of the distribution and no ties in the data. The aim is to be able to detect, among k treatments, those who lead to significant differencies in the values for a variable of interest.
}
\usage{higgins.fisher.kruskal.test(resp, grp, alpha=.05)}

\arguments{
  \item{resp}{vector containing the values for the variable of interest.}
  \item{grp}{vector specifying in which group is each observation.}
  \item{alpha}{level of the test.}
  }

\details{First, the Kruskal-Wallis test is used to test the equality of the distributions of each treatment. If the test is significant at the level \code{alpha}, the method can be applied.
}
\value{
  A matrix with two columns. Each row indicates a combinaison of two groups that have significant different distributions.
}

\references{ J.J. Higgins, (2004), \emph{Introduction to Modern Nonparametric Statistics}, Brooks/Cole, Cengage Learning.
}

\author{John E. Kolassa}

\examples{
#Randomly generated data
resp <- c(rnorm(15,0,1),rnorm(15,50,1), rnorm(15,50,1))
grp <- c(rep(1,15),rep(2,15),rep(3,15))
higgins.fisher.kruskal.test(resp, grp)
#We expect to have 1-3 and 1-2 significantly different 
#and not 2-3 since they are generated the same way (so same distribution).
}
