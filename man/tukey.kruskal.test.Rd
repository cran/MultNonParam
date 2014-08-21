\name{tukey.kruskal.test}
\alias{tukey.kruskal.test}

\title{Tukey HSD procedure}
\description{Rank-based method for controlling experiment-wise error. We make the hypothesis of the normality of the distribution for the variable of interest. 
}
\usage{tukey.kruskal.test(x, y, alpha=.05)}

\arguments{
  \item{x}{vector specifying in which group is each observation.}
  \item{y}{vector containing the values for the variable of interest.}
  \item{alpha}{level of the test.}
  }

\details{The original Tuckey HSD procedure is supposed to be applied for equal sample sizes. However, the \code{tukey.kruskal.test} function performs the Tukey-Kramer procedure that works for unequal sample sizes.
 }

\value{A logical vector for every combinaison of two groups. \code{TRUE} if the distribution in one group is significantly different from the distribution in the other group. 
}

\references{ J.J. Higgins, (2004), \emph{Introduction to Modern Nonparametric Statistics}, Brooks/Cole, Cengage Learning.
}

\author{John E. Kolassa}

\examples{
#Example using data from plant Qn1 from the CO2 data set.
tukey.kruskal.test(DNase[,3],DNase[,1])
#Here no group seem to have a different distribution from another one.
}