\name{wilding}
\alias{wilding}
\title{Two-stage Mann-Whitney test}
\description{Function that implements a two-stage Mann-Whitney test as described in the G. E. Wilding, G. Shan and A. D. Hutson paper (see references).}

\usage{wilding(u1, u2, m1, n1, m2, n2)}

\arguments{
  \item{u1, u2}{The asymptotic null distribution of first and second stage test statistics}
  \item{n1, n2}{Number of Y observations for stage 1 and 2.}
  \item{m1, m2}{Number of X observations for stage 1 and 2.}

}

\value{The exact two-stage Mann-Whitney probability distribution.}

\author{John E. Kolassa}
\references{G. E. Wilding, G. Shan and A. D. Hutson (2010), \emph{Exact Two-Stage Designs for Phase II Clinical Trials with Rank-Based Endpoints}.\cr
S.S Brown, M.J.R. Healy, M. Kearns  (1981), \emph{Report on the inter-laboratory trial of the reference method for the determination of the total calcium in serum}, J.Clin. Chem. Clin. Biochem. 19, 395-426.
}

\examples{
## From Brown-Healy-Kearns (1981), data concerning calcium measurments.
df.calcium <- data.frame(x1=c(1206,1202,1202,1201)
					, y1=c(1124, 1133, 1139, 1133)
					, x2=c(1017, 1017, 1012, 1020)
					, y2=c(965, 961, 960, 962))
attach(df.calcium)
u1 <- wilcox.test(x1,y1,exact=FALSE)$statistic
u2 <- wilcox.test(c(x1,x2),c(y1,y2),exact=FALSE)$statistic
wilding(u1[[1]],u2[[1]],length(x1),length(y1),length(x2), length(y2))
}
