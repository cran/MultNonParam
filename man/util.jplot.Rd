\name{util.jplot}
\alias{util.jplot}
\title{Smoothered plot}
\description{Graphical function that takes a noisy version of a piecewise constant function and plots a smoothered version.}

\usage{util.jplot(x, y, ...)}

\arguments{
  \item{x,y}{Vectors of ordinate and abscissa for the points to be ploted.}

  \item{\dots}{R graphical options to be added.}
}

\value{A plot where steps are only represented by segments instead of noisy scattered points.}

\author{John E. Kolassa}


\examples{

##Randomly generated data is used.
eps <- rnorm(300, 0,0.001) 
#eps is added in order to have slightly different values of y for each point.
y <- c(rep(1,100),rep(2,110),rep(3,90))+eps
x <- seq(from=0,to=1, length.out=300)

##Graphs
par(mfrow=c(2,1))
#Normal plot of x by y
plot(x,y)
#util.jplot of x by y
util.jplot(x,y)
}
