\name{sensitivity.plot}
\alias{sensitivity.plot}

\title{Sensitivity plot}

\description{Compare the sensitivity of different statistics.}

\usage{sensitivity.plot(y,sub,stats)}

\arguments{
\item{y}{vector of the data.}
\item{sub}{subtitle for the plot.}
\item{stats}{vector of functions to be plotted.}
}

\details{To compare the sensitivity, outliers are added to the original data. The shift of each statistics due to the new value is measured and plotted.}

%\references{}

\author{John E. Kolassa}

\examples{
location<-function(y){out<-c(mean(y),median(y),mean(y,trim=.1))
    names(out)<-c("mean","median","trim mean"); return(out)}
# Draw a sensitivity plot comparing the performance of mean, median,
# and trimmed mean for the DNase rat DNA density data.  The data set
# is in the R datasets package.
sensitivity.plot(DNase[DNase$Run==1,3],"DNase Run 1",location)
}