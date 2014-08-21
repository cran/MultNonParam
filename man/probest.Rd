\name{probest}
\alias{probest}

\title{Stratified Multivariate Kawaguchi Koch Wang  Estimators}
\description{
Function that return the estimators and their variance-covariance matrix calculated with the Kawaguchi - Koch - Wang method.
}
\usage{probest(ds, resp, grp, str, covs = NULL, delta = NA)}

\arguments{
  \item{ds}{The \code{data.frame} to be used.}
  \item{resp}{The vector of the response manifest variable. There can be more than one variable. It has to be the name of the variable as a character string.}
  \item{grp}{The vector of the variable that divides the population into groups. It has to be the name of the variable as a character string.}
  \item{str}{The vector of the variable used for the strata. It has to be the name of the variable as a character string.}
  \item{covs}{The covariates to be used in the model. It has to be the name of the variable as a character string. Default value is \code{NULL}}
  \item{delta}{Default value is \code{NA}.}
}

\details{
The function calls a Fortran code to calculate the estimators \code{b} and their variance-covariance matrix \code{Vb}.

}
\value{

\item{b}{vector of adjusted estimates from the method}
\item{Vb}{the corresponding estimated covariance matrix}
}

\references{A. Kawaguchi, G. G. Koch and X. Wang (2012), \emph{Stratified Multivariate Mann-Whitney Estimators for the Comparison of Two Treatments with Randomization Based Covariance Adjustment}.\cr
J. E. Kolassa and Y. Seifu (2013), \emph{Nonparametric Multivariate Inference on Shift Parameters }, Academic Radiology vol.20 n7, 883-888.

}

\author{John E. Kolassa}

\examples{
##We use here the data about breast cancer from the MultNonParam package.
data(sotiriou.data)
attach(sotiriou.data)

##First simple plot of the data
plot(AGE,TUMOR_SIZE,pch=(recur+1),
     main="Age and Tumor Size",
     sub="Breast Cancer Recurrence Data",
     xlab="Age (years)",ylab="Tumor Size",col=c("blue","darkolivegreen"))
legend(31,8,legend=c("Not Recurrent","Recurrent"),pch=1:2,col=c("blue","darkolivegreen"))


##AGE and TUMOR_SIZE are the response variables, recur is used for the groups, 
##TAMOXIFEN_TREATMENT for the stratum and ELSTON.ELLIS_GRADE is a covariable.
po<-probest(sotiriou.data
			,c("AGE","TUMOR_SIZE")
			,"recur"
			,"TAMOXIFEN_TREATMENT"
			,"ELSTON.ELLIS_GRADE"
			)
}
