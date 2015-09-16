\name{MultNonParam-package}
\alias{MultNonParam-package}
\alias{MultNonParam}
\docType{package}
\title{MultNonParam}
\description{A collection of multivariate nonparametric methods.}
\details{
\tabular{ll}{
Package: \tab MultNonParam\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-06-03\cr
License: \tab GPL-2\cr
}
% An overview of how to use the package, including the most important functions 
}
\author{
John E. Kolassa & Stephane Jankowski

Maintainer: John E. Kolassa <kolassa@stat.rutgers.edu>
}
\references{
A. Kawaguchi, G. G. Koch and X. Wang (2012), \emph{Stratified Multivariate Mann-Whitney Estimators for the Comparison of Two Treatments with Randomization Based Covariance Adjustment}.\cr
S. Madhavan, Y. Gusev, M. Harris, D. Tanenbaum, R. Gauba, K. Bhuvaneshwar, A. Shinohara, K. Rosso, L. Carabet, L. Song, R. Riggins, S. Dakshanamurthy, Y. Wang, S. Byers, R. Clarke, L. Weiner (2011), \emph{A systems medicine platform for personalized oncology}, Neoplasia 13. \cr
A. V. D'Amico, R. Whittington, S. B. Malkowicz, D. Schultz, K. Blank, G. A. Broderick, J. E. Tomaszewski, A. A. Renshaw, I. Kaplan, C. J. Beard, A. Wein (1998) , \emph{Biochemical outcome after radical prostatectomy, external beam radiation therapy, or interstitial radiation therapy for clinically localized prostate cancer}, JAMA : the journal of the American Medical Association 280 969-74. \cr
J. E. Kolassa and Y. Seifu (2013), \emph{Nonparametric Multivariate Inference on Shift Parameters }, Academic Radiology vol.20 n7, 883-888. \cr
C. Sotiriou, P. Wirapati, S. Loi, A. Harris, S. Fox, J. Smeds, H. Nordgren, P. Farmer, V. Praz, B. Haibe-Kains, C. Desmedt, D. Larsimont, F. Cardoso, H. Peterse, D. Nuyten, M. Buyse, M. Van de Vijver, J. Bergh, M. Piccart, M. Delorenzi  (2006), \emph{Gene expression profiling in breast cancer: understanding the molecular basis of histologic grade to improve prognosis}, Journal of the National Cancer Institute 98 262-72. \cr
G. E. Wilding, G. Shan and A. D. Hutson (2010), \emph{Exact Two-Stage Designs for Phase II Clinical Trials with Rank-Based Endpoints}.\cr
S.S Brown, M.J.R. Healy, M. Kearns  (1981), \emph{Report on the inter-laboratory trial of the reference method for the determination of the total calcium in serum}, J.Clin. Chem. Clin. Biochem. 19, 395-426. \cr
J.J. Higgins, (2004), \emph{Introduction to Modern Nonparametric Statistics}, Brooks/Cole, Cengage Learning.
}

\keyword{ package }

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
#Code below produces a confidence region.  
#Example below is crude; a better region would be produced with the commented-out dimensions.
#dl<-list(.5+(-60:60)/10,.6+(-50:50)/120)
dl<-list(.5+2*(-3:3),.6+(-3:3)/6)
statv<-array(NA,unlist(lapply(dl,"length")))
pct <- 0
a <- "0"
for(i in seq(length(dl[[1]]))) for(j in seq(length(dl[[2]])))\{
  pct <- pct +1
  if(paste(round(pct/((length(dl[[1]])*length(dl[[2]])))*100,digit=0), "\%") != a)\{
    a <- print(paste(round(pct/((length(dl[[1]])*length(dl[[2]])))*100,digit=0), "\%"))\}
  po<-probest(sotiriou.data,
				c("AGE","TUMOR_SIZE"),
				"recur",
				"TAMOXIFEN_TREATMENT",
				"ELSTON.ELLIS_GRADE",
				c(dl[[1]][i],dl[[2]][j]))
  statv[i,j]<-(t(po$b-.5)\%*\%solve(po$Vb)\%*\%(po$b-.5))[1,1]
\}

##Plot for the confidence regions at 90\%, 95\%and 97.5\%.
contour(dl[[1]],dl[[2]],statv,levels=qchisq(c(.9,.95,.975),2),
        xlab="Age (years)",ylab="Tumor Size",
        main="Confidence Region offsets of Age and Tumor Size",
        sub="Breast Cancer Recurrence Data", col=c("red","blue","darkolivegreen"))
}
