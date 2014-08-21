\name{prostate.data}
\docType{data}
\alias{prostate.data}
\title{Prostate cancer data set}
\description{
   221 prostate cancer patients are collected in this data set.
}
\usage{prostate.data}
\format{
hosp : Hospital in which the patient is hospitalized.\cr
stage : stage of the cancer.\cr
gleason score : used to help evaluate the prognosis of the cancer.\cr
psa : prostate-specific antigen.\cr
age : age of the patient.\cr
advanced : boolean. \code{TRUE} if the cancer is advanced.\cr
}
\source{A. V. D'Amico, R. Whittington, S. B. Malkowicz, D. Schultz, K. Blank, G. A. Broderick, J. E. Tomaszewski, A. A. Renshaw, I. Kaplan, C. J. Beard, A. Wein (1998) , \emph{Biochemical outcome after radical prostatectomy, external beam radiation therapy, or interstitial radiation therapy for clinically localized prostate cancer}, JAMA : the journal of the American Medical Association 280 969-74.}
\examples{
data(prostate.data)
attach(prostate.data)
	plot(age,psa,
     main="Age and PSA",
     sub="Prostate Cancer Data",
     xlab="Age (years)",ylab="PSA")
}

