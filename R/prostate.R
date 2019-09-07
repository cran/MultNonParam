#' @title prostate
#' @name prostate
#' @description 221 prostate cancer patients are collected in this data set.
#' @docType data
#' @format
#' \itemize{
#' \item hosp : Hospital in which the patient is hospitalized.
#' \item stage : stage of the cancer.
#' \item gleason score : used to help evaluate the prognosis of the cancer.
#' \item psa : prostate-specific antigen.
#' \item age : age of the patient.
#' \item advanced : boolean. \code{TRUE} if the cancer is advanced.
#'}
#' @references
#' A. V. D'Amico, R. Whittington, S. B. Malkowicz, D. Schultz, K. Blank, G. A. Broderick, J. E. Tomaszewski, A. A. Renshaw, I. Kaplan, C. J. Beard, A. Wein (1998) , \emph{Biochemical outcome after radical prostatectomy, external beam radiation therapy, or interstitial radiation therapy for clinically localized prostate cancer}, JAMA : the journal of the American Medical Association 280 969-74.
#'
#' @examples
#' data(prostate)
#' attach(prostate)
#' plot(age,psa,main="Age and PSA",sub="Prostate Cancer Data",
#'    xlab="Age (years)",ylab="PSA")
NULL
