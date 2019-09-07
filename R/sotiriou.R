#' @title Breast cancer data set
#' @name sotiriou
#' @description 187 breast cancer patients are collected in this data set.
#' @docType data
#' @usage data(sotiriou)
#' @format A data set with the following variables
#' \itemize{
#' \item AGE : Age of the patient
#' \item TUMOR_SIZE : The size of the tumor, numeric variable
#' \item recur : 1 if the patient has a recurent breast cancer, 0 if it is not reccurent.
#' \item ELSTON.ELLIS_GRADE : Elston Ellis grading system in order toclassify the breast cancers. It can be a low, intermediate or high grade (high being the worst prognosis)
#' \item TAMOXIFEN_TREATMENT : boolean. \code{TRUE} if the patient is treated with the Tamoxifen treatment.
#'}
#' @references
#' S. Madhavan, Y. Gusev, M. Harris, D. Tanenbaum, R. Gauba, K. Bhuvaneshwar, A. Shinohara, K. Rosso, L. Carabet, L. Song, R. Riggins, S. Dakshanamurthy, Y. Wang, S. Byers, R. Clarke, L. Weiner (2011), \emph{A systems medicine platform for personalized oncology}, Neoplasia 13. 
#'
#' C. Sotiriou, P. Wirapati, S. Loi, A. Harris, S. Fox, J. Smeds, H. Nordgren, P. Farmer, V. Praz, B. Haibe-Kains, C. Desmedt, D. Larsimont, F. Cardoso, H. Peterse, D. Nuyten, M. Buyse, M. Van de Vijver, J. Bergh, M. Piccart, M. Delorenzi  (2006), \emph{Gene expression profiling in breast cancer: understanding the molecular basis of histologic grade to improve prognosis}, Journal of the National Cancer Institute 98 262-72.
#' @examples 
#' data(sotiriou)
#' plot(sotiriou$AGE,sotiriou$TUMOR_SIZE,pch=(sotiriou$recur+1),
#'    main="Age and Tumor Size",
#'    sub="Breast Cancer Recurrence Data",
#'    xlab="Age (years)",ylab="Tumor Size",
#'    col=c("blue","darkolivegreen"))
#' legend(31,8,legend=c("Not Recurrent","Recurrent"),pch=1:2,
#'    col=c("blue","darkolivegreen"))
#' 
#' @source https://gdoc.georgetown.edu/gdoc/
"sotiriou"
