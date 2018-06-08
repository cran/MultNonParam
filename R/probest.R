#' @title  Stratified Multivariate Kawaguchi Koch Wang  Estimators
#' @description Function that return the estimators and their variance-covariance matrix calculated with the Kawaguchi - Koch - Wang method.
#' @param ds The data frame to be used.
#' @param resp The vector of the response manifest variable. There can be more than one variable. It has to be the name of the variable as a character string.
#' @param grp The vector of the variable that divides the population into groups. It has to be the name of the variable as a character string.
#' @param str The vector of the variable used for the strata. It has to be the name of the variable as a character string.
#' @param covs The covariates to be used in the model. It has to be the name of the variable as a character string.
#' @param delta Offeset for covariates.
#' @param correct Should the variance estimator be corrected as in Chen and Kolassa?
#' @details The function calls a Fortran code to calculate the estimators \code{b} and their variance-covariance matrix \code{Vb}
#' @return  A list with components b, the vector of adjusted estimates from the method, and Vb, the corresponding estimated covariance matrix.
#' @references
#' A. Kawaguchi, G. G. Koch and X. Wang (2012), "Stratified Multivariate Mann-Whitney Estimators for the Comparison of Two Treatments with Randomization Based Covariance Adjustment", \emph{Statistics in Biopharmaceutical Research} 3 (2) 217-231.
#'
#'J. E. Kolassa and Y. Seifu (2013), Nonparametric Multivariate Inference on Shift Parameters, \emph{Academic Radiology} 20 (7), 883-888.
#'
#' @examples
#'# Breast cancer data from the MultNonParam package.
#' data(sotiriou)
#' attach(sotiriou)
#'#First simple plot of the data
#'plot(AGE,TUMOR_SIZE,pch=(recur+1),main="Age and Tumor Size",
#'   sub="Breast Cancer Recurrence Data",xlab="Age (years)",
#'   ylab="Tumor Size",col=c("blue","darkolivegreen"))
#'legend(31,8,legend=c("Not Recurrent","Recurrent"),
#'   pch=1:2,col=c("blue","darkolivegreen"))
#'#AGE and TUMOR_SIZE are the response variables, recur is used for the groups,
#'#TAMOXIFEN_TREATMENT for the stratum and ELSTON.ELLIS_GRADE is a covariate.
#'po<-probest(sotiriou,c("AGE","TUMOR_SIZE"),"recur",
#'   "TAMOXIFEN_TREATMENT","ELSTON.ELLIS_GRADE")
#' @export
#' @importFrom graphics legend
#' @importFrom graphics plot
#' @importFrom stats model.matrix
#' @useDynLib MultNonParam probestf
probest<-function(ds,resp,grp,str=NULL,covs=NULL,delta=NA,correct=FALSE){
  r<-length(resp)
  if(length(delta)==1){if(is.na(delta)) delta<-rep(0,r)}
  if(length(delta)==0) delta<-rep(0,r)
  covariates<-if(length(covs)<2){
    if(is.null(covs)){
      as.double(NULL)
    }else{
      if(is.factor(ds[,covs])){
        model.matrix(~ds[,covs]-1)[,-1,drop=F]
      }else{
        ds[,covs,drop=F]
      }
    }
  }else{
    ds[,covs]
  }
  M<-as.integer(if(is.null(covs)) 0 else dim(covariates)[2])
  N<-dim(ds[,resp,drop=F])[1]
  grpv<-as.numeric(as.factor(ds[[grp]]))
  gn<-sort(unique(grpv))
  if(length(gn)!=2) cat("I am surprised that there are more than two groups.  Correction will not work.\n")
  if(is.null(str)){
     strv<-rep(1,N)
  }else{
     strv<-as.numeric(as.factor(ds[[str]]))
  }
  ustr<-unique(strv)
  out<-.Fortran(probestf,
                as.integer(length(resp)), M, as.integer(N), as.integer(grpv),
                as.integer(length(gn)), as.integer(gn),
                as.integer(strv), as.integer(ustr),as.integer(length(ustr)),
                as.double(as.matrix(ds[,resp])), as.double(as.matrix(covariates)),
                as.logical(!is.na(ds[,resp])),as.double(delta),
                b=as.double(rep(0,r)),Vb=as.double(rep(0,r^2)),
                correct=as.logical(correct))
  return(list(b=out$b,Vb=array(out$Vb,c(r,r))))
}
