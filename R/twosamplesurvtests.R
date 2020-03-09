#' @title  Two Sample Omnibus Tests of Survival Curves
#' @description Returns the Kolmogorov-Smirnov and Anderson-Darling test statistics for two right-censored data sets.
#' @param times Event and censoring times
#' @param delta Indicator of event (1) or censoring (0).
#' @param grp Variable that divides the population into groups.
#' @details The function calls a Fortran code to calculate the estimators \code{b} and their variance-covariance matrix \code{Vb}
#' @return  A vector of length two, with the Kolmogorov-Smirnov and Anderson-Darling statistics.
#' @export
#' @useDynLib MultNonParam tskmsurv
twosamplesurvtests<-function(times,delta,grp){
  outa<-.Fortran("tskmsurv",
     as.integer(length(times)),
     as.integer(rank(times,ties.method="min")), 
     as.integer(delta), 
     as.integer(length(table(grp))),
     as.integer(as.factor(grp)),
     as.integer(0),
     stats=as.double(0),
#    names=as.character("   "),
     PACKAGE="MultNonParam")
  out<-.Fortran("tskmsurv",
     as.integer(length(times)),
     as.integer(rank(times,ties.method="min")), 
     as.integer(delta), 
     as.integer(length(table(grp))),
     as.integer(as.factor(grp)),
     as.integer(outa[[6]]),
     stats=as.double(rep(0,outa[[6]])),
#    names=as.character(rep("   ",outa[[6]])),
     PACKAGE="MultNonParam")
  stats<-out$stats
# names(stats)<-c("KolmogorovSmirnov","AndersonDarling")
# names(stats)<-out$names
  names(stats)<-c("KS ","AD ","CM1","CM2")
  return(stats)
}
