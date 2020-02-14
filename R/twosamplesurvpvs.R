#' @title  Two Sample Omnibus Tests of Survival Curves
#' @description Returns the Kolmogorov-Smirnov and Anderson-Darling test statistics for two right-censored data sets.
#' @param times Event and censoring times
#' @param delta Indicator of event (1) or censoring (0).
#' @param grp Variable that divides the population into groups.
#' @param nmc Number of Monte Carlo samples for p value calculation
#' @param plotme logical; indicates whether to plot or not.
#' @param exact logical; indicates whether to use exhaustive enumeration of permutations or not.
#' @details The function calls a Fortran code to calculate the estimators \code{b} and their variance-covariance matrix \code{Vb}
#' @return  A vector of length two, with the Kolmogorov-Smirnov and Anderson-Darling statistics.
#' @export
#' @examples
#' twosamplesurvpvs(rexp(20),rbinom(20,1,.5),rbinom(20,1,.5))
#' @useDynLib MultNonParam tskmsurv
twosamplesurvpvs<-function(times,delta,grp,nmc=10000,plotme=TRUE,exact=FALSE){
   tobs<-twosamplesurvtests(times,delta,grp)
   if(exact){
      pvalsa<-.Fortran("tskmsurvpv",     as.integer(length(times)),
         as.integer(rank(times,ties.method="min")), as.integer(delta),
         as.integer(length(table(grp))), as.integer(as.factor(grp)),
         as.integer(0),pvals=as.double(c(0,0)),as.character("   "),
         as.integer(0), PACKAGE="MultNonParam")$pvals
      out<-.Fortran("tskmsurvpv",     as.integer(length(times)),
         as.integer(rank(times,ties.method="min")), as.integer(delta),
         as.integer(length(table(grp))), as.integer(as.factor(grp)),
         as.integer(pvalsa[[6]]),pvals=as.double(rep(0,pvalsa[[6]])), 
         as.character(rep("   ",pvalsa[[6]])),as.integer(0),
         PACKAGE="MultNonParam")
      pvals<-out$pvals
#     names(pvals)<-out[[8]]
      names(pvals)<-c("KS ","AD ","CM1","CM2")
   }else{
      count<-rep(0,length(tobs))
      out<-array(NA,c(nmc,length(tobs)))
      for(jj in seq(nmc)){
         onetests<-twosamplesurvtests(times,delta,sample(grp))
         if(jj==1){
            dimnames(out)<-list(NULL,names(onetests))
            names(count)<-names(onetests)
         }
         count<-count+((out[jj,]<-onetests)>=tobs)
      }
      if(plotme){
         plot(c(0,1),range(out),type="n")
         for(jj in seq(length(tobs))){
            lines(seq(nmc)/nmc,sort(out[,jj]),lty=jj)
            abline(h=tobs[jj],lty=jj)
         }
         legend(0,max(out),lty=seq(length(tobs)),legend=names(tobs))
      }
      names(count)<-names(tobs)
      pvals<-count/nmc
   }
   return(list(pvals=pvals,stats=tobs))
}
