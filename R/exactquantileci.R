#' @title Exact Quantile Confidence Interval
#' @description Calculates exact quanitle confidence intervals by inverting the generalization of the sign test.
#' @param xvec vector of observations
#' @param tau quantile to be estimated.  If this is a vector, separate intervals and tests for each value will be calculated.
#' @param alpha 1-confidence level.
#' @param md null value of quantile
#' @return A list with components cis, an array with two columns, representing lower and upper bounds, and a vector pvals, of p-values.
#' @export
#' @importFrom stats qbinom
#' @importFrom stats pbinom
exactquantileci<-function(xvec,tau=.5,alpha=0.05,md=0){
   n<-length(xvec)
   xvec<-sort(xvec)
   cis<-array(NA,c(2,length(tau)))
   pvals<-rep(NA,length(tau))
   tt<-sum(xvec<=md)
   for(j in seq(length(tau))){
      ii<-qbinom(alpha/2,n,tau[j])
      cis[1,j]<-if(ii>=1) xvec[ii] else -Inf
      ii<-n+1-qbinom(alpha/2,n,1-tau[j])
      cis[2,j]<-if(ii<=n) xvec[ii] else Inf
      pvals[j]<-2*min(c(pbinom(tt,n,tau), pbinom(n+1-tt,n,1-tau)))
   }
   return(list(cis=cis,pvals=pvals))
}
