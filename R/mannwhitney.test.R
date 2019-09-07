mannwhitney.test<-function(x,y,alternative=c("two.sided","less","greater")){
#' Perform the Mann Whitney two-sample test
#' @param x A vector of values from the first sample.
#' @param y A vector of values from the first sample.
#' @param alternative Specification of alternative hypothesis.
#' @return Test results of class htest
#' @examples
#' mannwhitney.test(rnorm(10),rnorm(10)+.5)
#' @export
#' @importFrom stats pnorm
   alternative<-match.arg(alternative)
   lx<-length(x);ly<-length(y)
   STATISTIC<-sum(outer(y,x,"-")>0)
   names(STATISTIC)="Mann-Whitney"
   z<-(STATISTIC-lx*ly/2)/sqrt(lx*ly*(lx+ly+1)/12)
   PVAL<-switch(alternative,
      two.sided=2*pnorm(-abs(z)),greater=pnorm(-z),less=pnorm(z))
   testout<-list(
      statistic=STATISTIC,
      p.value=PVAL,
      alternative=alternative,
      method="Mann-Whitney",
#     null.value=list(shift=0),
#     estimate=NA,
      dname=paste(deparse(substitute(x)),"and",deparse(substitute(y))))
   class(testout)<-"htest"
   return(testout)
}
