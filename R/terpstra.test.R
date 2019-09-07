terpstra.test<-function(x,g,alternative=c("two.sided","less","greater")){
#' Perform the Terpstra version of the multi-ordered-sample test
#' @param x A vector of values from all samples.
#' @param g A vector of group labels.
#' @param alternative Specification of alternative hypothesis.
#' @return Test results of class htest
#' @examples
#' terpstra.test(rnorm(15),rep(1:3,5))
#' @export
#' @importFrom stats pnorm
   alternative<-match.arg(alternative)
   STATISTIC<-0
   expect<-0
   var<-length(x)*(length(x)+1)*(2*length(x)+1)
   gg<-unique(sort(g))
   nreps<-rep(NA,length(gg))
   for(kk in seq(length(gg))){
      nreps[kk]<-sum(g==gg[kk])
      var<-var-nreps[kk]*(nreps[kk]+1)*(2*nreps[kk]+1)
   }
   for(kk in seq(length(gg)-1)){
      for(ll in (kk+1):length(gg)){
         STATISTIC<-STATISTIC+mannwhitney.test(x[g==gg[kk]],x[g==gg[ll]])$statistic
         expect<-expect+nreps[kk]*nreps[ll]/2
      }
   }
   var<-var/72
   names(STATISTIC)<-"terpstra"
   z<-(STATISTIC-expect)/sqrt(var)
   PVAL<-switch(alternative,two.sided=2*pnorm(-abs(z)),greater=pnorm(-z),less=pnorm(z))
      testout<-list(
      statistic=STATISTIC,
      p.value=PVAL,
      alternative=alternative,
      method="Terpstra",
#     null.value=list(shift=0),
#     estimate=NA,
      z=z,
      dname=paste(deparse(substitute(x)),"classified by",deparse(substitute(g))))
   class(testout)<-"htest"
   return(testout)
}
