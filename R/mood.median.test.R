#' @title Mood's Median test, extended to odd sample sizes.
#' @description Test whether two samples come from the same distribution.  This version of Mood's median test is presented for pedagogical purposes only.  Many authors successfully argue that it is not very powerful.  The name "median test" is a misnomer, in that the null hypothesis is equality of distributions, and not just equality of median.  Exact calculations are not optimal for the odd sample size case.
#' @param x First data set.
#' @param y Second data set.
#' @param exact Indicator for whether the test should be done exactly or approximately.
#' @details The exact case reduces to Fisher's exact test.
#' @return The two-sided p-value.
#' @export
#' @importFrom stats pnorm
#' @importFrom stats median
#' @importFrom stats fisher.test
mood.median.test<-function(x,y,exact=FALSE){
   mm<-median(c(x,y))
   a21<-sum(x>mm)
   a22<-sum(x<mm)
   a11<-sum(y>mm)
   a12<-sum(y<mm)
   if(exact){
      p.value<-fisher.test(cbind(c(a11,a12),c(a21,a22)))$p.value
   }else{
      rowtot<-c(a11+a12,a21+a22)
      coltot<-c(a11+a21,a12+a22)
      nn<-sum(rowtot)
      v<-prod(rowtot)*prod(coltot)/(nn^2*(nn-1))
      if(a11>rowtot[1]*coltot[1]/nn){
         a11<-a11-.5
      }else{
         a11<-a11+.5
      }
      z<-(a11-rowtot[1]*coltot[1]/nn)/sqrt(v)
      p.value<-2*pnorm(-abs(z))
   }
   return(list(p.value=p.value))
}
