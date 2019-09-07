#' @title Power for the Kruskal-Wallis test.
#' @description \code{kwpower} approximates power for the Kruskal-Wallis test, 
#' using a chi-square approximation under the null, and a non-central chi-square approximation under the alternative.  The noncentrality parameter is calculated using alternative means and the null variance structure.
#' @param nreps The numbers in each group.
#' @param shifts The offsets for the various populations, under the alternative hypothesis.
#' @param distname The distribution of the underlying observations; normal and logistic are currently supported.
#' @param level The test level.
#' @param mc 0 for asymptotic calculation, or positive for mc approximation.
#' @param taylor logical determining whether Taylor series approximation is used for probabilities.
#' @details The standard noncentral chi-square power formula, or Monte Carlo, is used.
#' @return A list with components power, giving the power approximation, ncp, giving the noncentrality parameter, cv, giving the critical value, probs, giving the intermediate output from pairwiseprobability, and expect, the quantities summed before squaring in the noncentrality parameter.
#' @importFrom stats qchisq
#' @importFrom stats pchisq
#' @importFrom stats rlogis
#' @export
#' @examples
#' kwpower(rep(10,3),c(0,1,2),"normal")
kwpower<-function(nreps,shifts,distname=c("normal","logistic"),level=0.05,mc=0,taylor=FALSE){
   vartheta<-pairwiseprobabilities(shifts,distname,taylor=taylor)
   expect<-(vartheta-.5)%*%nreps
   ncp<-12*sum(nreps*expect^2)/(sum(nreps)*(sum(nreps)+1))
   cv<-qchisq(1-level,length(nreps)-1)
   if(mc==0) power<-1-pchisq(cv,length(nreps)-1,ncp)
   if(mc>0){
      count<-0
      samptot<-sum(nreps)
      g<-rep(seq(length(nreps)),nreps)
      for(j in seq(mc)){
         x<-switch(distname,normal=rnorm(samptot),logistic=rlogis(samptot))+
            shifts[g]
         count<-count+(kruskal.test(x,g,alternative="greater")$p.value<level)
      }
      power<-count/mc
   }
   return(list(power=power,ncp=ncp,cv=cv,probs=vartheta,expect))
}
