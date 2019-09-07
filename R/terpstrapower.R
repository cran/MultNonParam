#' @title Power for the nonparametric Terpstra test for an ordered effect.
#' @description \code{terpstrapower} approximates power for the one-sided Terpstra test, 
#' using a normal approximation with expectations under the null and alternative, and using the null standard deviation.
#' @param nreps The numbers in each group.
#' @param shifts The offsets for the various populations, under the alternative hypothesis.
#' @param distname The distribution of the underlying observations; normal and logistic are currently supported.
#' @param level The test level.
#' @param mc Zero indicates asymptotic calculation.  Positive for MC calculation.
#' @details The standard normal-theory power formula is used.
#' @return A list with components power, giving the power approximation, expect, giving null and alternative expectations, var, giving the null variance, probs, giving the intermediate output from pairwiseprobability, and level.
#' @export
#' @importFrom stats qnorm
#' @examples
#' terpstrapower(rep(10,3),c(0,1,2),"normal")
#' terpstrapower(c(10,10,10),0:2,"normal",mc=1000)
terpstrapower<-function(nreps,shifts,distname=c("normal","logistic"),level=0.025,mc=0){
   vartheta<-pairwiseprobabilities(shifts,distname)
   totsamp<-sum(nreps)
   expect<-c(0,0)
   for(j in seq(length(nreps)-1)){
      for(k in (j+1):length(nreps)){
         expect<-expect+nreps[j]*nreps[k]*c(.5,vartheta[j,k])
      }
   }
   var<-totsamp*(totsamp+1)*(2*totsamp+1)
   for(k in seq(length(nreps))) var<-var-nreps[k]*(nreps[k]+1)*(2*nreps[k]+1)
   var<-var/72
   cv<-qnorm(1-level);#message(cv); message(-diff(expect)/sqrt(var))
   if(mc==0){
      power<-1-pnorm(-diff(expect)/sqrt(var)+cv)
      names(power)<-"terpstra"
   }
   if(mc>0){
      count<-0
      samptot<-sum(nreps)
      g<-rep(seq(length(nreps)),nreps)
      for(j in seq(mc)){
         x<-switch(distname,normal=rnorm(samptot),logistic=rlogis(samptot))+
            shifts[g]
         count<-count+(terpstra.test(x,g,alternative="greater")$p.value<level)
      }
      power<-count/mc
   }
   return(list(power=power,expect=expect,var=var,probs=vartheta,level=level))
}
