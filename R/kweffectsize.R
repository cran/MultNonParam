#' @title Sample Size for the Kruskal-Wallis test.
#' @description \code{kweffectsize} approximates effect size for the Kruskal-Wallis test, 
#' using a chi-square approximation under the null, and a non-central chi-square approximation under the alternative.  The noncentrality parameter is calculated using alternative means and the null variance structure.
#' @param totsamp sample size
#' @param shifts The offsets for the various populations, under the alternative hypothesis.  This is used for direction on input.
#' @param distname The distribution of the underlying observations; normal and logistic are currently supported.
#' @param targetpower The distribution of the underlying observations; normal and logistic are currently supported.
#' @param proportions The proportions in each group.
#' @param level The test level.
#' @details The standard noncentral chi-square power formula, or Monte Carlo, is used.
#' @return A list with components power, giving the power approximation, ncp, giving the noncentrality parameter, cv, giving the critical value, probs, giving the intermediate output from pairwiseprobability, and expect, the quantities summed before squaring in the noncentrality parameter.
#' @importFrom stats qchisq
#' @importFrom stats pchisq
#' @importFrom stats rlogis
#' @export
#' @examples
#' kwpower(rep(10,3),c(0,1,2),"normal")
kweffectsize<-function(totsamp,shifts,distname=c("normal","logistic"),
   targetpower=0.8,proportions=rep(1,length(shifts))/length(shifts),level=0.05){
   kappacirc<-probabilityderiv(distname)
   myncp<-solvencp(length(proportions)-1,level=level,targetpower=targetpower)
   denominator<-kappacirc*sqrt(totsamp*12*(sum(shifts^2*proportions)-sum(shifts*proportions)^2))
   return(sqrt(myncp)/denominator)
}
