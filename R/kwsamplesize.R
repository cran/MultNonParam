#' @title Sample Size for the Kruskal-Wallis test.
#' @description \code{kwsamplesize} approximates sample size for the Kruskal-Wallis test, 
#' using a chi-square approximation under the null, and a non-central chi-square approximation under the alternative.  The noncentrality parameter is calculated using alternative means and the null variance structure.
#' @param shifts The offsets for the various populations, under the alternative hypothesis.
#' @param distname The distribution of the underlying observations; normal and logistic are currently supported.
#' @param targetpower The distribution of the underlying observations; normal and logistic are currently supported.
#' @param proportions The proportions in each group.
#' @param level The test level.
#' @param taylor Logical flag forcing the approximation of exceedence probabilities using the first derivative at zero.
#' @details The standard noncentral chi-square power formula, is used.
#' @return A list with the total number of observations needed to obtain approximate power, as long as this number is split amomg groups according to argument proportion.
#' @importFrom stats qchisq
#' @importFrom stats pchisq
#' @importFrom stats rlogis
#' @export
#' @examples
#' #Calculate the sample size necessary to detect differences among three
#' #groups with centers at 0,1,2, from normal observations, using a test of
#' #level 0.05 and power 0.80.
#' kwsamplesize(c(0,1,2),"normal")
kwsamplesize<-function(shifts,distname=c("normal","logistic","cauchy"),targetpower=0.8,proportions=rep(1,length(shifts))/length(shifts),level=0.05,taylor=FALSE){
   vartheta<-pairwiseprobabilities(shifts,distname,taylor)
   expect<-(vartheta-.5)%*%proportions
   myncp<-solvencp(length(proportions)-1,level=level,targetpower=targetpower)
   denominator<-12*sum(proportions*expect^2)
   return(myncp/denominator)
}
