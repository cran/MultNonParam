theil<-function(x,y,conf=.9){
#' Perform the Theil nonparametric estimation and confidence interval for a slope parameter.
#' @param x A vector of values of the explanatory variable.
#' @param y A vector of values of the response variable.
#' @param conf Level of confidence interval.
#' @return A list with letters and numbers.
#' \itemize{
#'   \item est - An estimate, the median of pairwise slopes.
#'   \item ci - A vector of confidence interval endpoints.
#' }
#' @examples
#' a<-0:19;b<-a^2.5
#' theil(a,b)
#' @export
#' @importFrom stats median
   n<-length(x)
   out<-rep(NA,n*(n-1)/2)
   k<-0
   for(i in seq(n-1)) for(j in i:n){
      k<-k+1;out[k]<-(y[j]-y[i])/(x[j]-x[i])
   }
   outs<-sort(out)
#  message(c(length(outs),length(unique(outs))))
   est<-median(outs)
   b<-qconcordant((1-conf)/2,n)
#  message(b)
   ci<-outs[c(b,n*(n-1)/2+1-b)]
   return(list(est=est,ci=ci))
}
