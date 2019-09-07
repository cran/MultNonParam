#' @title Diagnosis for multivariate stratified Kawaguchi - Koch - Wang method
#' @description Diagnostic tool that verifies the normality of the estimates of the probabilities b with the Kawaguchi - Koch - Wang method. The diagnostic method is based on a Monte Carlo method.
#' @param n number of observations in the first group.
#' @param m number of observations in the second group.
#' @param k number of strata.
#' @param nsamp The number of estimates that will be calculated. Must be enough to be sure that the results are interpretable.
#' @param delta Offset that depends on group.
#' @param beta Correlation between x and y.
#' @param disc The Mann Whitney test is designed to handle continuous data, but this method applies to discretized data; \code{disc} adjusts the discreteness.
#' @details This functions serves as a diagnosis to prove that the Kawaguchi - Koch - Wang method gives Gaussian estimates for b. It generates random data sets, to which the Mann Whitney test gets applied.  \code{y} is the generated response variable and \code{x} the generated covariable related to \code{y} through a regression model.
#' @return Nothing is returned.  A QQ plot is drawn.
#' @references
#' A. Kawaguchi, G. G. Koch and X. Wang (2012), "Stratified Multivariate Mann-Whitney Estimators for the Comparison of Two Treatments with Randomization Based Covariance Adjustment", \emph{Statistics in Biopharmaceutical Research} 3 (2) 217-231.
#'
#'J. E. Kolassa and Y. Seifu (2013), Nonparametric Multivariate Inference on Shift Parameters, \emph{Academic Radiology} 20 (7), 883-888.
#'
#' @examples
#' testve(10,15,3,100,0.4)
#' @export
#' @importFrom graphics abline
#' @importFrom stats rnorm
#' @importFrom stats qqnorm
testve<-function(n,m,k,nsamp=100,delta=0,beta=0,disc=0){
  o<-rep(NA,nsamp)
  i<-rep(c(rep(0,m),rep(1,n)),k)
  str<-rep(1:k,rep(n+m,k))
  for(j in seq(nsamp)){
    x<-rnorm(k*(n+m))
    y<-rnorm(k*(n+m))
    if(disc!=0) y<-round(y/disc)*disc
    ds<-as.data.frame(list(str=str,x=x,y=y+(i-1)*delta+beta*x, i=i))
    out <-probest(ds,"y","i","str","x",0)
    o[j]<-(out$b-.5)/sqrt(out$Vb)
  }
  st<-paste("Delta=",delta,"beta=",beta,"m=",m,"n=",n)
  if(disc!=0) st<-paste(st,"made discrete")
  aaa<-qqnorm(o,plot.it=F)
  plot(aaa$x,aaa$y,main="Normal QQ Plot for Distribution of KKW Statistic",
       sub=st,type="p")
  abline(b=1,a=0)
  #  qqline(o)
}
