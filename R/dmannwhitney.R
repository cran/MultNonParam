#' @title Mann Whitney Probability Mass function
#' @description Calculates the Mann Whitney Probability Mass function recursively.
#' @param u Statistic value
#' @param m Group 1 size
#' @param n Group 2 size
#' @return Probability that the Mann-Whitney statistic takes the value u under H0
#' @export
dmannwhitney<-function(u,m,n){
  if((u<0)|(u>(m*n))){
    val<-0
  }else{
    if(m==1) val<-1
    if(n==1) val<-1
    if((m>1)&(n>1)) val<-(m*dmannwhitney(u,m-1,n)+n*dmannwhitney(u-m,m,n-1))/(m+n)
  }
  return(val)
}

