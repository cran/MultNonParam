#' Calculate the quantiles of the count of concordant pairs among indpendent pairs of random variables.
#' @param qq Desired quantile
#' @param nn number of pairs
#' @param exact flag to trigger exact calculation when possible.
#' @return Integer quantile
#' @export
qconcordant<-function(qq,nn,exact=TRUE){
   if((nn<=20)&exact){
      if((qq<0)|(qq>1)) dc<-NA
      if(qq==0) dc<-0
      if(qq==1) dc<-nn*(nn-1)/2
      if((qq>0)&(qq<1)) dc<-.Fortran("qconcordant", as.double(qq), as.integer(nn), dc=as.integer(0), PACKAGE="MultNonParam")$dc
   }else{
      if((qq<0)|(qq>1)) dc<-NA
      if(qq==0) dc<-0
      if(qq==1) dc<-nn*(nn-1)/2
      v<-nn*(nn-1)*(5+2*nn)/72
      if((qq>0)&(qq<1)) dc<-round(nn*(nn-1)/4+ qnorm(qq)*sqrt(v))
   }
   return(dc)
}
