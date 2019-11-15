#' Calculate the cumulative distribution of the count of concordant pairs among indpendent pairs of random variables.
#' @param ss Integer number of pairs
#' @param nn number of pairs
#' @return real probability
#' @export
pconcordant<-function(ss,nn){
   if(nn<=20){
      if(ss<0) dc<-0
      if(ss>(nn*(nn-1)/2)) dc<-1
      if((ss>=0)&(ss<=(nn*(nn-1)/2))) dc<-.Fortran("pconcordant", as.integer(ss), as.integer(nn), dc=as.double(0.0), PACKAGE="MultNonParam")$dc
   }else{
     dc<-NA
   }
   return(dc)
}
