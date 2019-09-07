#' @title Permutation test of assication
#' @description Calculate the p-value for the test of association between two variables using the permutation method.
#' @param x First vector to be associated.
#' @param y First vector to be associated.
#' @return p-value
#' @examples
#' #Example using data from plant Qn1 from the CO2 data set.^M
#' betatest(CO2[CO2$Plant=="Qn1",4],CO2[CO2$Plant=="Qn1",5])
#' @export
#' @useDynLib MultNonParam betatestf
betatest<-function(x,y){
  if(length(x)>10) cat("Warning: this will take forever.")
  out<-.Fortran("betatestf",as.integer(length(x))
                ,as.double(x)
                ,as.double(y)
                ,pval=as.double(0),PACKAGE="MultNonParam")
  return(out$pval)
}
