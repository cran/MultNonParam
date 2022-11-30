#' @title Nonparametric Confidence Region for a Vector Shift Parameter
#' @description Inversion of a one-sample bivariate rank test is used to produce a confidence region.  The region is constructed by building a grid of potential parameter values, evaluating the test statistic on each grid point, collecting the p-values, and then drawing the appropriate countour of the p-values.  The grid is centered at the bivariate median of the data set.
#' @param xm A two-column matrix of bivariate data whose two location parameters are to be estimated.
#  @param hpts The number of grid points on either side of the median.
#' @return nothing
#' @export
#' @importFrom stats median IQR
#' @importFrom graphics contour
#' @importFrom ICSNP rank.ctest
shiftcr<-function(xm,hpts=50){
   bmo<-rep(1,dim(xm)[1])
   th1<-2*IQR(xm[,1])*(-hpts:hpts)/(2*hpts)+median(xm[,1])
   th2<-2*IQR(xm[,2])*(-hpts:hpts)/(2*hpts)+median(xm[,2])
   pvo<-array(NA,c(length(th1),length(th2)))
   for(i in seq(length(th1))) for(j in seq(length(th2)))
      pvo[i,j]<-rank.ctest(
         xm-outer(rep(1,dim(xm)[1]),c(th1[i],th2[j])))$p.value
   contour(th1,th2,pvo,levels=.05,
      main="Median Blood Pressure Change Confidence Region",
      xlab="Diastolic",ylab="Systolic")
}
#shiftcr(bp[,c("dpd","spd")])
