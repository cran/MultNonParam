#' @title Nonparametric Confidence Region for a Vector Shift Parameter
#' @description Inversion of a one-sample bivariate rank test is used to produce a confidence region.
#' @param xm A bivariate data vector whosse two location parameters are to be estimated.
#' @return nothing
#' @export
#' @importFrom stats median sd
#' @importFrom graphics contour
#' @importFrom ICSNP rank.ctest
shiftcr<-function(xm){
   bmo<-rep(1,dim(xm)[1])
   xv<-2.5*sd(xm[,1])*(-50:50)/100+median(xm[,1])
   yv<-2.5*sd(xm[,2])*(-50:50)/100+median(xm[,2])
   pvo<-array(NA,c(length(xv),length(yv)))
   for(i in seq(length(xv))) for(j in seq(length(yv)))
      pvo[i,j]<-rank.ctest(
         xm-outer(rep(1,dim(xm)[1]),c(xv[i],yv[j])))$p.value
   contour(xv,yv,pvo,levels=.05,
      main="Median Blood Pressure Change Confidence Region",
      xlab="Diastolic",ylab="Systolic")
}
#shiftcr(bp[,c("dpd","spd")])
