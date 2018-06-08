#' @title Compare the sensitivity of different statistics.
#' @description Compare the sensitivity of different statistics.
#' @param y vector of the data.
#' @param sub subtitle for the plot.
#' @param stats vector of functions to be plotted.
#' @details To compare the sensitivity, outliers are added to the original data. The shift of each statistics due to the new value is measured and plotted.
#' @export
#' @importFrom graphics plot
#' @importFrom graphics legend
sensitivity.plot<-function(y,sub,stats){
  ra<-range(y)
  xr<-mean(ra)+c(-1,1)*diff(ra)
  outlier<-xr[1]+(0:100)*diff(xr)/100
  base<-stats(y)
  sens<-array(NA,c(length(base),length(outlier)))
  dimnames(sens)[[1]]<-names(base)
  for(j in seq(length(outlier))){
    sens[,j]<-stats(c(y,outlier[j]))-base
  }
  plot(xr,range(sens),type="n",main="Sensitivity",sub=sub,
       ylab="Change in Statistic Value",xlab="New Observation")
  inds<-seq(length(base))
  for(i in inds) lines(outlier,sens[i,],lty=i,col=i)
  legend(min(ra),max(sens),lty=inds,col=inds,
         legend=names(base))
}
