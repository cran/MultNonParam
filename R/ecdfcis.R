#' @title Confidence Intervals for Empirical Cumulative Distribution Functions
#' @param data vector of observations
#' @param alpha 1-confidence level.
#' @param dataname Name of variable for use in axis labeling
#' @param exact logical value controlling whether confidence intervals are exact or asymlptotic.
#' @param newplot logical value controlling whether the estimate is added to an existing plot, or whether a new plot should be constructed.
#' @export
#' @importFrom graphics segments
#' @importFrom stats binom.test quantile
ecdfcis<-function(data,alpha=0.05,dataname=NA,exact=TRUE,newplot=TRUE){
   startlty<-if(newplot) 0 else 2
   x<-sort(unique(data))
   x<-c(min(x)-1,x,max(x)+1)
   nn<-length(x)
   y<-rep(NA,nn)
   for(j in seq(nn)) y[j]<-sum(data<=x[j])
   f<-y/length(data)
   if(alpha>0){
      if(exact){
         uci<-lci<-rep(NA,length(y))
         for(j in seq(length(y))){
            out<-binom.test(y[j],length(data),conf.level=1-alpha)$conf.int
            lci[j]<-out[1]
            uci[j]<-out[2]
         }
      }else{
         se<-sqrt(f*(1-f)/length(data))
         lci<-f-qnorm(1-alpha/2)*se
         uci<-f+qnorm(1-alpha/2)*se
      }
   }
# Next four lines fill out line from just upper left corners.
# I should be able to do this automatically from step function,
# but I can't figure it out.
   xx<-rep(x,rep(2,length(x)))[-c(1,2*length(f))]
   ff<-rep(f,rep(2,length(f)))[-2*length(f)+c(0,1)]
   if(alpha>0){
      llci<-rep(lci,rep(2,length(f)))[-2*length(f)+c(0,1)]
      uuci<-rep(uci,rep(2,length(f)))[-2*length(f)+c(0,1)]
   }
   mmain<-if(alpha>0) "Empirical CDF and Confidence Bounds" else "Empirical CDF"
   xxlab<-"Data values"
   if(!is.na(dataname)) {
      mmain<-paste(mmain,"for",dataname)
      xxlab<-dataname
   }
   if(newplot){
      if(alpha>0){
         plot(range(xx),c(0,1),type="n",
            main=mmain,xlab=xxlab,ylab="Probability",
            sub=paste("Confidence Level",1-alpha,
               if(exact) "Exact" else "Naive Normal")
         )
      }else{
         plot(range(xx),c(0,1),type="n",
            main=mmain,xlab=xxlab,ylab="Probability")
      }
   }
   for(jj in seq(length(xx)/2)) 
      segments(xx[2*jj-1],ff[2*jj-1],xx[2*jj],ff[2*jj],lty=startlty+1)
   for(ii in seq(length(data))){
      points(xx[2*ii+1],ff[2*ii+1],pch=16)
      points(xx[2*ii],ff[2*ii],pch=1)
   }
   if(alpha>0){
      lines(xx,llci,lty=startlty+2)
      lines(xx,uuci,lty=startlty+2)
      legend(quantile(x,.75),.25,lty=startlty+(1:2),
         legend=c("Empirical CDF","Confidence Bounds"))
   }
}
