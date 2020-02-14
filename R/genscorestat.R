#' @title Normal-theory two sample scorestatistic.
#' @description Calculates the p-value from the normal approximation to the permutation distribution of a two-sample score statistic.
#' @param scores scores of the data.
#' @param group numeric or character vector of group identities.
#' @param correct half the minimal distance between two potential values of the score statistic.
#' @return Object of class htest containing the p-value.
#' @export
genscorestat<-function(scores,group,correct=0){
   N<-length(group); MV<-table(group)
   refg<-names(MV)[1]
   if(is.numeric(group)) refg<-as.numeric(refg)
   if(length(MV)!=2){
      message("genscorestat works only for two groups")
      out<-NA
   }else{
      abar<-mean(scores)
      ahat<-mean(scores^2)
      vv<-MV[2]*MV[1]*(ahat-abar^2)/(N-1)
      ee<-MV[1]*abar
      stat<-structure(sum(scores[group==refg]),.Names= "Gaussian")
      correct<-abs(correct)*sign(stat-ee)
      parameters<-structure(c(ee,vv),.Names=c("mean","variance"))
      out<-list(null.value=structure(0, .Names = "median difference"),
         alternative="two-sided",method="General Score Test",data.name=NULL,
         statistic=stat,parameters=parameters,
         p.value=2*pnorm(-abs(stat-ee-correct)/sqrt(vv)))
      class(out)<-"htest"
   }
   return(out)
}
