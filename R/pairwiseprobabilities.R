#' @title Pairwise probabilities of Exceedence
#' @description \code{pairwiseprobabilities} calculates probabilities of one variable exceeding another,
#' where the variables are independent, and with identical distributions  except for a location shift.
#' This calculation is useful for power of Mann-Whitney-Wilcoxon, Jonckheere-Terpstra, and Kruskal-Wallis testing.
#' @param shifts The offsets for the various populations, under the alternative hypothesis.
#' @param distname The distribution of the underlying observations; normal and logistic are currently supported.
#' @param taylor Logical flag forcing the approximation of exeedence probabilities using a Taylor series.
#' @details Probabilities of particular families must be calculated analytically.
#' @return A matrix with as many rows and colums as there are shift parameters.  Row i and column j give the probability of an observation from group j exceeding one from group i.
#' @export
#' @examples
#' pairwiseprobabilities(c(0,1,2),"normal")
pairwiseprobabilities<-function(shifts,distname=c("normal","logistic"),taylor=FALSE){
   dname<-match.arg(distname)
   out<-array(.5,rep(length(shifts),2))
   if(!taylor){
      for(i in seq(length(shifts)-1)){
         for(j in (i+1):length(shifts)){
            out[i,j]<-switch(distname,
               normal=pnorm((shifts[j]-shifts[i])/sqrt(2)),
               logistic=exp(shifts[j]-shifts[i])*(exp(shifts[j]-shifts[i])-(shifts[j]-shifts[i])-1)/(exp(shifts[j]-shifts[i])-1)^2)
            out[j,i]<-1-out[i,j]
         }
      }
   }else{
      for(i in seq(length(shifts)-1)){
         for(j in (i+1):length(shifts)){
            out[i,j]<-.5+probabilityderiv(distname)*(shifts[j]-shifts[i])
            out[j,i]<-1-out[i,j]
         }
      }
   }
   return(out)
}
#' @title Derivative of pairwise probabilities of Exceedence
#' @description \code{probabilityderiv} calculates derivatives probabilities of one variable exceeding another,
#' where the variables are independent, and with identical distributions  except for a location shift, at the null hypothesis.
#' This calculation is useful for power of Mann-Whitney-Wilcoxon, Jonckheere-Terpstra, and Kruskal-Wallis testing.
#' @param distname The distribution of the underlying observations; normal and logistic are currently supported.
#' @details Probabilities of particular families must be calculated analytically, and then differentiated.
#' @return The scalar derivative.
#' @export
probabilityderiv<-function(distname=c("normal","logistic")){
   return(switch(distname, normal=1/(2*sqrt(pi)),logistic=1/6))
}
