#' @title Noncentrality Parameter for a Given Level and Power
#' @description This function calculates the noncentrality parameter required to give a test whose null distribution is central chi-square and whose alternative distribution is noncentral chi-square the required level and power.
#' @param df Common degrees of freedom for null and alternative distributions.
#' @param level Level (that is, type I error rate) for the test.
#' @param targetpower Desired power
#' @return required noncentrality parameter.
#' @examples
#' solvencp(4)
#' @export
#' @importFrom stats qchisq
#' @importFrom stats pchisq
#' @importFrom stats uniroot
solvencp<-function(df,level=0.05,targetpower=0.8){
   f<-function(ncp,level,df,targetpower){
      out<-pchisq(qchisq(1-level,df),df,ncp)-(1-targetpower)
      return(out)
   }
   ii<-0
   done<-FALSE
   while((!done)&(ii<8)){
      if(f(2^ii-1,level,df,targetpower)<0){
         done<-TRUE
      }else{
         ii<-ii+1
      }
   }
   out<-if(ii>=8){
      NA
   }else{
      uniroot(f,c(2^(ii-1)-1,2^(ii)-1),level,df,targetpower)$root
   }
   return(out)
}
#solvencp(4)
