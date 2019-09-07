#' @title Tukey HSD procedure
#' @description 
#' Rank-based method for controlling experiment-wise error. Assume normality of the distribution for the variable of interest.
#' @param resp vector containing the values for the variable of interest.
#' @param grp vector specifying in which group is each observation.
#' @param alpha level of the test.
#' @details
#' The original Tuckey HSD procedure is supposed to be applied for equal sample sizes. However, the \code{tukey.kruskal.test} function performs the Tukey-Kramer procedure that works for unequal sample sizes.
#' @return A logical vector for every combinaison of two groups. \code{TRUE} if the distribution in one group is significantly different from the distribution in the other group.
#' @references
#'J.J. Higgins, (2004), \emph{Introduction to Modern Nonparametric Statistics}, Brooks/Cole, Cengage Learning.
#'
#' @export
#' @importFrom stats qtukey
tukey.kruskal.test<-function(resp,grp,alpha=.05){
  ld<-split(resp,grp)
  nv<-sapply(ld,length)
  rm<-sapply(split(rank(resp),grp),mean)
  ng<-length(rm)
  N<-sum(nv)
  out<-array(NA,c(ng*(ng-1)/2,2))
  dimnames(out)<-list(rep("",ng*(ng-1)/2),NULL)
  rr<-rank(resp)
  count<-0
  # This is similar to higgins.fisher.kruskal.test, except that
  # here a loop is substituted for vector calc.,
  # and we calculate confidence intervals instead of p-values.
  for(i in 1:(ng-1)) for(j in (i+1):ng){
    count<-count+1
    out[count,]<-rm[j]-rm[i]+c(-1,1)*qtukey(1-alpha,ng,N-ng)*
      sqrt(N*(N+1)/24)*sqrt(1/nv[i]+1/nv[j])
    dimnames(out)[[1]][count]<-paste(i,j,sep="-")
  }
  different<-(out[,1]>0)|(out[,2]<0)
  names(different)<-dimnames(out)[[1]]
  return(different)
}
