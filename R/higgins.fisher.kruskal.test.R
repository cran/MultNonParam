#' @title Fisher's LSD method applied to the Kruskal-Wallis test
#' @description This function applies a rank-based method for controlling experiment-wise error. Two hypothesis have to be respected: normality of the distribution and no ties in the data. The aim is to be able to detect, among k treatments, those who lead to significant differencies in the values for a variable of interest.
#' @param resp vector containing the values for the variable of interest.
#' @param grp vector specifying in which group is each observation.
#' @param alpha level of the test.
#'
#' @details First, the Kruskal-Wallis test is used to test the equality of the distributions of each treatment. If the test is significant at the level \code{alpha}, the method can be applied.
#'
#' @return  A matrix with two columns. Each row indicates a combinaison of two groups that have significant different distributions.
#'
#' @references
#'  J.J. Higgins, (2004), \emph{Introduction to Modern Nonparametric Statistics}, Brooks/Cole, Cengage Learning.
#' @export
#' @importFrom stats pnorm
#' @importFrom stats kruskal.test
higgins.fisher.kruskal.test<-function(resp,grp,alpha=.05){
  ld<-split(resp,grp)
  grpl<-unique(grp)
  nv<-sapply(ld,length)
  rm<-sapply(split(rank(resp),grp),mean)
  ng<-length(rm)
  N<-sum(nv)
  out<-NULL
  if(kruskal.test(resp,grp)$p.value<alpha){
    pvfs<-(1-pnorm(abs(outer(rm,rm,"-"))/
                     sqrt(N*(N+1)*outer(1/nv,1/nv,"+")/12)))*2+
      (outer(seq(ng),seq(ng),"-") <=0)
    i1<-outer(seq(ng),rep(1,ng))
    out<-cbind(t(i1)[pvfs<alpha],i1[pvfs<alpha])
  }
  return(out)
}
