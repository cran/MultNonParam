#' @title Generalization of Wilcoxon signed rank test
#' @description This function returns either exact or asymptotic p-values for score tests of the null hypothesis of univariate symmetry about 0.
#' @param y Vector of data on which test will be run.
#' @param scores Scores to be used for the test.  Defaults to integers 1:length(y).
#' @param exact Logical variable indicating whether the exact p-value should be calculate.  Default is false.
#' @param sides Integer; 1 for one sided test rejecting for large values of the statistic, and 2 for the two-sided test.  Defaults to 1.
#' @details The statistic considered here is the sum of scores corresponding to those entries in y that are positive.  If exact=T, the function calls a Fortran code to cycle through all permutations.  If exact=F, the expectation of the statistic is calculated as half the sum of the scores, the variance is calculated as one quarter the sum of squares of scores about their mean, and the statistic is compared to its approximating normal distribution.
#' @return A list with components pv, the p-value obtained with the permutation tests, and tot, the total number of rearrangements of the data considred in calculating the p-value.
#' @examples
#' symscorestat(y=c(1,-2,3,-4,5),exact=TRUE)
#' @references 
#' J.J. Higgins, (2004), \emph{Introduction to Modern Nonparametric Statistics}, Brooks/Cole, Cengage Learning.
#' @export
#' @importFrom stats pnorm
#' @useDynLib MultNonParam signtestperm
symscorestat<-function(y,scores=NULL,exact=F,sides=1){
   if(is.null(scores)) scores<-seq(length(y))
   y<-y[order(abs(y))]
   scores<-sort(scores)
   stat<-sum(scores[y>0])
   if(exact){
      out<-.Fortran("signtestperm",
         as.double(y),
         as.double(scores),
         as.integer(length(y)),
         out=as.integer(0),
         verbose=as.logical(FALSE)
         )$out
#     cat("out",out)
      pvo<-out*2^(-length(y))
   }else{
      zstat<- -(stat-sum(scores)/2)/
         sqrt(sum(scores^2)/4)
      pvo<-pnorm(-zstat)
   }
   if(sides==2) pvo<-2*min(pvo,1-pvo)
   out<-c(pvo,2^length(y),stat)
   names(out)<-c("pv","tot","stat")
   return(out)
}
