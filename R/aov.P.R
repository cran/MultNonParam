#' @title One-way ANOVA using permutation tests
#' @description \code{aov.P} uses permutation tests instead of classic theory tests to run a one-way or two-way ANOVA.
#' @param dattab The table on which the ANOVA has to be done, or a vector of responses.
#' @param permi If dattab is a table, ignored.  If dattab is a vector, a vector of treatment labels.
#' @param be If dattab is a table, ignored.  If dattab is a vector, a vector of end points of blocks.  In this case, blocks must form contiguous subvectors of dattab.  If null, no blocking.
#' @details The function calls a Fortran code to perform the permutation tests and the ANOVA.  The function has to be applied directly on a cross-table of two variables.
#' @return A list with fields pv, the p-value obtained with the permutation tests, and tot, the total number of permutations. 
#' @export
#' @useDynLib MultNonParam aovp
aov.P<-function(dattab,permi=NULL,be=NULL){
# Rows represent treatments, columns represent blocks
  bad<-F
  if(is.matrix(dattab)){
     n<-prod(dim(dattab))
     nb<-dim(dattab)[2]
     ng<-dim(dattab)[1]
     permi<-rep(1:ng,nb)
     be<-seq(nb)*ng
     xx<-as.vector(dattab)
  }else{
     n<-length(dattab)
     if(is.null(permi)){
        cat("\n Error in aov.P.  permi cannot be null if dattab is a vector\n")
        bad<-T
     }
     permi<-as.numeric(as.factor(permi))
     xx<-dattab
     if(is.null(be)) be<-n
     ng<-length(unique(permi))
     nb<-length(be)
  }
  if(!bad){
     out<-.Fortran(aovp,
                as.integer(n),
                as.integer(permi),
                as.integer(nb),
                as.integer(be),
                as.double(xx),
                tot=as.double(0),pv=as.double(0))
     return(list(pv=out$pv,tot=out$tot))
  }else{
     return(NULL)
  }
}
