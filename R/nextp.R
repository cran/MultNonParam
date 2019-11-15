#' @title Next permutation
#' @description Cycles through permutations of first argument
#' @param perm indices to be permutedj
#' @param b number to begin at.  Set equal to 1.
#' @return The next permutation
#' @export
nextp<-function(perm,b=1){
      out<-.Fortran("nextp",perm=as.integer(perm),as.integer(length(perm)),as.integer(b),PACKAGE="MultNonParam")
#    out<-1
     return(out$perm)
}
