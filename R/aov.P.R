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
     out<-.Fortran("aovp",
                as.integer(n),
                as.integer(permi),
                as.integer(nb),
                as.integer(be),
                as.double(xx),
                tot=as.double(0),pv=as.double(0)
                ,PACKAGE="MultNonParam")
     return(list(pv=out$pv,tot=out$tot))
  }else{
     return(NULL)
  }
}
