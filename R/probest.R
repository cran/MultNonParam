probest<-function(ds,resp,grp,str=NULL,covs=NULL,delta=NA,correct=FALSE){
  r<-length(resp)
  if(length(delta)==1){if(is.na(delta)) delta<-rep(0,r)}
  if(length(delta)==0) delta<-rep(0,r)
  covariates<-if(length(covs)<2){
    if(is.null(covs)){
      as.double(NULL) 
    }else{
      if(is.factor(ds[,covs])){
        model.matrix(~ds[,covs]-1)[,-1,drop=F]
      }else{
        ds[,covs,drop=F]
      }
    }
  }else{
    ds[,covs]
  }
  M<-as.integer(if(is.null(covs)) 0 else dim(covariates)[2])
  N<-dim(ds[,resp,drop=F])[1]
  grpv<-as.numeric(as.factor(ds[[grp]]))
  gn<-sort(unique(grpv))
  if(length(gn)!=2) cat("I am surprised that there are more than two groups.  Correction will not work.\n")
  if(is.null(str)){
     strv<-rep(1,N)
  }else{
     strv<-as.numeric(as.factor(ds[[str]]))
  }
  ustr<-unique(strv)
  out<-.Fortran("probest",
                as.integer(length(resp)), M, as.integer(N), as.integer(grpv),
                as.integer(length(gn)), as.integer(gn),
                as.integer(strv), as.integer(ustr),as.integer(length(ustr)),
                as.double(as.matrix(ds[,resp])), as.double(as.matrix(covariates)),
                as.logical(!is.na(ds[,resp])),as.double(delta),
                b=as.double(rep(0,r)),Vb=as.double(rep(0,r^2)),
                correct=as.logical(correct),
                PACKAGE="MultNonParam"
  )
  return(list(b=out$b,Vb=array(out$Vb,c(r,r))))
}
