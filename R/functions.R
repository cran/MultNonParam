sensitivity.plot<-function(y,sub,stats){
  ra<-range(y)
  xr<-mean(ra)+c(-1,1)*diff(ra)
  outlier<-xr[1]+(0:100)*diff(xr)/100
  base<-stats(y)
  sens<-array(NA,c(length(base),length(outlier)))
  dimnames(sens)[[1]]<-names(base)
  for(j in seq(length(outlier))){
    sens[,j]<-stats(c(y,outlier[j]))-base
  }
  plot(xr,range(sens),type="n",main="Sensitivity",sub=sub,
       ylab="Change in Statistic Value",xlab="New Observation")
  inds<-seq(length(base))
  for(i in inds) lines(outlier,sens[i,],lty=i,col=i)
  legend(min(ra),max(sens),lty=inds,col=inds,
         legend=names(base))
}




dmannwhitney<-function(u,m,n){
  if((u<0)|(u>(m*n))){
    val<-0
  }else{
    if(m==1) val<-1
    if(n==1) val<-1
    if((m>1)&(n>1)) val<-dmannwhitney(u,m-1,n)+dmannwhitney(u-m,m,n-1)
  }
  return(val)
}


betatest<-function(x,y){
  if(length(x)>10) cat("Warning: this will take forever.")
  out<-.Fortran("betatest"
                ,as.integer(length(x))
                ,as.double(x)
                ,as.double(y)
                ,pval=as.double(0)
                ,PACKAGE="MultNonParam"
                )
  return(out$pval)
}



# Tukey HSD method applied to Kruskal-Wallis test.
# x is the response variable, and y indicates group.
tukey.kruskal.test<-function(x,y,alpha=.05){
  ld<-split(x,y)
  yl<-unique(y)
  nv<-sapply(ld,length)
  rm<-sapply(split(rank(x),y),mean)
  ng<-length(rm)
  N<-sum(nv)
  out<-array(NA,c(ng*(ng-1)/2,2))
  dimnames(out)<-list(rep("",ng*(ng-1)/2),NULL)
  rr<-rank(x)
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

# Fisher's LSD method applied to the Kruskal-Wallis test, as presented by Higgins.
higgins.fisher.kruskal.test<-function(x,y,alpha=.05){
  ld<-split(x,y)
  yl<-unique(y)
  nv<-sapply(ld,length)
  rm<-sapply(split(rank(x),y),mean)
  ng<-length(rm)
  N<-sum(nv)
  out<-NULL
  if(kruskal.test(x,y)$p.value<alpha){
    pvfs<-(1-pnorm(abs(outer(rm,rm,"-"))/
                     sqrt(N*(N+1)*outer(1/nv,1/nv,"+")/12)))*2+
      (outer(seq(ng),seq(ng),"-") <=0)
    i1<-outer(seq(ng),rep(1,ng))
    out<-cbind(t(i1)[pvfs<alpha],i1[pvfs<alpha])
  }
  return(out)
}




aov.P<-function(dattab){
  # Rows represent treatments, columns represent blocks
  n<-prod(dim(dattab))
  nb<-dim(dattab)[2]
  ng<-dim(dattab)[1]
  permi<-rep(1:ng,nb)
  be<-seq(nb)*ng
  out<-.Fortran("aovp",
                as.integer(n),
                as.integer(permi),
                as.integer(nb),
                as.integer(be),
                as.double(as.vector(dattab)),
                tot=as.double(0),pv=as.double(0)
                ,PACKAGE="MultNonParam")
  return(list(pv=out$pv,tot=out$tot))
}



wilding<-function(u1,u2,m1,n1,m2,n2){
  out<-.Fortran("wildings",u1=as.integer(u1)
                ,u2=as.integer(u2)
                ,m1=as.integer(m1)
                ,n1=as.integer(n1)
                ,m2=as.integer(m2)
                ,n2=as.integer(n2)
                ,out=as.double(0.0)
                ,PACKAGE="MultNonParam")
  return(out$out)
}


util.jplot<-function(x,y,...){
  newx<-newy<-rep(NA,length(y))
  newy[1]<-y[1]; newx[1:2]<-x[1]
  ry<-diff(range(y))
  begin<-1
  for(j in 2:length(y)){
    if(abs(y[j]-newy[begin])>(ry*.001)){ 
      begin<-begin+1; 
      newy[begin]<-y[j]; 
      newx[2*begin-(1:0)]<-x[j]
    }else{
      newx[2*begin]<-x[j]
    }
  }
  newx<-newx[seq(2*begin)]
  newy<-rep(newy[seq(begin)],rep(2,begin))
  plot(newx,newy,type="n",...)
  for(j in seq(length(newy)/2)){
    if(newx[2*j]==newx[2*j-1]){
      points(newx[2*j],newy[2*j])
    }else{
      lines(newx[2*j-(1:0)],newy[2*j-(1:0)])
    }
  }
  return(invisible(list(x=newx,y=newy)))
}


testve<-function(n,m,k,nsamp=100,delta=0,beta=0,disc=0){
  o<-rep(NA,nsamp)
  i<-rep(c(rep(0,m),rep(1,n)),k)
  str<-rep(1:k,rep(n+m,k))
  for(j in seq(nsamp)){
    x<-rnorm(k*(n+m))
    y<-rnorm(k*(n+m))
    if(disc!=0) y<-round(y/disc)*disc
    ds<-as.data.frame(list(str=str,x=x,y=y+(i-1)*delta+beta*x, i=i))
    out <-probest(ds,"y","i","str","x",0)
    o[j]<-(out$b-.5)/sqrt(out$Vb)
  }
  st<-paste("Delta=",delta,"beta=",beta,"m=",m,"n=",n)
  if(disc!=0) st<-paste(st,"made discrete")
  browser()
  aaa<-qqnorm(o,plot.it=F)
  plot(aaa$x,aaa$y,main="Normal QQ Plot for Distribution of KKW Statistic",
       sub=st,type="p")
  abline(b=1,a=0)
  #  qqline(o)
}


#ds :data set
#resp : response manifest variable vector
#grp  : vector of the variable used to form groups (Advanced/non-advanced for the prostate example, Recurrent/non-recurrent for the breast cancer example)
#str : strata variable vector
probest<-function(ds,resp,grp,str,covs=NULL,delta=NA){
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
  strv<-as.numeric(as.factor(ds[[str]]))
  ustr<-unique(strv)
  out<-.Fortran("probest",
                as.integer(length(resp)), M, as.integer(N), as.integer(grpv),
                as.integer(length(gn)), as.integer(gn),
                as.integer(strv), as.integer(ustr),as.integer(length(ustr)),
                as.double(as.matrix(ds[,resp])), as.double(as.matrix(covariates)),
                as.logical(!is.na(ds[,resp])),as.double(delta),
                b=as.double(rep(0,r)),Vb=as.double(rep(0,r^2))
                ,PACKAGE="MultNonParam"
  )
  return(list(b=out$b,Vb=array(out$Vb,c(r,r))))
}
