exactquantileci<-function(xvec,tau=.5,alpha=0.05){
   n<-length(xvec)
   xvec<-sort(xvec)
   out<-array(NA,c(2,length(tau)))
   for(j in seq(length(tau))){
      ii<-qbinom(alpha/2,n,tau[j])
      out[1,j]<-if(ii>=1) xvec[ii] else -Inf
      ii<-n+1-qbinom(alpha/2,n,1-tau[j])
      out[2,j]<-if(ii<=n) xvec[ii] else Inf
   }
   return(out)
}
