page.test.unbalanced<-function(x,trt,blk,sides=2){
#' Perform Page test for unbalanced two-way design
#' @param x A vector of responses
#' @param trt A vector of consecutive integers starting at 1 indicating treatment
#' @param blk A vector of consecutive integers starting at 1 indicating block
#' @param sides A single integer indicating sides.  Defaults to 2.
#' @return P-value for Page test.
#' @examples
#' page.test.unbalanced(rnorm(15),rep(1:3,5),rep(1:5,rep(3,5)))
   makecmat<-function(trt,blk){
      nmat<-table(trt,blk)
      cmat<-array(NA,rep(length(unique(trt)),2))
      nt<-length(unique(trt))
      mvec<-rep(NA,nt)
      nb<-length(unique(blk))
      nsl<-apply(nmat,2,sum)
      nsk<-apply(nmat,1,sum)
      mvec<-rep(0,nt)
      for(k in seq(nt)){
         for(l in seq(nb)) mvec[k]<-mvec[k]+nmat[k,l]*(nsl[l]+1)
      }
      mvec<-mvec/(2*nsk)
      for(k in seq(nt-1)) for(m in ((k+1):nt)) {
         cmat[k,m]<-0
         for(l in seq(nb))
            cmat[k,m]<-cmat[k,m]-nmat[k,l]*nmat[m,l]*(nsl[l]+1)
         cmat[k,m]<-cmat[k,m]/(12*nsk[k]*nsk[m])
         cmat[m,k]<-cmat[k,m]
      }
      for(k in seq(nt)){
         cmat[k,k]<-0
         for(l in seq(length(unique(blk))))
            cmat[k,k]<-cmat[k,k]+nmat[k,l]*sum(nmat[-k,l])*(nsl[l]+1)
         cmat[k,k]<-cmat[k,k]/(12*nsk[k]^2)
      }
      return(list(mvec=mvec,cmat=cmat))
   }
   pagestatu<-function(x,trt,blk){
      ranksums<-rep(0,length(unique(trt)))
      names(ranksums)<-unique(trt)
      for(oneblk in unique(blk)){
         rr<-rank(x[blk==oneblk])
         thesetrt<-trt[blk==oneblk]
         for(onet in unique(trt)){
            ranksums[onet]<-ranksums[onet]+rr[thesetrt==onet]
         }
      }
      return(ranksums/table(trt))
   }
   fun.checkvar<-function(trt,blk,nsamp=10000){
      ranksums<-rep(0,length(unique(trt)))
      names(ranksums)<-unique(trt)
      rankmeans<-array(NA,c(nsamp,length(unique(trt))))
      for(m in 1:nsamp){
         x<-rnorm(length(trt))
         rankmeans[m,]<-pagestatu(x,trt,blk)
      }
      return(list(variance=var(rankmeans),mean=apply(rankmeans,2,mean)))
   }
# First executable
   rm<-pagestatu(x,trt,blk)
   nt<-length(rm)
   moms<-makecmat(trt,blk)
   z<-seq(nt)%*%(rm-moms$mvec)/sqrt(seq(nt)%*%moms$cmat%*%seq(nt))
   pval<-pnorm(z,lower.tail=F)
   if(sides==2) pval<-min(pval,1-pval)*2
   return(list(z=z,pval=pval))
}
#fun.checkvar(trt,blk)
