symscorestat<-function(y,scores=NULL,exact=F,sides=1){
   if(is.null(scores)) scores<-seq(length(y))
   y<-y[order(abs(y))]
   scores<-sort(scores)
   stat<-sum(scores[y>0])
   if(exact){
      out<-.Fortran("signtestperm",as.double(y),as.double(scores),as.integer(length(y)),out=as.integer(0),verbose=as.logical(FALSE),PACKAGE="MultNonParam")$out
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
