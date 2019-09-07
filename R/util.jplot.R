#' @title Plot a curve, skipping bits where there is a large jump.
#' @description Plot a curve, skipping bits where there is a large jump.
#' @param x Ordinates to be plotted.
#' @param y Abcissas to be plotted.
#' @param ... Arguents passed directly to plot.
#' @export
#' @importFrom graphics points
#' @importFrom graphics lines
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
