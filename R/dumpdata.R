f <- function(ncp,level,df,targetpower) pchisq(qchisq(1-level,df),df,ncp)-(1-targetpower)
