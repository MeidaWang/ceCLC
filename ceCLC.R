
#-----------------------------------------------------------------
#' In the ceCLC function,
#' x: represents the genotype (number of minor alleles) vector.          
#' y: represents the trait value matrix. Each row represents             
#'    an individual and each column represents a phenotype.              
#' L0: represents the number of clusters.                                
#-----------------------------------------------------------------

###score test####################################
score=function(x,y)
{
  n=nrow(y)
  yvar=apply(y,2,var)
  Tstat=n*cov(y,x)/sqrt(n*yvar*as.numeric(var(x)))
  return(Tstat)
}

###ceCLC test####################################
ceCLC=function(x,y,L0)
{
  Tstat=score(x,y)
  Sigma=cor(y)
  dist=1-Sigma
  hc=hclust(as.dist(dist))
  pv0=rep(9999,L0)
  U=list()
  for (L in c(1:L0))
  {	
    index=cutree(hc,L)
    B=sapply(1:L, function(t) as.numeric(index==t))
    W=t(B)%*%ginv(Sigma)
    U[[L]]=t(W)%*%ginv(W%*%Sigma%*%t(W))%*%W
    CLC=t(Tstat)%*%U[[L]]%*%Tstat
    pv0[L]=1-pchisq(CLC,L)
  }
  ACAT=sum(tan((0.5-pv0[1:L0])*pi))/L0
  pvalue=pcauchy(ACAT, lower.tail = F)
  #pvalue=0.5-atan(ACAT)/pi
  return(pvalue)
}


