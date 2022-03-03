
### Data simulation 

#---------- -------------------------------------------------------------
# Generate data under null hypothesis. All phenotypes are quantitative. #
# K: represents the number of phenotypes.                               #
# model: represents the model which is described in the paper.          #
# N: represents sample size.                                            #
#------------------------------------------------------------------------
ceCLC_generate_typeI=function(K,model,N)
{
  rho=0.6
  if (model==1) 
  {
    
    c=0.5
    GAMMA=matrix(rep(1,K),ncol=1)
    F=rnorm(N,0,1)
  }
  if (model==2)  
  {
    
    c=0.5
    R=2
    GAMMA=bdiag(rep(list(rep(1,K/R)),R))
    mu=rep(0,R)
    Sigma=(1-rho)*diag(R)+rho*matrix(1,R,R)
    F=mvrnorm(N,mu,Sigma)
  }
  if (model==3) 
  {
    
    c=0.5
    R=5
    GAMMA=bdiag(rep(list(rep(1,K/R)),R))
    mu=rep(0,R)
    Sigma=(1-rho)*diag(R)+rho*matrix(1,R,R)
    F=mvrnorm(N,mu,Sigma)
  }
  if (model==4) 
  {
    
    c=0.5
    R=5
    GAMMA=bdiag(rep(list(rep(1,K/R)),R))
    mu=rep(0,R)
    Sigma=(1-rho)*diag(R)+rho*matrix(1,R,R)
    F=mvrnorm(N,mu,Sigma)	
  }
  E=mvrnorm(N,rep(0,K),diag(K))
  y=c*GAMMA%*%t(F)+sqrt(1-c^2)*t(E)
  as.matrix(t(y))
}


#---------- --------------------------------------------------------------
# Generate data under null hypothesis. Half phenotypes are quantitative, #
# the other half phenotypes are qualitative(binary).                     #
#-------------------------------------------------------------------------
ceCLC_generate_typeI_quan_qual=function(K,model,N)
{
  rho=0.6
  n0=N*0.16
  if (model==1) 
  {
    c=0.5#0.1
    GAMMA=matrix(rep(1,K),ncol=1)
    F=rnorm(N,0,1)
    E=mvrnorm(N,rep(0,K),diag(K))
    y=t(c*GAMMA%*%t(F)+sqrt(1-c^2)*t(E))
    index=sample(K,K/2)
    y1=y[,-index]
    y2=y[,index]
    a0=apply(y2,2,function(z) sort(z,decreasing = FALSE))[n0,]
    y3=t(apply(y2,1, function(z) ifelse(z>=a0,1,0)))
    y4=cbind(y1,y3)
  }
  if (model==2) 
  {
    c=0.5
    R=2
    GAMMA=bdiag(rep(list(rep(1,K/R)),R))
    mu=rep(0,R)
    Sigma=(1-rho)*diag(R)+rho*matrix(1,R,R)
    F=mvrnorm(N,mu,Sigma)
    E=mvrnorm(N,rep(0,K),diag(K))
    y=t(c*GAMMA%*%t(F)+sqrt(1-c^2)*t(E))
    index=sample(K,K/2)
    y1=y[,-index]
    y2=y[,index]
    a0=apply(y2,2,function(z) sort(z,decreasing = FALSE))[n0,]
    y3=t(apply(y2,1, function(z) ifelse(z>=a0,1,0)))
    y4=cbind(y1,y3)
  }
  if (model==3) 
  {
    
    c=0.5
    R=5
    GAMMA=bdiag(rep(list(rep(1,K/R)),R))
    mu=rep(0,R)
    Sigma=(1-rho)*diag(R)+rho*matrix(1,R,R)
    F=mvrnorm(N,mu,Sigma)
    E=mvrnorm(N,rep(0,K),diag(K))
    y=t(c*GAMMA%*%t(F)+sqrt(1-c^2)*t(E))
    index=sample(K,K/2)
    index1=subset(index, index>(K/5*3) & index<=(K/5*4))
    y[,index1]=-y[,index1]
    y1=y[,-index]
    y2=y[,index]
    a0=apply(y2,2,function(z) sort(z,decreasing = FALSE))[n0,]
    y3=t(apply(y2,1, function(z) ifelse(z>=a0,1,0)))
    y4=cbind(y1,y3)
  }
  if (model==4) 
  {
    c=0.5
    R=5
    GAMMA=bdiag(rep(list(rep(1,K/R)),R))
    mu=rep(0,R)
    Sigma=(1-rho)*diag(R)+rho*matrix(1,R,R)
    F=mvrnorm(N,mu,Sigma)
    E=mvrnorm(N,rep(0,K),diag(K))
    y=t(c*GAMMA%*%t(F)+sqrt(1-c^2)*t(E))
    index=sample(K,K/2)
    index1=subset(index, index>(K/5*2) & index<=(K/5*4))
    y[,index1]=-y[,index1]
    y1=y[,-index]
    y2=y[,index]
    a0=apply(y2,2,function(z) sort(z,decreasing = FALSE))[n0,]
    y3=t(apply(y2,1, function(z) ifelse(z>=a0,1,0)))
    y4=cbind(y1,y3)	
  }
  as.matrix(y4)
}


#---------- --------------------------------------------------------------------
# Generate data under alternative hypothesis. All phenotypes are quantitative. #
# K: representss the number of phenotypes.                                     #
# N: represents  sample size.                                                  #
# model: represents the model which is described in the paper.                 #
# beta1: represents the effect size of genotype.                               #
#-------------------------------------------------------------------------------
ceCLC_generate_data=function(K,model,N,beta1)
{
	rho=0.6
	G=as.matrix(rbinom(N,2,0.3),ncol=1)
	if (model==1) 
	{
		
		c=0.5
		GAMMA=matrix(rep(1,K),ncol=1)
		LAMBDA=beta1*seq(K)
		F=rnorm(N,0,1)
	}
	if (model==2) 
	{
		
		c=0.5
		R=2
		GAMMA=bdiag(rep(list(rep(1,K/R)),R))
		LAMBDA=c(rep(0,K-K/R),rep(beta1,K/R))
		mu=rep(0,R)
		Sigma=(1-rho)*diag(R)+rho*matrix(1,R,R)
		F=mvrnorm(N,mu,Sigma)
	}
	if (model==3) 
	{
		
		c=0.5
		R=5
		GAMMA=bdiag(rep(list(rep(1,K/R)),R))
		LAMBDA=beta1*c(rep(0,K-2*K/R),rep(-1,K/R),2*seq(K/R)/(K/R+1))
		mu=rep(0,R)
	        Sigma=(1-rho)*diag(R)+rho*matrix(1,R,R)
     		F=mvrnorm(N,mu,Sigma)
	}
	if (model==4) 
	{
		
		c=0.5
		R=5
		GAMMA=bdiag(rep(list(rep(1,K/R)),R))
		LAMBDA=beta1*c(rep(0,K/R),rep(-1,K/R),rep(1,K/R),-2*seq(K/R)/(K/R+1),2*seq(K/R)/(K/R+1))
		mu=rep(0,R)
	        Sigma=(1-rho)*diag(R)+rho*matrix(1,R,R)
     		F=mvrnorm(N,mu,Sigma)	
	}
	E=mvrnorm(N,rep(0,K),diag(K))
	y=LAMBDA%*%t(G)+c*GAMMA%*%t(F)+sqrt(1-c^2)*t(E)
	as.matrix(cbind(t(y),G))
}


#---------- ---------------------------------------------------------------------
# Generate data under alternative hypothesis. Half phenotypes are quantitative, #
# the other half phenotypes are qualitative(binary).                            #
#--------------------------------------------------------------------------------
ceCLC_generate_data_quan_qual=function(K,model,N,beta1)
{
	rho=0.6
	G=as.matrix(rbinom(N,2,0.3),ncol=1)
	n0=N*0.16
	if (model==1) 
	{
		c=0.5
		GAMMA=matrix(rep(1,K),ncol=1)
		LAMBDA=beta1*seq(K)
		F=rnorm(N,0,1)
		E=mvrnorm(N,rep(0,K),diag(K))
		y=t(LAMBDA%*%t(G)+c*GAMMA%*%t(F)+sqrt(1-c^2)*t(E))
		index=c(1:(K/2)) 
		y1=y[,-index]
		y2=y[,index]
		a0=apply(y2,2,function(z) sort(z,decreasing=TRUE))[n0,]
		y3=t(apply(y2,1, function(z) ifelse(z>=a0,1,0)))
		y4=cbind(y3,y1)
	}
	if (model==2) 
	{
		c=0.5
		R=2
		GAMMA=bdiag(rep(list(rep(1,K/R)),R))
		LAMBDA=c(rep(0,K-K/R),rep(beta1,K/R))
		mu=rep(0,R)
		Sigma=(1-rho)*diag(R)+rho*matrix(1,R,R)
		F=mvrnorm(N,mu,Sigma)
		E=mvrnorm(N,rep(0,K),diag(K))
		y=t(LAMBDA%*%t(G)+c*GAMMA%*%t(F)+sqrt(1-c^2)*t(E))
		index=c(1:(K/2)) 
		y1=y[,-index]
		y2=y[,index]
		a0=apply(y2,2,function(z) sort(z,decreasing=TRUE))[n0,]
		y3=t(apply(y2,1, function(z) ifelse(z>=a0,1,0)))
		y4=cbind(y3,y1)
	}
	if (model==3) 
	{
		
		c=0.5
		R=5
		GAMMA=bdiag(rep(list(rep(1,K/R)),R))#
		LAMBDA=beta1*c(rep(0,K-2*K/R),rep(-1,K/R),2*seq(K/R)/(K/R+1))
		mu=rep(0,R)
	    Sigma=(1-rho)*diag(R)+rho*matrix(1,R,R)
     	F=mvrnorm(N,mu,Sigma)
     	E=mvrnorm(N,rep(0,K),diag(K))
		y=t(LAMBDA%*%t(G)+c*GAMMA%*%t(F)+sqrt(1-c^2)*t(E))
		index=c(1:(K/2)) 
		index1=subset(index, index>(K/5*3) & index<=(K/5*4))
		y[,index1]=-y[,index1]
		y1=y[,-index]
		y2=y[,index]
		a0=apply(y2,2,function(z) sort(z,decreasing=TRUE))[n0,]
		y3=t(apply(y2,1, function(z) ifelse(z>=a0,1,0)))
		y4=cbind(y3,y1)
	}
	if (model==4) 
	{
		c=0.5
		R=5
		GAMMA=bdiag(rep(list(rep(1,K/R)),R))
		LAMBDA=beta1*c(rep(0,K/R),rep(1,K/R),rep(-1,K/R),-2*seq(K/R)/(K/R+1),2*seq(K/R)/(K/R+1))
		mu=rep(0,R)
	    Sigma=(1-rho)*diag(R)+rho*matrix(1,R,R)
     	F=mvrnorm(N,mu,Sigma)
     	E=mvrnorm(N,rep(0,K),diag(K))
		y=t(LAMBDA%*%t(G)+c*GAMMA%*%t(F)+sqrt(1-c^2)*t(E))
		index=c(1:(K/2))
		index1=subset(index, index>(K/5*2) & index<=(K/5*4))
		y[,index1]=-y[,index1]
		y1=y[,-index]
		y2=y[,index]
		a0=apply(y2,2,function(z) sort(z,decreasing=TRUE))[n0,]
		y3=t(apply(y2,1, function(z) ifelse(z>=a0,1,0)))
		y4=cbind(y3,y1)	
	}
	as.matrix(cbind(y4,G))
}

