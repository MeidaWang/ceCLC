# ceCLC
A computationally efficient clustering linear combination approach to jointly analyze multiple phenotypes for GWAS.

## Overview
ceCLC tests the association between multiple phenotypes and a genetic variant of interest, which uses the Cauchy combination test to combine all p-values of the CLC test statistics obtained from each possible number of clusters. The test statistic of ceCLC approximately follows a standard Cauchy distribution, so the p-value can be obtained from the cumulative density function without the need for the simulation procedure. 


## An example
We illustrate the usage of ceCLC using simulated data.

```
library(stats)
library(MASS)
library(Matrix)

K <- 20 
L0 <- 20
N <- 5000
model <- 1
beta1 <- 0.0032

# simulate data under alternative
data <- ceCLC_generate_data(K,model,N,beta1)
y <- data[,-ncol(data), drop=FALSE]
x <- data[,ncol(data), drop=FALSE]
pv <- ceCLC(x,y,L0)


# simulate data under the null
y <- ceCLC_generate_typeI(K,model,N)
x <- as.matrix(rbinom(N,2,0.3), ncol=1)
pv <- ceCLC(x,y,L0)

```

## References:
[1] Wang M, Zhang S, Sha Q. A computationally efficient clustering linear combination approach to jointly analyze multiple phenotypes for GWAS.

[2] Sha Q, Wang Z, Zhang X, Zhang S. A clustering linear combination approach to jointly analyze multiple phenotypes for GWAS. Bioinformatics. 2019; 35(8):1373-9.

