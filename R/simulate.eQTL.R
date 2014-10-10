logistic <- function(x) 1/(1+exp(-x)) 


generate.epieQTL.logistic <- function ( numSNP = 512, avg.eSNP=20,numGenes = 1, numSamp=20, numFeature=9, featureStat= list(enhancer=list(mean=rep(2,numFeature), cov= 2* diag(numFeature)), random=list(mean=rep(1,numFeature), cov=4* diag(numFeature))), seed=NULL,   evSNP =0.025){
    if (!is.null(seed)) set.seed(seed)
    #numSNP  <- numSNP
    numSNPAll <- ceiling(numSNP *( numGenes/2 + .5))
    ldblock =2 
    indep.SNP <- ceiling(numSNPAll/ldblock)
    x.indep <- matrix(sample(x=c(0,1,2,NA) ,size= indep.SNP*numSamp, replace=T , prob=c(.6, .25, .15, .002)), ncol=indep.SNP, nrow=numSamp)
    x <- matrix( rbind(x.indep,x.indep), nrow=numSamp)[,1:numSNPAll]  
    x.nona = x
    #cat("here")
    x.nona[is.na(x)]=0
    ev <-rep(.75, numGenes)  # rbeta(n=numGenes, 3.5, 5)
    nenhancer	<- ceiling(ev/evSNP )
    mask2 <- matrix(0, numGenes*numSNP,2)
    xx <- matrix(0, numSamp, numGenes*numSNP)
    enhancer.out <- enhancer  <- enhancerLDblock <- list()
    beTa <- gamMa <- numeric(numGenes* numSNP)
    y <- matrix(0, numSamp, numGenes)
    for(gene in seq(numGenes)){
	gene.start  <- numSNP * (gene -1) + 1
	gene.end  <- gene.start +  numSNP - 1
	xinx.start  <- ceiling(numSNP/2 * (gene -1) + 1)
	xinx.end  <- xinx.start +  numSNP - 1
	enh <- sample(seq(1, numSNP, ldblock), size=nenhancer[gene])
	beTa1 <- numeric(numSNP)
	beTa1[enh] <- rnorm(n=length(enh),sd=2) 
	xx[,gene.start:gene.end] = x[, xinx.start:xinx.end]
	y.predict <- x.nona[, xinx.start:xinx.end] %*% beTa1
	fac <- c(sqrt(ev[gene] / var(y.predict) ))
	y[,gene] = y.predict * fac + rnorm(n=numSamp, sd=sd(y.predict)) * sqrt(1-ev[gene] )
	beTa1 <- beTa1*fac
	beTa[gene.start:gene.end] = beTa1
	gamMa[gene.start + enh -1]  <- 1
	mask2[gene.start:gene.end,] = cbind(rep(gene, numSNP), gene.start:gene.end)
	enhancer[[gene]]  <- xinx.start + enh -1 
	enhancer.out[[gene]] <- gene.start + enh -1
	enhancerLDblock[[gene]]  <- c( xinx.start + enh -1,  xinx.start + enh ) 
    }
    #browser()
    #noneeSNP.enhancer  <- sample(seq(1, numSNPAll,2), size=nenhancer*numGenes*5, replace=F)  
    #enhancer1  <-  unique(c(unlist(enhancer), noneeSNP.enhancer)) 
    enhancerLDblock  <-  c(unlist(enhancerLDblock)); 
    snpinx <- 1:numSNPAll;
    temp1 <- snpinx[!(snpinx %in% enhancerLDblock)]
    if (length(temp1) < sum(nenhancer) ) stop("increase number of SNP per gene, too many enhancers")
    noneeSNP.enhancer  <- sample(temp1, size=sum(nenhancer), replace=F)  
    enhancer1  <-  unique(c(unlist(enhancer), noneeSNP.enhancer)) 
    
    nonenhancer <- which( !(seq(numSNPAll) %in% enhancer1))
    feature = matrix(0, nrow=numSNPAll, ncol=(numFeature))
    feature[enhancer1, ] <- mvrnormArma(n=length(enhancer1), mu= featureStat$enhancer$mean, sigma= featureStat$enhancer$cov )
    feature[nonenhancer, ] <- mvrnormArma(n=length(nonenhancer), mu= featureStat$random$mean, sigma= featureStat$random$cov )
    #feature[nonenhancer, ] <- sapply(1:numFeature, function(ee) rnorm(n=length(nonenhancer), mean=featureStat$random$mean[ee], sd=featureStat$random$sd[ee]) ) 
    #feature[enhancer1, ] <- sapply(1:numFeature, function(ee) { rnorm(n=length(enhancer1), mean=featureStat$enhancer$mean[ee], sd=featureStat$enhancer$sd[ee])} ) 
    #feature[nonenhancer, ] <- rnorm(n=length(nonenhancer)*(numFeature), mean=1, sd=4)  
    #feature[enhancer1, ] <- rnorm(n=length(enhancer1)*(numFeature), mean=2, sd=2)
    featureAll <- matrix(0, numSNP*numGenes, numFeature)
    theta.temp  <-  numeric(numSNPAll)
    theta.temp[enhancer1]  <- 1  
    theta  <-  numeric(numSNP*numGenes)
    for(gene in seq(numGenes)){
	gene.start  <- numSNP * (gene -1) + 1
	gene.end  <- gene.start +  numSNP - 1
	xinx.start  <- ceiling(numSNP/2 * (gene -1) + 1)
	#print(xinx.start )
	xinx.end  <- xinx.start +  numSNP - 1
	featureAll[gene.start: gene.end, ] = feature[xinx.start:xinx.end, ] 
	theta[gene.start: gene.end]  <-  theta.temp[xinx.start:xinx.end]
    }
    #browser()
    en.inx = which(theta > 0)
    if (length(en.inx) >  nrow(featureAll)/2) stop("increase number of SNP per gene, too many enhancers")
    nenh =which(!( (1:nrow(featureAll) ) %in% en.inx))
    inx = sample(nenh, length(en.inx))
    #en.inx  <-  en.inx[1:200]
    #inx <- inx[1:200]
    yy = c(rep(1, length(en.inx)) , rep(0, length(en.inx)))
    XX = rbind(featureAll[en.inx,], featureAll[inx,]) 
    p.fit <- BayesLogit::logit(y=yy, X = cbind(1, XX), P0=1/100*diag(numFeature +1))
    alpHa <- apply(p.fit$beta,2, median)
    const <- cbind(1, featureAll) %*% matrix( alpHa, ncol=1)
    gamMa.prob <- logistic(const)
    #confusion.matrix 
    object = c(list(x=xx, y=y, y.predict=y.predict, feature=featureAll, alpHa=alpHa, gamMa=gamMa, gamMa.prob=gamMa.prob, beTa=beTa, enhancer = enhancer.out, theta=theta, mask2 = mask2))
    class(object) = "eeSNP.model5"
    return(object)
}
  
