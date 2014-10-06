center <- function(x){(x-mean(x))}
normalize.0m1s <- function(x){(x-mean(x))/sd(x)}
eeSNP.regress <- function( x, y, feature, alpHa=NULL, beTa=rep(0, dim(x)[2]),  estimate.alpha=T, estimate.beta=T, B.inv.alpHa=rep(1/100,  dim(feature)[2]+1),itermax = 5000, thin=40, burnIn=1000, oversample=F, Debug =F,...)
{
    # oversample is currently not working. Avoid it. 
    #if(Debug) set.seed(42)
  NG <- FALSE
    require(MASS); require(BayesLogit); require(geoR); require(SamplerCompare)
    x[is.na(x)] <- 0; y <- center(y); 
    feature <- t(t(feature) - colMeans(feature)); feature <-cbind(1, feature)
    x <- t(t(x) - colMeans(x)) 
    p=numSNP=dim(x)[2]; n=numSamp=dim(x)[1]; numSamp.logit= floor(numSNP/2); f <- dim(feature)[2]
    stopifnot(!(is.null(alpHa) && (!estimate.alpha))); stopifnot(!(is.null(beTa) && (!estimate.beta)))
    if(is.null(alpHa)) alpHa <- rnorm(f) 
    gamMa.prob <- c(logistic(feature %*% alpHa ))
    #gamMa <-  sapply(gamMa.prob, function(pp) rbinom(n=1, size=1, prob=pp) );
    gamMa <- numeric(p); 
    #gamMa[gamMa.prob >= .5] <- 1
    if(is.null(B.inv.alpHa)) B.inv.alpHa <- rep(1/100,  f) 

    beTa.mat <- matrix(0, ncol=p, nrow=ceiling((itermax - burnIn)/thin));
    alpHa.mat <- matrix(0, ncol=length(alpHa), nrow=ceiling((itermax - burnIn)/thin));
    gamMa.mat <- matrix(0, ncol=p, nrow=ceiling((itermax - burnIn)/thin));
    gamMa.prob.mat <- matrix(0, ncol=p, nrow=ceiling((itermax - burnIn)/thin));

    if(Debug) browser()
    w <- rep(0, p)
    varY <- var(y)
    lambdal  <- 4*varY/9; # estimate of sigma2
    seqv  <- seq(0,50 , .05) 
    seqT <- sapply(seqv, function(vv) pinvgamma(varY, shape=vv/2, scale=lambdal*vv/2) ) - .9
    seqT[seqT < 0]  <- 100
    nu <- seqv[which.min(seqT)]
   
   
     
    count  <- 1; tau2 <- .003; c2 <- 10000;
    sigma2 <- c(varY)
    #v1 <- (max(y) - min(y))/(2*3); #delta(y)/3*maxdel(x) 
    #v1 <-c(var(y))/(9*var(x)); #delta(y)/3*maxdel(x) 

    if(NG) v1 <- varY/(36* mean(var(x)))
    else v1  <- 1/(16*mean(var(x))) #v1 <- var(y)/(9* mean(var(x)) # on expected 1/6 of variance of y will be explained by a SNP 
    delta <- .05*sd(y)/2; # mindelta(y)/maxdelta(x)
    #eqn 8 of Geroge et. al. (approaches for BVS)
    minR <- 10; maxR <- 1e12
    delc <-function(a0, a1) log( log(a1/a0) / (1/a0 - 1/a1))
    seqc  <-exp( log(minR) + seq(0,1,.0001)* log(maxR/(minR)))
    seqF <- sapply(seqc,function(cc) delc(v1/cc, v1)) - 2*log(delta) 
    v0  <- v1/seqc[which.min(abs(seqF))]

 
    logit.freq <- rep(1,p); gamMa.inclusion <- gamMa; zz  <- rep(0,p); sample.prob <-  rep(1,p)
    x.covar <- colSums(x^2);
    #x.invcovar <- 1 /(x.covar + v1);
    xTy <- t(x) %*% y
    xbeTa <- x %*% beTa
    TT1 <- TT0 <- TT  <- 0


    #initialize
   #prior <-  sum(log(1-gamMa.prob) ) 
    cc <- 10 # g-priors
    epiDiscount  <- .1
    
    for(iter in seq(itermax)){
	print(iter)
    if(NG){
	if(estimate.beta){
	    #gibbs sampling
	    if(Debug) browser()
	    localmaxiter <- 1
	    if(oversample) logit.freq[] <- 0; 
	    gamMa.inclusion[] <- 0
	    for(localiter in seq(localmaxiter)){
		zero <- which(gamMa==0); nonzero <- which(gamMa!=0); k  <- length(nonzero); x.nonzero <- x[,nonzero, drop=F]
		#Chol <- chol(t(x.nonzero) %*% x.nonzero + diag(1/sigma2A, k,k))
		#xty <-  t(x.nonzero) %*% y
		#Ctxty <- solve( t(Chol), xty) 
		if(k >0){ 
		    #omega <- ginv(t(x.nonzero) %*% x.nonzero + diag(1/sigmaA[nonzero]^2))
		    omega <- ginv(t(x.nonzero) %*% x.nonzero + diag(v1, k, k))
		    omegaxy <- omega %*% t(x.nonzero) %*% y
		    beTa[nonzero]  <- mvrnorm(1, mu = omegaxy, Sigma = omega*sigma2)
		    #beTa[nonzero]  <- omegaxy;
		}
		#if(k<p) beTa[zero] <- 0
		if(k<p) beTa[zero] <- rnorm(length(zero), mean=0, sd = sqrt(tau2))
		#browser()
		#sapply(sigmaA[zero], function(ss) rnorm(1, mean =0, sd = ss/10))  
		xbeTa <- x.nonzero %*% beTa[nonzero] 
		residue <- y - xbeTa
		#sigma2 <- rinvchisq( n = 1, df= n - 2 , scale = mean(residue^2) )
		sigma2 <- 1/rgamma( n = 1, shape= n/2 , scale = 2/sum(residue^2) ) 
		for(snp in seq(p)){
		  #if(Debug) browser()
		  # lambda (ridge regression )  = 1/v1 * sigma2 
		    #beTacurr <-  (xTy[snp]- (x[,snp] %*% xbeTa)+ x.covar[snp] * beTa[snp]) /(x.covar[snp] + sigma2/v1 )
		    #beTacurr <-  (xTy[snp]- x[,snp] %*% xbeTa + x.covar[snp] * beTa[snp])*x.invcovar[snp]
		    #beTacurr <- rnorm(1 , mean=beTacurr, sd= sqrt(sigma2 * x.invcovar[snp]))
		    #sigmaA[snp] <- sqrt(rinvchisq( n = 1, df=vv + 1, scale = S + sum(beTa[snp]^2) )  ) 
		    #residue.snp <- residue  + x[,snp] * gamMa[snp] * beTa[snp]
		    #like1 <- sum(dnorm(residue - x[,snp]*beTa[snp] , mean=0, sd=sigma2, log=T))
		    prior1 <- dnorm(beTa[snp], mean=0, sd=sqrt(v1), log=T)
		    #prior1 <- dnorm(beTa[snp], mean=0, sd=sqrt(c2*tau2), log=T)
		    #prior1 <- dnorm(beTacurr, mean=0, sd=sqrt(v1), log=T)
		    post1 <-  prior1 + log(gamMa.prob[snp])
		    #like0 <- sum(dnorm(residue, mean=0, sd=sigma2, log=T))
		    prior0 <- dnorm(beTa[snp], mean=0, sd=sqrt(v0), log=T)
		    #prior0 <- dnorm(beTa[snp], mean=0, sd=sqrt(tau2), log=T)
		    #prior0 <- dnorm(beTacurr, mean=0, sd=sqrt(v0), log=T)
		    post0 <-  prior0 + log(1 - gamMa.prob[snp])
		    llratio <-  post1 - post0
		    #gamMa.prob[snp] <-  logistic( llratio)
		    rate <-  logistic( llratio)
		    gamMa[snp]  <-  rbinom(n=1, size=1, prob=rate)
		    #if(gamMa[snp] == 0) beTacurr = 0
		    #deltacurr <- beTacurr - beTa[snp]
		}
		#if(localiter==10)
		#browser()

		# sample for logistic regression
		zero <- which(gamMa==0); nonzero <- which(gamMa!=0); k  <- length(nonzero); 
		#if(oversample) gamMa.inclusion[nonzero] <- gamMa.inclusion[nonzero] + 1 ;
		##gamMa.inclusion <- gamMa.inclusion + gamMa.prob;
		gamMa.inclusion <- gamMa.inclusion + gamMa;
	    }
	    if(oversample){
		sample.prob[] <- (localmaxiter*p - sum(gamMa.inclusion)) * gamMa.inclusion +  sum(gamMa.inclusion) *(localmaxiter- gamMa.inclusion) 
		logit.freq[] <- as.numeric(tabulate(sample.int(n=p, size=1000, replace=T, prob=sample.prob), p))
		#logit.freq <-  logit.freq + freq.curr
		#logit.freq[] = rep(1, p);
	    }
	    gamMa.inclusion <- gamMa.inclusion/localmaxiter
	}
    }else{
	zero <- which(gamMa==0); nonzero <- which(gamMa!=0); qq  <- length(nonzero); 
	gamMa.prob  <- gamMa.prob * (1-2*epiDiscount) + epiDiscount; 
	prior <- sum( log( gamMa.prob[nonzero] ))  + sum( log(1 - gamMa.prob[zero] )) 
	gg.last  <-  -(n + nu)/2 * log(nu * lambdal + t(y) %*% y ) + prior
	gamMa.inclusion[] <- 0;
       #browser()	
	for(snp in seq(p)){
	    if(Debug) browser()
	    if(gamMa[snp] ==1){
		TT1 <- TT;
	       	#TT0 <- chdd(TT, x[,snp])
		nonzero <- nonzero[!(nonzero == snp)]
		qq <- qq-1
		gg1 <- gg.last
		prior1 <- prior
		prior0 <- prior + log(gamMa.prob[snp]) - log(1 - gamMa.prob[snp])
		if(qq > 0 ){
		    TT0 <- chol( t(x[, nonzero] )%*% x[, nonzero])
		    WW <- backsolve(TT0, xTy[nonzero], transpose=T)
		    gg0 <- -qq/2 * log(1 + cc ) - (n + nu)/2 * log (nu * lambdal + (t(y) %*% y)  - (t(WW) %*% WW)/(1 + 1/cc) ) + prior0
		}else{
		    TT0 <- 0
		    gg0 <- -qq/2 * log(1 + cc ) - (n + nu)/2 * log (nu * lambdal + (t(y) %*% y) ) + prior0
		}	    
	    }else{
		TT0 <- TT; 
		#TT1 <- chud(TT, x[,snp])
		TT1  <-  chol(t(x[, c(nonzero, snp)] )%*% x[, c(nonzero, snp)]) 
		WW <- backsolve(TT1, xTy[c(nonzero, snp)], transpose=T)
		gg0 <- gg.last
		prior0 <- prior
		prior1 <- prior + log(gamMa.prob[snp]) - log(1 - gamMa.prob[snp])
		gg1 <- -(qq+1)/2 * log(1 + cc ) - (n + nu)/2 * log (nu * lambdal + (t(y) %*% y)  - (t(WW) %*% WW)/(1 + 1/cc) ) + prior1
	    }
	    #gg1 <- -qq/2 * log(1 + cc ) - (n + nu)/2 log (nu * lambdal + t(y) %*% y  - (t(WW) %*% W)/(1 + 1/cc) ) + log(piGamma)
	   stopifnot(!is.nan(c(gg0,gg1))) 
	    gamMa.curr <- rbinom(n=1, size=1, prob = logistic(gg1 - gg0))
	   #gamMa.curr <- 1; if(logistic(gg1 - gg0) < .5) gamMa.curr <- 0
	    gamMa[snp]   <- gamMa.curr
	    if(gamMa.curr==1){
		#if(gg.last > gg1) browser()
		nonzero <- c(nonzero, snp)
		prior <- prior1
		gg.last <- gg1
		TT <- TT1
		qq <- qq +1
	    }else{
		#if(gg.last > gg0) browser()
		#stopifnot(gg.last <= gg0)
		prior <- prior0
		gg.last <- gg0
		TT <- TT0
	    }
	    #if(length(nonzero) ==0) browser()

	} 

	gamMa.inclusion <- gamMa.inclusion + gamMa;
	
	##Gibbs sampling for beta and sigma
	zero <- which(gamMa==0); nonzero <- which(gamMa!=0); qq  <- length(nonzero);
	#VV =  1/(1/sigma2 + 1/cc)(X'X)^-1,  mm = VV X'Y    
	TT.inv = backsolve(TT, diag(qq)) 
	VV <- 1/(1/sigma2 + 1/cc) * (TT.inv %*% t(TT.inv) )
	mm <- VV %*% xTy[nonzero]/sigma2
	beTa[nonzero] <- mvrnorm(n=1, mu=mm, Sigma=VV)
	beTa[zero] <- 0
	if(any(abs(beTa) > 50)) browser()	
	xbeTa  <-  (x[, nonzero, drop=F]%*% beTa[nonzero])	
	residue <- y - xbeTa 
	sigma2 <- rinvgamma(n=1, shape= (qq + nu+ n +1), scale= 1/2 * (sum(residue^2) + nu* lambdal))     
	cc <- rinvgamma(n=1, shape=(qq+1)/2, scale = (sum(xbeTa^2) + n) /2)
	#print(cc)	
    }

	#browser();
	if(estimate.alpha){
	    if(Debug) browser()
	    #sample wi
	    freq.nonzero <- which(logit.freq > 0);p.nonzero <- length(freq.nonzero); zz <- logit.freq*(gamMa.inclusion - 1/2)
	    fTalpHa <- feature[freq.nonzero,] %*% alpHa
	    # subsampling is required here due imbalanced dataset
	    #w[freq.nonzero] <- sapply(1:p.nonzero, function(pp)  rpg(num=1, h=logit.freq[freq.nonzero[pp]], z=fTalpHa[pp])) # there can be variable -- fix the length
	    w[freq.nonzero] <- rpg(num=p.nonzero, h=logit.freq[freq.nonzero], z=fTalpHa[1:p.nonzero]) # there can be variable -- fix the length
	    ##calculate m.w V.w
	    diag.w  <- diag(w[freq.nonzero])
	    V.w <- ginv((t(feature[freq.nonzero,]) %*% diag.w %*% feature[freq.nonzero,])+  diag( B.inv.alpHa)) 
	    m.w <- V.w %*% (t(feature[freq.nonzero,]) %*% (zz[freq.nonzero])) # this is like omegafg
	    alpHa[] <- mvrnorm(n=1, mu = m.w, Sigma=V.w)
	    gamMa.prob[] <- logistic(feature %*% alpHa)
	}

	if((iter %%thin==0)  && (iter > burnIn)){
	    gamMa.mat[count,] <- gamMa
	    beTa.mat[count,] <- beTa
	    gamMa.prob.mat[count,] <- gamMa.prob
	    alpHa.mat[count,] <- alpHa
	    count = count +1;
	    
	}

    }
    list(beTa.mat = beTa.mat, alpHa.mat = alpHa.mat, gamMa.mat = gamMa.mat, gamMa.prob.mat=gamMa.prob.mat)
}


enhancer.selection.regress <- function( x, y, feature,   beTa=matrix(0, 3 , dim(feature)[2] + 1), gamMa.prob = matrix(1/3, nrow=dim(x)[2], ncol=3), itermax = 5000, thin=40, burnIn=1000, ...)
{
    require(MASS); require(geoR) 
    x[is.na(x)] <- 0; y <- center(y); feature <-cbind(1, feature)
    p=numSNP=dim(x)[2]; n=numSamp=dim(x)[1]; numSamp.logit= floor(numSNP/2); f <- dim(feature)[2]
    stopifnot(!(is.null(beTa) && (!estimate.beta)))
    gamMa <-  apply(gamMa.prob, 1,  function(pp) sample(c(-1,0,1), size=1, prob=pp) ); 

    numF <- 3*f; 
    beTa.mat <- matrix(0, ncol=numF, nrow=ceiling((itermax - burnIn)/thin));
    gamMa.mat <- matrix(0, ncol=p, nrow=ceiling((itermax - burnIn)/thin));

    w <- rep(0, n)
    count  <- 1; tau2 <- .003; c2 <- 10000;
    logit.freq <- rep(1,p); gamMa.inclusion <- gamMa; zz  <- rep(0,p); sample.prob <-  rep(1,p)
    sigma2 <- sd(y)
    fxg <- matrix(0, nrow=numSamp, ncol=numF)
    fxb <- matrix(0, nrow=numSamp, ncol=numSNP)
    beTa.flat <- numeric(numF); post <- numeric(3)
    for(iter in seq(itermax)){
	print(iter)
	#gibbs sampling
	#browser()
	fxg[] <- 0
	fg <- (feature * gamMa) 
	#fxg <- rowSums(sapply(1:p, function(ss) fg[ x[,ss] + 1,ss]))
	for(samp in seq(numSamp)){
	    for(snp in 0:2){
		fxg[samp, snp*f+seq(f)] <- colSums(fg[which(x[samp,]==snp), ]) 
	    }
	}
	zero <- which(gamMa==0); nonzero <- which(gamMa!=0); k  <- length(nonzero); 
	#Chol <- chol(t(x.nonzero) %*% x.nonzero + diag(1/sigma2A, k,k))
	#xty <-  t(x.nonzero) %*% y
	#Ctxty <- solve( t(Chol), xty) 
	if(k >0){ 
	    #omega <- ginv(t(x.nonzero) %*% x.nonzero + diag(1/sigmaA[nonzero]^2))
	    omega <- ginv(t(fxg) %*% fxg + diag(c2*tau2, numF, numF))
	    omegaxy <- omega %*% t(fxg) %*% y
	    beTa.flat[]  <- mvrnorm(1, mu = omegaxy, Sigma = omega*sigma2)
	}
	beTa[] <- matrix(beTa.flat, ncol=f, byrow=T)
	fb <- beTa %*% t(feature)
	fxb[] <-0 	
	for(samp in seq(numSamp)){
	    for(snp in 0:2){
		inx <- which(x[samp,] == snp)
		fxb[samp, inx ] <- fb[snp+1, inx]
	    }
	}
	#if(k<p) beTa[zero] <- 0
	#sapply(sigmaA[zero], function(ss) rnorm(1, mean =0, sd = ss/10))  
	residue <- y - fxg %*% beTa.flat
	sigma2 <- sqrt(rinvchisq( n = 1, df= n - 2 , scale = mean(residue^2) ))  
	for(snp in seq(p)){
	    residue.curr <-  residue + gamMa[snp] * fxb[, snp]
	    for(include in c(-1,0,1)){
		inx1 <- include + 2 
		post[inx1] <- sum( dnorm(residue.curr - include*fxb[,snp], mean=0, sd=sigma2, log=T)) + log(gamMa.prob[snp,inx1])
	    }
	    post <- exp(post - max(post))
	    gamMa[snp]  <-  sample(c(-1,0,1), size=1, prob=post)
	}

	if((iter %%thin==0)  && (iter > burnIn)){
	    gamMa.mat[count,] <- gamMa
	    beTa.mat[count,] <-c( t(beTa))
	    count = count +1;
	}

    }
    list(beTa.mat = beTa.mat,  gamMa.mat = gamMa.mat)
}
