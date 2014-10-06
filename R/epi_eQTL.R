find.correlation  <- function(mask, x, y)
{
  yLen <- nrow(y); xLen <- nrow(x);
  corr <- list()
  for(yy in seq(yLen)){
    if(length(mask[[yy]] )> 0) corr[[yy]] <- apply(x[mask[[yy]], ,drop=F] , 1, cor, y=y[yy,]) 
  }
  corr
}

center <- function(x){t(t(x)-colMeans(x))}
normalize.0m1s <- function(x){(x-mean(x))/sd(x)}
epi.eQTL <- function( x, y, feature, pairFeature, mask,    estimate.alpha=T, estimate.beta=T, alpHa=NULL, beTa=NULL, gamMa=NULL, B.inv.alpHa=NULL, itermax = 5000, thin=40, burnIn=1000, oversample=F, Debug =F,...)
{
    # oversample is currently not working. Avoid it. 
    require(MASS); require(BayesLogit); require(geoR); require(SamplerCompare)
    x[is.na(x)] <- 0; y <- center(y); 
    feature <-cbind(1, feature)
    x <- t(t(x) - colMeans(x))
   pairf <-ifelse(is.null(pairFeature),0, ncol(pairFeature[[(which(lapply(mask, length) >0)[1])]]))
   numSNP=dim(x)[2]; n=numSamp=dim(x)[1]; f <- ncol(feature)+pairf; numGenes <- dim(y)[2]; gene.range <- 1:numGenes;
    stopifnot((!is.null(alpHa)) || estimate.alpha);    if(is.null(alpHa)) alpHa <- rnorm(f) 
    x.curr  <- xTy  <-   TT  <- list(); gg <- numeric(numGenes)
    if(Debug) browser()
    defined.gamMa <- defined.beTa <- T;
    if (is.null(gamMa) ){  gamMa <- list(); defined.gamMa <- F}
   if (is.null(beTa) ) {beTa <- list(); defined.beTa <- F} 
    if(is.null(B.inv.alpHa)) B.inv.alpHa <- rep(1/100,  f) 

    numSamples <- ceiling((itermax - burnIn)/thin)
     gamMa.prob.mat <- gamMa.mat <- beTa.mat <- list();
    alpHa.mat <- matrix(0, ncol=length(alpHa), nrow=numSamples);

   
   # needed only  once 
    small <- .01;
    varY <- apply(y , 2, var)
    sigma2  <- lambdal  <- 1*varY/9; # estimate of sigma2
    seqv  <-c( seq(0, 50, .05), exp(seq(-20, 20, .01)))
    nu <-  numeric(numGenes); p <- 0; 
    w <- rep(0, p)
    for(gene in gene.range){
	seqT <- sapply(seqv, function(vv) pinvgamma(varY[gene], shape=vv/2, scale=lambdal[[gene]]*vv/2) ) - .9
	seqT[seqT < 0]  <- 100
	nu[[gene]] <- seqv[which.min(seqT)]
	p <- p+ length(mask[[gene]])
    }

    #initialize
    logit.freq <- rep(1,p); gamMa.inclusion <- zz  <-  rep(0,p); featureCat <- matrix(0, nrow=p, ncol=f)
    cc <- numeric(numGenes); cc[] <- 10 # g-priors
    epiDiscount  <- .02; count <- 1

    for(gene in gene.range){
	if(gene ==1) gene.end <- 0
	gene.start <- gene.end + 1; gene.end <- gene.start+ length(mask[[gene]]) - 1
	p.curr <- length(mask[[gene]])
	if(p.curr > 0){ 
	    #separate x 
	    x.curr[[gene]] <- x[, mask[[gene]], drop=F]
	    if(!defined.gamMa) gamMa[[gene]] <-  numeric(p.curr)
	    if(!defined.beTa) beTa[[gene]] <- numeric(p.curr)
	    beTa.mat[[gene]] <- gamMa.prob.mat[[gene]] <- 	gamMa.mat[[gene]] <-  matrix(0, ncol=p.curr, nrow=numSamples);
	    rownames(beTa.mat[[gene]]) <- rownames(gamMa.mat[[gene]]) <- rownames(x.curr) 
	    featureCat[gene.start:gene.end, ] <- cbind(feature[mask[[gene]],, drop=F] , pairFeature[[gene]])
	    if(!estimate.beta) zz[gene.start:gene.end] <- gamMa[[gene]] - .5;
	    #separate y
	    #y.curr[[gene]] <- y[,gene]
	    xTy[[gene]] <- t(x.curr[[gene]]) %*%y[,gene,drop=F] 
	    gamMa.prob  <- rep(min(.5, 20/p.curr), p.curr)
	    TT[[gene]]  <- 0
	    k <- length(gamMa[[gene]]); zero <- which(gamMa[[gene]]==0); nonzero <- which(gamMa[[gene]]!=0); qq  <- length(nonzero);
	    #prior <- sum( log( gamMa.prob[nonzero] ))  + sum( log(1 - gamMa.prob[zero] )) 
	    gg[gene]  <-  -(n + nu[[gene]])/2 * log(nu[[gene]] * lambdal[[gene]] + t(y[,gene,drop=F]) %*% y[,gene,drop=F] ) 
	    #browser()
	}
    }
    for(iter in seq(itermax)){
	print(iter)
	if(estimate.beta){
	    for(gene in gene.range){
		if(gene ==1) gene.end <- 0
		gene.start <- gene.end + 1; gene.end <- gene.start+ length(mask[[gene]]) - 1
		p.curr <- length(mask[[gene]])
		if(p.curr > 0){ 
		    gamMa.prob <- logistic(featureCat[gene.start:gene.end, ] %*% alpHa)
		    gamMa.prob  <- gamMa.prob * (1-2*epiDiscount) + epiDiscount;

		    k <- length(gamMa[[gene]]); zero <- which(gamMa[[gene]]==0); nonzero <- which(gamMa[[gene]]!=0); qq  <- length(nonzero);
		    #prior <- sum( log( gamMa.prob[nonzero] ))  + sum( log(1 - gamMa.prob[zero] )) 
		    gamMa.rate <- numeric(k)
		    #gamMa.inclusion[[gene]] <- numeric(k); 
		    #browser()
		    for(snp in seq(p.curr)){
			#if(Debug) browser()
			if(gamMa[[gene]][snp] ==1){
			    TT1 <- TT[[gene]];
			    #TT0 <- chdd(TT, x[,snp])
			    nonzero.curr <- nonzero[!(nonzero == snp)]
			    qq.curr <- qq-1
			    gg1 <- gg[gene]
			    #prior1 <- prior
			    #prior0 <- prior + log(gamMa.prob[snp]) - log(1 - gamMa.prob[snp])
			    if(qq.curr > 0 ){
				TT0 <- chol( (t(x.curr[[gene]][, nonzero.curr] )%*% x.curr[[gene]][, nonzero.curr]) + small * diag(qq.curr))
				WW <- backsolve(TT0,xTy[[gene]][nonzero.curr], transpose=T)
				gg0 <- -qq.curr/2 * log(1 + cc[[gene]] ) - (n + nu[[gene]])/2 * log (nu[[gene]] * lambdal[[gene]] + (t(y[,gene,drop=F]) %*% y[,gene,drop=F])  - (t(WW) %*% WW)/(1 + 1/cc[[gene]]) ) 			    }else{
				TT0 <- 0
				gg0 <- -qq.curr/2 * log(1 + cc[[gene]] ) - (n + nu[[gene]])/2 * log (nu[[gene]] * lambdal[[gene]] + (t(y[,gene,drop=F]) %*% y[,gene,drop=F]) ) 
			    }	    
			}else{
			    TT0 <- TT[[gene]]; 
			    nonzero.curr <- sort(c(snp, nonzero) )
			    qq.curr <- qq+1
			    TT1  <-  chol(t(x.curr[[gene]][, nonzero.curr] )%*% x.curr[[gene]][,nonzero.curr]+ small * diag(qq.curr)) 
			    WW <- backsolve(TT1,xTy[[gene]][nonzero.curr], transpose=T)
			    gg0 <- gg[gene]
			    #prior0 <- prior
			    #prior1 <- prior + log(gamMa.prob[snp]) - log(1 - gamMa.prob[snp])
			    gg1 <- -(qq.curr)/2 * log(1 + cc[[gene]] ) - (n + nu[[gene]])/2 * log (nu[[gene]] * lambdal[[gene]] + (t(y[,gene,drop=F]) %*% y[,gene,drop=F])  - (t(WW) %*% WW)/(1 + 1/cc[[gene]]) )
			}
		    #browser()
			stopifnot(!is.nan(c(gg0,gg1)))
			print(cat(gene,  snp, c(WW) ))
			if(gene==1) {
			  print("nonzero")
			print(nonzero.curr)
			}
			gamMa.rate[snp] <- logistic(gg1 - gg0 +  log(gamMa.prob[snp]) - log(1 - gamMa.prob[snp]))	
			#gamMa.curr <- rbinom(n=1, size=1, prob = gamMa.rate[snp])
			gamMa.curr <- 1; if(logistic(gg1 - gg0) < .5) gamMa.curr <- 0
			if(gamMa.curr != gamMa[[gene]][snp] ) {
			    #print("changing")
			    gamMa[[gene]][snp]   <- gamMa.curr
			    if(gamMa.curr==1){
				#if(gg.last > gg1) browser()
				nonzero <- nonzero.curr 
				#prior <- prior1
				gg[gene] <- gg1
				TT[[gene]] <- TT1
				qq <- qq.curr
			    }else{
				#if(gg.last > gg0) browser()
				#stopifnot(gg.last <= gg0)
				nonzero <- nonzero.curr 
				#prior <- prior0
				gg[gene] <- gg0
				TT[[gene]] <- TT0
				qq <- qq.curr
			    }
			}
			#print(nonzero)

			if(qq > 0) {
			    if ( ! identical( TT[[gene]], chol( (t(x.curr[[gene]][, nonzero] )%*% x.curr[[gene]][, nonzero]) + small * diag(qq)))) browser()
			}
			#if(length(nonzero) ==0) browser()
		    }
		   browser() 
		    ##Gibbs sampling for beta and sigma
		    zero <- which(gamMa[[gene]]==0); 
		    #VV =  1/(1/sigma2 + 1/cc)(X'X)^-1,  mm = VV X'Y    
		    if(qq > 0) {
			TT.inv = backsolve(TT[[gene]], diag(qq)) 
			VV <- 1/(1/sigma2[[gene]] + 1/cc[[gene]]) * (TT.inv %*% t(TT.inv) )
			mm <- VV %*% xTy[[gene]][nonzero]/sigma2[[gene]]
			beTa[[gene]][nonzero] <- mvrnorm(n=1, mu=mm, Sigma=VV)
			beTa[[gene]][zero] <- 0
			xbeTa  <-  x.curr[[gene]][, nonzero, drop=F]%*% beTa[[gene]][nonzero]
		    }else{
			beTa[[gene]][zero] <- 0; xbeTa <- 0
		    }
		    if( any(abs(beTa[[gene]]) > 100) ) browser()
		    residue <- y[,gene,drop=F] - xbeTa 
		    sigma2[[gene]] <- rinvgamma(n=1, shape= (qq + nu[[gene]]+ n +1), scale= 1/2 * (sum(residue^2) + nu[[gene]]* lambdal[[gene]]))     
		    cc[[gene]] <- rinvgamma(n=1, shape=(qq+1)/2, scale = (sum(xbeTa^2) + n) /2)

		    zz[gene.start:gene.end] <- gamMa.rate - .5; # TODO change if multiple iteration are used
		    if((iter %%thin==0)  && (iter > burnIn)){
			gamMa.mat[[gene]][count,] <- gamMa[[gene]]
			gamMa.prob.mat[[gene]][count,] <- gamMa.prob
			beTa.mat[[gene]][count,] <- beTa[[gene]]
		    }
		}

	    }
	}
	#browser();
	if(estimate.alpha){
	    if(Debug) browser()
	    #sample wi
	    freq.nonzero <- which(logit.freq > 0); p.nonzero <- length(freq.nonzero);
	    fTalpHa <- featureCat %*% alpHa
	    # subsampling is required here due imbalanced dataset
	    #w[freq.nonzero] <- rpg(num=p.nonzero, h=logit.freq, z=fTalpHa) 
	    w[freq.nonzero] <- sapply(1:p.nonzero, function(pp)  rpg(num=1, h=logit.freq[freq.nonzero[pp]], z=fTalpHa[pp])) # there can be variable -- fix the length
	    ##calculate m.w V.w
	    diag.w  <- diag(w[freq.nonzero])
	    V.w <- ginv((t(featureCat) %*% diag.w %*% featureCat)+  diag( B.inv.alpHa)) 
	    m.w <- V.w %*% t(featureCat) %*% zz # this is like omegafg
	    alpHa[] <- mvrnorm(n=1, mu = m.w, Sigma=V.w)
	}
	if((iter %%thin==0)  && (iter > burnIn)){
	    alpHa.mat[count,] <- alpHa
	    count = count +1;

	}
	if((iter%*% 1200) ==0) save(out, file="out.200gene") 

    }
    list(beTa.mat = beTa.mat, alpHa.mat = alpHa.mat, gamMa.mat = gamMa.mat, gamMa.prob.mat=gamMa.prob.mat)
}
