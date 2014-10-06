logistic <- function(x) 1/(1+exp(-x)) 
fit.logit <- function(yy, xx, B.inv=rep(1/100, dim(xx)[2]),  maxiter = 2000, burnIn = 400, thin = 40)
{
   require(BayesLogit); require(MASS)  
  num <- length(yy); dd <- dim(xx)[2]; wlocal <- rep(0, num);    bb <- rep(0,dd) 
    bb.mat <- matrix(0, nrow=dd, ncol=floor((maxiter - burnIn)/thin))
    count <- 0
    for(iter in seq(maxiter)){
	xxTbb <- xx %*% bb
	# ubsampling is required here due imbalanced dataset
	wlocal[] <- sapply(xxTbb, function(xtb)  rpg(num=1, h=1, z=xtb)) # there can be variable -- fix the length
	#calculate m.w V.w
	diag.w  <- diag(wlocal, num, num)
	#V.w <- ginv((t(xx) %*% diag.w %*% xx)+  diag( B.inv)) # change this to cholesky decomposition 
	L <- chol(t(xx) %*% diag.w %*% xx +  diag( B.inv))
	Lt.inv <- solve(a=t(L), b=diag(1,length(B.inv)))
	V.w <- solve(a=L, b=Lt.inv)
	m.w <- V.w %*% (t(xx) %*% (yy-.5)) # this is like omegafg
	bb <- mvrnorm(n=1, mu = m.w, Sigma=V.w)
	if((iter %%thin==0)  && (iter > burnIn)){
	    bb.mat[,count] <- bb; count <- count + 1;
	}
    }
    #apply(bb.mat,1, median)
    bb.mat
}


simulate.eQTL <- function (object,   numSNP = 512, numGenes = 5, numSamp=20, nsim=1, seed=NULL, interaction=T, sparse=T, ...){
    if (!is.null(seed)) set.seed(seed)
    #   is inefficient having two loops 1:nsim
    #------------------------
    if(exists("object$x")) x = object$x 
    else x <- matrix(rbinom(n=numSNP*numSamp, prob=0.3, size=1), numSNP, numSamp)
    if(exists("object$pp")) pp = object$pp
    else {
	if(interaction)
	    pp <- array(0,c(numSNP, numGenes, 4))
	else 
	    pp <- array(0,c(numSNP, numGenes, 3))
    pp[,,1] = 1
    pp[,,2:3] = array(rbeta(n=numSNP*numGenes*2, shape1=0.2, shape2=0.5), c(numSNP, numGenes, 2))
    if (interaction) pp[,,4] = pp[,,2] * pp[,,3]
    }
    #w = matrix(rbeta(n=numSNP*numGenes, shape1=0.2, shape2=0.5), numSNP, numGenes)
    w = matrix(rnorm(n=numSNP*numGenes, mean=0.8,sd=0.2), numSNP, numGenes)
    if(sparse){
	wnonzero1 = ceiling(min(.1*numSNP, 2))
	wnonzero2= ceiling(min(numSNP, 30))
	wnonzero = sample(x=seq(wnonzero1, wnonzero2), replace=T, size=numGenes)
	wmask = sapply(wnonzero , function(x) 
		       sample(c(-1,1,0), replace=T, prob=c(x/(2*numSNP), x/(2*numSNP), (numSNP-x)/numSNP), size=numSNP))
	w = w*wmask 
    } 
    if(exists("object$alpha")) alpha = object$alpha 
    else {
	if (interaction) alpha = matrix( runif(4), ncol=1) 
	else  alpha = matrix( runif(3), ncol=1)
    } 
    findCC <- function(pp,alpha) apply(pp, c(1,2), function(x) sum(x*alpha))
    cc = findCC(pp,alpha)
    y =  t(cc*w) %*% x + matrix(rnorm(n=numGenes*numSamp, sd=.001),numGenes, numSamp )  
    #x = seq(0,1, length=100) 
    #plot(x, dbeta(x,shape1=0.2, shape2=0.5) )

    object$x <- x
    object$y <- y
    object$w  <- w
    object$pp <- pp
    object$cc  <- cc
    object$alpha  <- alpha
    return(object)
}

generate.eeSNP <- function (object,   numSNP = 512,
			    numGenes = 5, numSamp=20, numFeature=9, nsim=1, seed=NULL, interaction=T, sparse=T, prefix , ev =0.2, logistic=F, includeAlpha=T, ...){
    if (!is.null(seed)) set.seed(seed)
    #   is inefficient having two loops 1:nsim
    #------------------------
    if(!is.null(object$x)) x = object$x
    else x <- matrix(sample(x=c(0,1,2,NA) ,size=numSNP*numSamp, replace=T , prob=c(.6, .25, .15, .002)), numSNP, numSamp)
    feature <- matrix(rnorm(n=numSNP*(numFeature) )^2, nrow=numSNP, ncol=(numFeature)) 
    x.nona = x
    x.nona[is.na(x)]=0
    if(!is.null(object$pp)) pp = object$pp
    else {
	if(interaction)
	    pp <- array(0,c(numSNP, numGenes, 4))
	else 
	    pp <- array(0,c(numSNP, numGenes, 3))
    pp[,,1] = 1
    temp = rbeta(n=numSNP, shape1=0.2, shape2=0.5)
    pp[,,2] = temp
    pp[,,3] = array(rbeta(n=numSNP*numGenes*1, shape1=0.5, shape2=0.2), c(numSNP, numGenes))
    #pp[,,2] = array(rbeta(n=numSNP*numGenes*2, shape1=0.2, shape2=0.5), c(numSNP, numGenes))
    #pp[,,3] = 0
    if (interaction) pp[,,4] = pp[,,2] * pp[,,3]
    }
    w = matrix(rnorm(n=numSNP*numGenes, mean=1.5, sd=1), numSNP, numGenes)
    #w = t( .8*pp[,,2] +  t(matrix(rnorm(n=numSNP*numGenes, mean=0.8,sd=0.2), numSNP, numGenes)))
    if(sparse){
	wnonzero1 = ceiling(min(.1*numSNP, 10))
	wnonzero2= ceiling(min(numSNP, 40))

	wnonzero = sample(x=seq(wnonzero1, wnonzero2), replace=T, size=numGenes)
	#wmask = sapply(wnonzero , function(x) 
#sample(c(-1,1,0), replace=T, prob=c(x/(2*numSNP), x/(2*numSNP), (numSNP-x)/numSNP), size=numSNP))
	if(F){	
	  wmask = sapply(wnonzero , function(x){
			 temp1 = sample(1:numSNP ,size=x, replace=F, prob=pp[,1,2]) 
			 temp2 = rep(0,numSNP)
			 temp2[temp1] = sample(c(1,-1), size=x, replace=T)
			 temp2})   
	}else{
	  wmask = sapply(wnonzero , function(x){
			 temp1 = sample(1:numSNP ,size=x, replace=F) 
			 temp2 = rep(0,numSNP)
			 temp2[temp1] = sample(c(1,-1), size=x, replace=T)
			 temp2})


	}
	w = w*wmask
    }
   if(!logistic) alpha = matrix(runif(3*(numFeature+1)),nrow=3)
   else alpha = matrix(rnorm(3*(numFeature+1), mean=0, sd=.5),nrow=3)
 
    z = x
   if(includeAlpha){ for(ii in 1:numSNP){
     featureCurr = c(1,feature[ii,])
     for(jj in 1:numSamp){
       if(!logistic) z[ii,jj]=sum(alpha[x.nona[ii,jj]+1,]*featureCurr)
       else z[ii,jj]= logistic(sum(alpha[x.nona[ii,jj]+1,]*featureCurr))
     }
    }
   }else
     z <- matrix(rnorm(n=numSNP*numSamp, mean=0, sd=2), numSNP, numSamp)
    
    #if(!is.null(object$alpha)) {
	#alpha = object$alpha
    #} 
    #else {
	#if (interaction) alpha = matrix( runif(4), ncol=1) 
	#else  alpha = matrix( runif(3), ncol=1)
    #} 
    #findCC <- function(pp,alpha) apply(pp, c(1,2), function(x) sum(x*alpha))
    #cc = findCC(pp,alpha)
    ##y =  t(cc*w) %*% x.nona + matrix(rnorm(n=numGenes*numSamp, sd=.01),numGenes, numSamp ) 
    z.nona = z; z.nona[is.na(z)]=0
    y.predict =  t(w) %*% z.nona
    y =t( apply( y.predict, 1,   function(x)  x+ rnorm(n=numSamp, sd=sd(x)*sqrt((1-ev)/ev) )))
    alpha = cbind(alpha, rep(0,3))
   
			  
    #x = seq(0,1, length=100) 
    #plot(x, dbeta(x,shape1=0.2, shape2=0.5) )

    #writing to file
    pedfile=paste(prefix, ".raw", sep="")
    genefile=paste(prefix, ".gene", sep="")
    enhfile=paste(prefix, ".enhancer", sep="")
    dimfile=paste(prefix, ".dim", sep="")
    dhsfile=paste(prefix, ".dhs", sep="")
    rnafile=paste(prefix, ".rna", sep="")
    eqtlfile=paste(prefix, ".eqtl", sep="")
    featurefile=paste(prefix, ".feature", sep="")
    ypredictfile=paste(prefix, ".ypredict", sep="")
    alphafile=paste(prefix, ".alpha", sep="")
    truealphafile=paste(prefix, ".truealpha", sep="")
    truezfile=paste(prefix, ".truez", sep="")
    print("here")
    pedinfo = data.frame(id = 1:numSamp, a=1, b=0, c=0, d=2, e=2)
    temp = data.frame(t(x))
    ped = cbind(pedinfo,temp)
    write.table(file=pedfile,x = ped , row.names = F, col.names =T,  sep=" ", quote=F )
    write.table(file=genefile,x =y  , row.names = T, col.names =T,  sep="\t", quote=F )
    write.table(file=enhfile,x =data.frame( b=pp[,1,2], c=0) ,
		row.names = F, col.names =T,  sep=" ", quote=F )
    write.table(file=dimfile,x =data.frame(a=numSNP, b=numSamp, c=numGenes),
		row.names = F, col.names =F,  sep=",", quote=F )
   info=data.frame(a="chr1", b=100006705,c=100007205, d="SNP_A-2257965")
  dhsmatrix = data.frame(matrix(rnorm(numSNP*15), ncol=15, nrow=numSNP))
    write.table(file=dhsfile,x= cbind(info, dhsmatrix),
		row.names = F, col.names =F,  sep=",", quote=F )
   info=data.frame(a="chr1", b=100006705,c=100007205,e="7896744", d="ENGL2257965")
  rnamatrix = data.frame(matrix(rnorm(numGenes*15), ncol=15, nrow=numGenes))
    write.table(file=rnafile,x= cbind(info, rnamatrix),
		row.names = F, col.names =F,  sep=",", quote=F )
   info=data.frame(a="chr1", b=100006705,c=100007205,e="SNP_A01909444", d=1:numSNP)
    write.table(file=featurefile,x= cbind(info, feature),
		row.names = F, col.names =F,  sep="\t", quote=F )
    require(Matrix)
    writeMM(file=eqtlfile,obj=Matrix(w ,sparse=T))
    write.table(file=ypredictfile,x = y.predict , row.names = F, 
    	col.names =F,  sep="\t", quote=F )
    if(!logistic) temp = cbind(c(0,1,2), matrix(0,ncol=10,nrow=3)) 
    else  temp = cbind(c(-6,0,6), matrix(0,ncol=10,nrow=3));  
    write.table(file=alphafile,x = temp , row.names = F, 
    	col.names =F,  sep=",", quote=F )
    write.table(file=truealphafile,x = alpha, row.names = F, 
    	col.names =F,  sep=",", quote=F )
    write.table(file=truezfile,x = z, row.names = F, 
    	col.names =F,  sep=",", quote=F )

    object$x <- x
    object$y <- y
    object$y.predict <- y.predict
    object$w  <- w
    object$z <- z 
    #object$pp <- pp
    #object$cc  <- cc
    object$alpha  <- alpha
    object$feature <- cbind(1, feature, 0)
    return(object)
}

generate.eeSNP.model5 <- function ( numSNP = 512, avg.eSNP=20,
				   numGenes = 1, numSamp=20, numFeature=9, seed=NULL,  prefix, ev =0.2,  ...){
    if (!is.null(seed)) set.seed(seed)
    #x <- matrix(sample(x=c(0,1,2,NA) ,size=numSNP*numSamp, replace=T , prob=c(.6, .25, .15, .002)), ncol=numSNP, nrow=numSamp)
    x <- matrix(rnorm(n=numSNP*numSamp), ncol=numSNP, nrow=numSamp)
    #feature <- matrix(rexp(n=numSNP*(numFeature), rate=1), nrow=numSNP, ncol=(numFeature)) 
    feature <- matrix(rnorm(n=numSNP*(numFeature), sd=.3), nrow=numSNP, ncol=(numFeature)) 
    x.nona = x
    x.nona[is.na(x)]=0
    # alpHa N(0,1) 
    #alpHa0 = log(avg.eSNP/(numSNP- avg.eSNP))
    alpHa = rnorm(numFeature, sd=4, mean=0)
    const <- feature %*% matrix( alpHa, ncol=1)
    alpHa0 <- quantile(const,probs= (1 - avg.eSNP/numSNP)) 
    const1 <- const - alpHa0
    alpHa <- c(-alpHa0, alpHa)
    gamMa.prob <- logistic(const1)
    #gamMa.prob <- sapply(1:numSNP, function(ii) logistic(sum(alpHa*c(1, feature[ii,]))) )
    gamMa <- sapply(gamMa.prob, function(prb) sample(x=0:1, size=1, prob=c(1-prb, prb)) )
    #gamMa <- rep(0, numSNP)
   #gamMa[gamMa.prob > .5] <- 1 
    zero = which(gamMa==0); nonzero = which(gamMa!=0)
    beTa <- rep(0, numSNP)
    beTa[nonzero] <- rnorm(n=length(nonzero)) # change it later to take higher values for gamMa.prob
    y.predict <- x.nona %*% beTa 
    y = y.predict + rnorm(n=numSamp, sd=sd(y.predict)*sqrt((1-ev)/ev) )
    object = c(list(x=x, y=y, y.predict=y.predict, feature=feature, alpHa=alpHa, gamMa=gamMa, gamMa.prob=gamMa.prob, beTa=beTa))
    class(object) = "eeSNP.model5"
    return(object)
}


generate.eeSNP.logistic <- function ( numSNP = 512, avg.eSNP=20,
				   numGenes = 1, numSamp=20, numFeature=9, seed=NULL,  prefix, ev =0.2, LD=FALSE, ...){
    if (!is.null(seed)) set.seed(seed)
if(!LD){    
  x <- matrix(sample(x=c(0,1,2,NA) ,size=numSNP*numSamp, replace=T , prob=c(.6, .25, .15, .002)), ncol=numSNP, nrow=numSamp)
    #x <- matrix(rnorm(n=numSNP*numSamp), ncol=numSNP, nrow=numSamp)
   enhancer <- sample.int(n=numSNP, size=avg.eSNP) # only one gene for now
}else{
  indep.SNP <- ceiling(numSNP/2)
  x.indep <- matrix(sample(x=c(0,1,2,NA) ,size= indep.SNP*numSamp, replace=T , prob=c(.6, .25, .15, .002)), ncol=indep.SNP, nrow=numSamp)
  x <- matrix( rbind(x.indep,x.indep), nrow=numSamp)[,1:numSNP]  
   enhancer <- sample(seq(1, numSNP, 2), size=avg.eSNP)
}

   nonenhancer <- which( !(seq(numSNP) %in% enhancer))
   feature = matrix(0, nrow=numSNP, ncol=(numFeature))
   feature[nonenhancer, ] <- rnorm(n=length(nonenhancer)*(numFeature), sd=1)  
   feature[enhancer, ] <- rnorm(n=length(enhancer)*(numFeature), mean=.5, sd=1)
    gamMa.ini <- rep(0, numSNP)
gamMa.ini[enhancer] <- 1
   #alpHa.mat <- fit.logit(yy=gamMa.ini, xx=cbind(1, feature), maxiter=10000, burnIn=1000, thin=100)
   temp.fit <- BayesLogit::logit(y=gamMa.ini, X=cbind(1, feature), P0 = 1/100*diag(numFeature +1) )
 
#alpHa <- apply(alpHa.mat,1, median)
alpHa <- apply(temp.fit$beta,2, median)
   x.nona = x
    x.nona[is.na(x)]=0
    # alpHa N(0,1) 
    #alpHa0 = log(avg.eSNP/(numSNP- avg.eSNP))
    #alpHa = rnorm(numFeature, sd=4, mean=0)
    const <- cbind(1, feature) %*% matrix( alpHa, ncol=1)
    #alpHa0 <- quantile(const,probs= (1 - avg.eSNP/numSNP)) 
    #const1 <- const - alpHa0
    #alpHa <- c(-alpHa0, alpHa)
    gamMa.prob <- logistic(const)
    #gamMa.prob <- sapply(1:numSNP, function(ii) logistic(sum(alpHa*c(1, feature[ii,]))) )
    #gamMa <- sapply(gamMa.prob, function(prb) sample(x=0:1, size=1, prob=c(1-prb, prb)) )
    gamMa <- numeric(numSNP)
   gamMa[gamMa.prob > .5] <- 1 
    zero = which(gamMa==0); nonzero = which(gamMa!=0)
    beTa <- rep(0, numSNP)
    beTa[nonzero] <- rnorm(n=length(nonzero)) # change it later to take higher values for gamMa.prob
    y.predict <- x.nona %*% beTa 
    y = y.predict + rnorm(n=numSamp, sd=sd(y.predict)*sqrt((1-ev)/ev) )
    object = c(list(x=x, y=y, y.predict=y.predict, feature=feature, alpHa=alpHa, gamMa=gamMa, gamMa.prob=gamMa.prob, beTa=beTa))
    class(object) = "eeSNP.model5"
    return(object)
}


enhancer.selection <- function ( numSNP = 512, avg.eSNP=20,
				numGenes = 1, numSamp=20, numFeature=9, seed=NULL,  prefix, ev =0.2, LD=FALSE, ...){
    if (!is.null(seed)) set.seed(seed)
    if(!LD){    
	x <- matrix(sample(x=c(0,1,2,NA) ,size=numSNP*numSamp, replace=T , prob=c(.6, .25, .15, .002)), ncol=numSNP, nrow=numSamp)
	#x <- matrix(rnorm(n=numSNP*numSamp), ncol=numSNP, nrow=numSamp)
	#enhancer <- sample.int(n=numSNP, size=avg.eSNP) # only one gene for now
	#gamMa <- sample(c(-1,1), size=avg.eSNP, replace=T)
	gamMa <- numeric(numSNP)
	enhancer <- sample(seq(1, numSNP, 2), size=avg.eSNP)
	effect <- sample(c(1,-1), size=length(enhancer), replace=T)
	gamMa[enhancer] <- effect  
    }else{
	indep.SNP <- ceiling(numSNP/2)
	x.indep <- matrix(sample(x=c(0,1,2,NA) ,size= indep.SNP*numSamp, replace=T , prob=c(.6, .25, .15, .002)), ncol=indep.SNP, nrow=numSamp)
	x <- matrix( rbind(x.indep,x.indep), nrow=numSamp)[,1:numSNP]  
	gamMa <- numeric(numSNP)
	enhancer <- sample(seq(1, numSNP, 2), size=avg.eSNP)
	effect <- sample(c(1,-1), size=length(enhancer), replace=T)
	gamMa[enhancer] <- effect  
    }
    beTa <- matrix(0, 3, numFeature + 1, 3);
    beTa[1,] = rnorm(numFeature + 1, sd=1);  
    beTa[2,] = rnorm(numFeature + 1, sd=2);  
    beTa[3,] = rnorm(numFeature + 1, sd=3);  
    nonenhancer <- which( !(seq(numSNP) %in% enhancer))
    feature = matrix(0, nrow=numSNP, ncol=(numFeature))
    feature[nonenhancer, ] <- rnorm(n=length(nonenhancer)*(numFeature), sd=1)  
    feature[enhancer, ] <- rnorm(n=length(enhancer)*(numFeature), mean=3, sd=1)
    x.nona = x
    x.nona[is.na(x)]=0
    
    beTaF  <- beTa %*% t(cbind(1,feature))  
    beTaxF <- sapply(1:numSNP, function(ss) beTaF[ x.nona[,ss]+1,ss]) 
    y.predict <-  beTaxF %*% gamMa 
    y = y.predict + rnorm(n=numSamp, sd=sd(y.predict)*sqrt((1-ev)/ev) )
    object = c(list(x=x, y=y, y.predict=y.predict, feature=feature, gamMa=gamMa, beTa=beTa))
    class(object) = "eeSNP.model5"
    return(object)
}

