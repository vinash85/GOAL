fit.logit <- function(yy, xx, B.inv=rep(1/100, dim(xx)[2]),  maxiter = 2000, burnIn = 400, thin = 40)
{
  require(BayesLogit); require(MASS)  
  num <- length(yy); dd <- dim(xx)[2]; wlocal <- rep(0, num);    bb <- rep(0,dd) 
  bb.mat <- matrix(0, nrow=dd, ncol=floor((maxiter - burnIn)/thin))
  count <- 0
  for(iter in seq(maxiter)){
    xxTbb <- c(xx %*% bb)
    # ubsampling is required here due imbalanced dataset
    #wlocal[] <- sapply(xxTbb, function(xtb)  rpg(num=1, h=1, z=xtb)) # there can be variable -- fix the length
    wlocal[] <- rpg(num=length(xxTbb), h=rep(1, length(xxTbb), z=xxTbb)) 
      #sapply(xxTbb, function(xtb)  rpg(num=1, h=1, z=xtb)) # there can be variable -- fix the length
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

