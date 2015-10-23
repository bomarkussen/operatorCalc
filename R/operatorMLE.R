operator.MLE <- function(matData,components=2,lambda.max=0.1,step.factor=3/2,gamma.min=0,gamma.max=0.1,verbose=TRUE) {
  # matData: (n,p,N)-array with data points.
  # assume equidistant sampling on [0,1], i.e. t_j = (j-1)/(p-1) for j=1,...,p
  
  # simulate data for test: TO BE REMOVED LATER!!!
  n <- 100
  p <- 10
  N <- 2
  matData <- aperm(apply(array(rnorm(n*p*N)/sqrt(p),dim=c(n,p,N)),c(1,3),cumsum),c(2,1,3))+
    rnorm(n*p*N)*1
  matplot(t(matData[,,1]),type="l")
  components <- 2
  
  # take dimensions
  n <- dim(matData)[1]
  p <- dim(matData)[2]
  N <- dim(matData)[3]
  
  # embed data
  X <- approximate(seq(0,1,length.out=p),aperm(matData,c(2,3,1)))
  #plot(transpose(tilde.X))
  
  # remove mean structure and make singular value decomposition
  my.svd <- svd(matrix(matData-rep(apply(matData,c(2,3),mean),each=n),n,p*N),0,components)
  my.data <- array(my.svd$v%*%diag(my.svd$d[1:components],rep(components,2)),
                   dim=c(p,N,components))
  delta1 <- approximate(seq(0,1,length.out=p),my.data*sqrt(p/n))
  alpha1 <- approximate(seq(0,1,length.out=p),
                       aperm(array(apply(matData,2,cov),dim=c(N,N,p)),c(3,1,2))
                       )-(1/p)*delta1*transpose(delta1)
  epsilon1 <- delta1
  
  alpha2 <- inverse(alpha1)
  delta2 <- -inverse(alpha1,delta1)*solve(diag(dim(delta1@f)[2])+integral(transpose(epsilon1)*inverse(alpha1,delta1)))
  epsilon2 <- inverse(transpose(alpha1),epsilon1)
  delta12   <- alpha1*delta2 + epsilon1*forward(transpose(delta1)*delta2) -
    delta1*forward(transpose(epsilon1)*delta2)
  epsilon12 <- transpose(alpha2)*epsilon1 + delta2*backward(transpose(epsilon2)*epsilon1) + 
    epsilon2*forward(transpose(delta2)*epsilon1)
  
  #####################################
  
  plot(alpha1)
  plot(delta1)
  plot(delta1*transpose(epsilon1))
  plot(delta1*transpose(epsilon12)+delta12*transpose(epsilon2))
  
  # true values
  #alpha <- approximate(seq(0,1,length.out=p),aperm(array(c(1,0,0,1),dim=c(N,N,p)),c(3,1,2)))
  #epsilon <- approximate(seq(0,1,length.out=p),aperm(array(c(1,0,0,1),dim=c(N,N,p)),c(3,1,2)))
  #delta   <- approximate(seq(0,1,length.out=p),array(c(seq(0,1,length.out=p),seq(0,0,length.out=p),
  #                                                     seq(0,0,length.out=p),seq(0,1,length.out=p)),dim=c(p,N,N)))
  #summary(t(apply(evaluate(seq(0,1,0.01),alpha),3,function(x) eigen(x)$values)))

  neg.gradient <- function(x) {
    # assume variables: n, X, S, h (list with elements alpha, delta, epsilon)
    tmpS <- S
    if (x!=0) {
      tmpS@alpha@f   <- S@alpha@f+x*h$alpha
      tmpS@delta@f   <- S@delta@f+x*h$delta
      tmpS@epsilon@f <- S@epsilon@f+x*h$epsilon
    }
    invS <- inverse(tmpS)
    return(list(alpha=0*triangular(p*invS@alpha+(1/2)*(invS@delta*transpose(invS@epsilon)+invS@epsilon*transpose(invS@delta))-(p/n)*X*transpose(X)),
                delta=triangular(invS@alpha*tmpS@epsilon+2*invS@delta*backward(transpose(invS@epsilon)*tmpS@epsilon)-(2*p/n)*X*backward(transpose(X)*tmpS@epsilon)),
                epsilon=triangular(invS@alpha*tmpS@delta+2*invS@epsilon*forward(transpose(invS@delta)*tmpS@delta)-(2*p/n)*X*forward(transpose(X)*tmpS@delta))
    ))
  }
  
  S <- inverse(as.operator(symmOperator(alpha=alpha,delta=(sqrt(p))*delta,epsilon=(sqrt(p))*epsilon)))
  
  lambda.start <- 1/(3*N*p*components)
  lambda.max <- 1e-4
  
  maxIter <- 100
  alphaConv <- matrix(0,maxIter,5)
  deltaConv <- matrix(0,maxIter,5)
  epsilonConv <- matrix(0,maxIter,5)
  
  g <- neg.gradient(0)
  h <- g
  #lambda <- lambda.start
  for (iter in 1:maxIter) {
    lambda <- lambda.max
#    for (up.steps in 0:2) {
#      g.new <- neg.gradient(lambda)
#      if (sum(g.new$alpha*h$alpha) +
#          sum(g.new$delta*h$delta) +
#          sum(g.new$epsilon*h$epsilon) < 0) break
#      lambda <- lambda*step.factor
#    }
    for (down.steps in 0:9) {
#      lambda <- lambda*(step.factor+1)/(2*step.factor)
      g.new <- neg.gradient(lambda)
      if (sum(g.new$alpha*h$alpha) +
          sum(g.new$delta*h$delta) +
          sum(g.new$epsilon*h$epsilon) >= 0) break
      lambda <- lambda/step.factor
    }
    
    # update coefficients
    S@alpha@f   <- S@alpha@f   + lambda*h$alpha
    S@delta@f   <- S@delta@f   + lambda*h$delta
    S@epsilon@f <- S@epsilon@f + lambda*h$epsilon
    
    # conjugate gradients
    gamma <- (sum((g.new$alpha-g$alpha)*g.new$alpha) +
              sum((g.new$delta-g$delta)*g.new$delta) +
              sum((g.new$epsilon-g$epsilon)*g.new$epsilon)) /
             (sum(g$alpha*g$alpha)+sum(g$delta*g$delta)+sum(g$epsilon*g$epsilon))
    gamma <- min(gamma.max,max(gamma.min,gamma))
    h$alpha   <- g.new$alpha   + gamma*h$alpha
    h$delta   <- g.new$delta   + gamma*h$delta
    h$epsilon <- g.new$epsilon + gamma*h$epsilon
    g <- g.new
    
    # stepest descent
    #    h <- g.new

    # delta and epsilon of same size: rescale
    #delta.vs.epsilon <- sqrt(sqrt(sum(S@delta@f^2)/sum(S@epsilon@f^2)))
#    delta.vs.epsilon <- sqrt(sqrt(sum(g$delta^2)/sum(g$epsilon^2)))
#    S@delta@f   <- S@delta@f   / delta.vs.epsilon
#    S@epsilon@f <- S@epsilon@f * delta.vs.epsilon
#    h$delta     <- h$delta     / delta.vs.epsilon
#    h$epsilon   <- h$epsilon   * delta.vs.epsilon
#    g$delta     <- g$delta     / delta.vs.epsilon
#    g$epsilon   <- g$epsilon   * delta.vs.epsilon
    
#    invS <- inverse(S)
#    dAlpha <- triangular(-n*1*invS@alpha-n/2*(invS@delta*transpose(invS@epsilon)+invS@epsilon*transpose(invS@delta))+1*X*transpose(X))
#    dDelta <- triangular(-n*invS@alpha*S@epsilon-2*n*invS@delta*backward(transpose(invS@epsilon)*S@epsilon)+2*1*X*backward(transpose(X)*S@epsilon))
#    dEpsilon <- triangular(-n*invS@alpha*S@delta-2*n*invS@epsilon*forward(transpose(invS@delta)*S@delta)+2*1*X*forward(transpose(X)*S@delta))

    # diagnostics
    if (verbose) {
      cat("iteration",iter,": steps down=",down.steps,": lambda=",lambda,", conjugate gamma=",gamma,"\n")
      delta.vs.epsilon <- sqrt(sqrt(sum(g$delta^2)/sum(g$epsilon^2)))
      tmp <- rbind(quantile(g$alpha),quantile(h$alpha),
                   quantile(g$delta / delta.vs.epsilon),quantile(h$delta / delta.vs.epsilon),
                   quantile(g$epsilon * delta.vs.epsilon),quantile(h$epsilon * delta.vs.epsilon))
      rownames(tmp) <- c("g.alpha","h.alpha","g.delta","h.delta","g.epsilon","h.epsilon")
      if (gamma!=0) print(tmp) else print(tmp[c(1,3,5),])
      alphaConv[iter,] <- quantile(g$alpha)
      deltaConv[iter,] <- quantile(g$delta / delta.vs.epsilon)
      epsilonConv[iter,] <- quantile(g$epsilon * delta.vs.epsilon)
    }
  }

  # convergence diagnostics
  if (verbose) {
    par(mfrow=c(1,1))
    matplot(alphaConv[,],type="b",xlab="iteration",ylab="5-point statistics of gradient",main="convergence statistics for alpha")
    matplot(deltaConv[,],type="b",xlab="iteration",ylab="5-point statistics of gradient",main="convergence statistics for delta")
    matplot(epsilonConv[,],type="b",xlab="iteration",ylab="5-point statistics of gradient",main="convergence statistics for epsilon")
  }

  # diagnostics
  invS <- inverse(S)
  plot(invS@alpha); title("inverse: alpha")
  plot(invS@delta*transpose(invS@epsilon)); title("inverse: delta*transpose(epsilon)")
  plot(approximate(seq(0,1,length.out=p),aperm(g$alpha,c(3,1,2)))); title("dAlpha")
  plot(approximate(seq(0,1,length.out=p),aperm(g$delta,c(3,1,2)))); title("dDelta")
  plot(approximate(seq(0,1,length.out=p),aperm(g$epsilon,c(3,1,2)))); title("dEpsilon")

  # image
  im <- matrix(0,p,p)
  for (x in 1:p) for (y in 1:p) {
    im[x,y] <- (1/p)*(evaluate(max(x,y)/p,invS@epsilon)[,,1] %*% evaluate(min(x,y)/p,transpose(invS@delta))[,,1])[1,1]-
      cov(matData[,x,1],matData[,y,1])+(x==y)
  }
  im.start <- matrix(0,p,p)
  for (x in 1:p) for (y in 1:p) {
    im.start[x,y] <- (evaluate(x/p,delta)[,,1] %*% evaluate(y/p,transpose(delta))[,,1])[1,1]-
      cov(matData[,x,1],matData[,y,1])+(x==y)
  }

  image(im,main="lattice misfit")
  summary(c(im))
  image(im.start,main="pca misfit")
  summary(c(im.start))
  par(mfrow=c(1,1))
  plot(sort(c(im.start)),sort(c(im)))
  abline(0,1)
}