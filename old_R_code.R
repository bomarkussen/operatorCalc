setMethod("inverse",
          signature=c(e1="matFct",e2="missing"),
          definition=function(e1) {
            # Take properties from matFct object
            order <- dim(e1@g)[4]
            ifelse(order>1, super <- 2*order+1, super <- 0)
            continuous <- e1@continuous
            if (length(continuous)!=dim(e1@g)[3]) continuous <- rep(continuous[1],dim(e1@g)[3])
            # Call Cpp
            res <- .Call('operatorCalc_inverse', PACKAGE = 'operatorCalc', e1@f, e1@g, order, super, continuous)
            # fix Rcpp result: piecewise constant function?
            if (order==0) {
              res$g <- array(res$g[,1:dim(e1@g)[3]],dim=c(dim(e1@g)[1:2],1,dim(e1@g)[3]))
              continuous <- rep(FALSE,dim(e1@g)[3])
            }
            # fix Rcpp result: piecewise linear function?
            if (order==1) {
              res$g <- array(0,dim=c(dim(e1@g)[1:2],0,dim(e1@g)[3]))
              continuous <- TRUE
            }
            # fix Rcpp result: piecewise polynomial function?
            if (order>1) {
              attr(res$g,"dim") <- c(dim(e1@g)[1:2],order,dim(e1@g)[3])
              continuous <- e1@continuous
            }
            # return result
            return(matFct(mesh=e1@mesh,f=res$f,g=aperm(res$g,c(1,2,4,3)),continuous=continuous))
          })

setMethod("inverse",
          signature=c(e1="matFct", e2="matFct"),
          definition=function(e1, e2) {
            # sanity check
            if (dim(e1@f)[1]!=dim(e1@f)[2]) stop("first matrix must be square")
            if (dim(e1@f)[2]!=dim(e2@f)[1]) stop("mismatch between matrix dimensions")
            # if e1 and e2 have different meshes then super sample both functions
            if (!setequal(e1@mesh,e2@mesh)) {
              e1 <- superSample(e2@mesh,e1)
              e2 <- superSample(e1@mesh,e2)
            }
            # find order: to be sophisticated lated!?
            order <- max(2,dim(e1@g)[4]+dim(e2@g)[4])
            super <- 2*order+1
            # find continuity property
            if (length(e1@continuous)==1) {
              if (length(e2@continuous)==1) {
                res.cont <- e1@continuous & e2@continuous
                continuous <- rep(res.cont,dim(e1@g)[3])
              } else {
                res.cont <- rep(e1@continuous,dim(e1@g)[3]) & e2@continuous
                continuous <- res.cont
              }
            } else {
              if (length(e2@continuous)==1) {
                res.cont <- e1@continuous & rep(e2@continuous,dim(e2@g)[3])
                continuous <- res.cont
              } else {
                res.cont <- e1@continuous & e2@continuous
                continuous <- res.cont
              }
            }
            # pointwise solution of linear equation
            res <- .Call('operatorCalc_solve', PACKAGE = 'operatorCalc', e1@f, e1@g, e2@f, e2@g, order, super, continuous)
            attr(res$g,"dim") <- c(dim(e1@g)[1],dim(e2@g)[2],order,dim(e1@g)[3])
            # return result
            return(matFct(mesh=e1@mesh,f=res$f,g=aperm(res$g,c(1,2,4,3)),continuous=continuous))
          })

setMethod("inverse",
          signature=c(e1="operator"),
          # update of old implementation: so far the best! worst error 1e-6
          # but should be improved on accuracy and speed!
          definition=function(e1) {
            # take dimensions
            p <- length(e1@alpha@mesh)
            N <- dim(e1@alpha@f)[1]
            K <- dim(e1@beta@f)[2]
            L <- dim(e1@delta@f)[2]
            
            # Take properties from matFct object and global options
            order <- max(options()$"minimal.polynomial.order",min(options()$"maximal.polynomial.order",dim(e1@alpha@g)[4]))
            
            # super sampling according to order of polynomial interpolation
            t.obs <- as.numeric(e1@alpha@mesh)
            if (order < 2) {
              t.super <- t.obs
              index   <- 1:length(t.obs)
            } else {
              super   <- 2*order+1
              t.super <- cumsum(c(t.obs[1],rep(diff(t.obs)/(super+1),each=super+1)))
              index   <- seq(1,length(t.super),super+1)
            }
            
            # ODE's for triangular part
            if (K>0) {
              # coefficient function for ODE
              dy <- function(tt,state,parms) {
                with(as.list(parms),{list(c(matrix(evaluate(tt,aa),K,N)%*%
                                              base::solve(matrix(evaluate(tt,bb),N,N),
                                                          matrix(evaluate(tt,cc),N,K))%*%
                                              matrix(state,K,K)))})
              }
              # solve ODE for kappa.transpose
              parms <- list(N=N,K=K,aa=-transpose(e1@gamma),bb=e1@alpha,cc=e1@beta)
              kappa.transpose <- approximate(t.super,
                                             array(deSolve::ode(y=c(diag(K)),times=t.super,func=dy,parms=parms,rtol=1e-10,atol=1e-10)[,-1],
                                                   dim=c(length(t.super),K,K)))
              # solve ODE for kappa.inverse
              parms <- list(N=N,K=K,aa=transpose(e1@beta),bb=transpose(e1@alpha),cc=e1@gamma)
              kappa.inverse <- approximate(t.super,
                                           array(deSolve::ode(y=c(diag(K)),times=t.super,func=dy,parms=parms,rtol=1e-10,atol=1e-10)[,-1],
                                                 dim=c(length(t.super),K,K)))
            }
            
            # solve multiplication part
            inv.alpha <- inverse(superSample(t.super,e1@alpha))
            
            # solve triangular part
            if (K==0) {
              inv.beta  <- e1@beta
              inv.gamma <- e1@gamma
            } else {
              inv.beta  <- -inv.alpha*e1@beta*kappa.transpose
              inv.gamma <- transpose(inv.alpha)*e1@gamma*kappa.inverse
            }
            
            # solve projection part
            if ((K==0) || (L==0)) {
              inv.delta    <- e1@delta
              if (L==0) inv.epsilon <- e1@epsilon else {
                epsilon.star <- transpose(inv.alpha)*e1@epsilon
                inv.epsilon  <- -epsilon.star*base::solve(diag(L)+integral(transpose(e1@delta)*epsilon.star))
              }
            } else {
              inv.epsilon <- transpose(inv.alpha)*e1@epsilon+inv.gamma*backward(transpose(inv.beta)*e1@epsilon)
              inv.delta   <- -(inv.alpha*e1@delta+inv.beta*forward(transpose(inv.gamma)*e1@delta))*
                base::solve(diag(L)+integral(transpose(e1@delta)*inv.epsilon))
            }
            
            # return result
            operator(alpha=inv.alpha,beta=inv.beta,gamma=inv.gamma,delta=inv.delta,epsilon=inv.epsilon)
          })
