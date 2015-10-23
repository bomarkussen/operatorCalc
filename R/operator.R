# ------------------
# Class definition
# ------------------

operator <- setClass("operator", slots=c(alpha="matFct",beta="matFct",gamma="matFct",delta="matFct",epsilon="matFct"))

setValidity("operator", function(object) {
  msg <- NULL
  valid <- TRUE
  if (length(unique(c(object@alpha@mesh[1],
                      object@beta@mesh[1],
                      object@gamma@mesh[1],
                      object@delta@mesh[1],
                      object@epsilon@mesh[1])))!=1) {
    valid <- FALSE
    msg <- c(msg,"Left endpoints of the domains of the coefficient functions must coinside")
  }
  if (length(unique(c(object@alpha@mesh[length(object@alpha@mesh)],
                      object@beta@mesh[length(object@beta@mesh)],
                      object@gamma@mesh[length(object@gamma@mesh)],
                      object@delta@mesh[length(object@delta@mesh)],
                      object@epsilon@mesh[length(object@epsilon@mesh)])))!=1) {
    valid <- FALSE
    msg <- c(msg,"Right endpoints of the domains of the coefficient functions must coinside")
  }
  if (dim(object@alpha@f)[1]!=dim(object@alpha@f)[2]) {
    valid <- FALSE
    msg <- c(msg,"Value and domain dimension of the multiplication operator must coinside")
  }
  if (dim(object@beta@f)[1]!=dim(object@alpha@f)[1]) {
    valid <- FALSE
    msg <- c(msg,"Triangular and multiplication operators must have the same value dimension")
  }
  if (dim(object@gamma@f)[1]!=dim(object@alpha@f)[2]) {
    valid <- FALSE
    msg <- c(msg,"Triangular and multiplication operators must have the same domain dimension")
  }
  if (dim(object@beta@f)[2]!=dim(object@gamma@f)[2]) {
    valid <- FALSE
    msg <- c(msg,"K dimension for value and domain of the triangular operator must coinside")
  }
  if (dim(object@delta@f)[1]!=dim(object@alpha@f)[1]) {
    valid <- FALSE
    msg <- c(msg,"Projection and multiplication operators must have the same value dimension")
  }
  if (dim(object@epsilon@f)[1]!=dim(object@alpha@f)[2]) {
    valid <- FALSE
    msg <- c(msg,"Projection and multiplication operators must have the same domain dimension")
  }
  if (dim(object@delta@f)[2]!=dim(object@epsilon@f)[2]) {
    valid <- FALSE
    msg <- c(msg,"L dimension for value and domain of the projection operator must coinside")
  }
  if (valid) TRUE else msg
})

# ---------
# Methods
# ---------

setMethod("show",
          signature="operator",
          definition=function(object) {
            cat("An object of class ",class(object),"\n",sep="")
            cat(" Lattice operator on the space of functions: [",object@alpha@mesh[1],",",
                object@alpha@mesh[length(object@alpha@mesh)],"] -> ",
                dim(object@alpha@f)[1],"-dimensional space\n",sep="")
            cat(" Number of triangular and projection components = (",dim(object@beta@f)[2],
                ",",dim(object@delta@f)[2],")\n",sep="")
            invisible(NULL)
          })

setMethod("superSample",
          signature=c(x="numeric",e1="operator"),
          definition=function(x,e1) {
            operator(alpha=superSample(x,e1@alpha),
                     beta=superSample(x,e1@beta),
                     gamma=superSample(x,e1@gamma),
                     delta=superSample(x,e1@delta),
                     epsilon=superSample(x,e1@epsilon))
          })

setMethod("[",
          signature=c(x="operator",i="numeric"),
          definition=function(x,i) {
            k <- dim(x@beta@f)[2];  k <- min(1,k):k
            l <- dim(x@delta@f)[2]; l <- min(1,l):l
            operator(alpha=x@alpha[i,i],beta=x@beta[i,k],gamma=x@gamma[i,k],delta=x@delta[i,l],epsilon=x@epsilon[i,l])
          })

setMethod("-",
          signature=c(e1="operator"),
          definition=function(e1) {
            operator(alpha=-e1@alpha,
                     beta=-e1@beta,
                     gamma=e1@gamma,
                     delta=-e1@delta,
                     epsilon=e1@epsilon)
          })

setMethod("+",
          signature=c(e1="operator",e2="operator"),
          definition=function(e1,e2) {
            operator(alpha=e1@alpha+e2@alpha,
                     beta=e1@beta|e2@beta,
                     gamma=e1@gamma|e2@gamma,
                     delta=e1@delta|e2@delta,
                     epsilon=e1@epsilon|e2@epsilon)
          })

setMethod("-",
          signature=c(e1="operator",e2="operator"),
          definition=function(e1,e2) {
            operator(alpha=e1@alpha-e2@alpha,
                     beta=e1@beta|-e2@beta,
                     gamma=e1@gamma|e2@gamma,
                     delta=e1@delta|-e2@delta,
                     epsilon=e1@epsilon|e2@epsilon)
          })

setMethod("*",
          signature=c(e1="operator",e2="operator"),
          definition=function(e1,e2) {
            operator(alpha=e1@alpha*e2@alpha,
                     beta=e1@beta|e1@alpha*e2@beta-e1@beta*backward(transpose(e1@gamma)*e2@beta),
                     gamma=transpose(e2@alpha)*e1@gamma+e2@gamma*backward(transpose(e2@beta)*e1@gamma)|e2@gamma,
                     delta=e1@delta|e1@alpha*e2@delta+e1@beta*forward(transpose(e1@gamma)*e2@delta),
                     epsilon=transpose(e2@alpha)*e1@epsilon+e2@gamma*backward(transpose(e2@beta)*e1@epsilon)+e2@epsilon*integral(transpose(e2@delta)*e1@epsilon)|e2@epsilon)
          })

setMethod("*",
          signature=c(e1="operator",e2="matFct"),
          definition=function(e1,e2) {
            e1@alpha*e2+e1@beta*forward(transpose(e1@gamma)*e2)+e1@delta*integral(transpose(e1@epsilon)*e2)
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
            
            # Take properties from matFct object and global options: c++ needs order > 1
            order <- max(2,max(options()$"minimal.polynomial.order",min(options()$"maximal.polynomial.order",dim(e1@alpha@g)[4])))
            super <- 2*order+1
            
            # super sampling according to order of polynomial interpolation
            t.obs <- as.numeric(e1@alpha@mesh)
            t.super <- cumsum(c(t.obs[1],rep(diff(t.obs)/(super+1),each=super+1)))
            
            # ODE for triangular part
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
              kappa.transpose <- aperm(array(deSolve::ode(y=c(diag(K)),times=t.super,func=dy,parms=parms,rtol=1e-10,atol=1e-10)[,-1],
                                             dim=c(length(t.super),K,K)),c(2,3,1))
            }
            
            # solve multiplication part
            inv.alpha <- inverse(e1@alpha)
            
            # solve triangular part
            if (K==0) {
              inv.beta  <- e1@beta
              inv.gamma <- e1@gamma
            } else {
              # approximate tilde(beta)
              res <- .Call('operatorCalc_inv_beta', PACKAGE = 'operatorCalc', e1@alpha@f, e1@alpha@g, e1@beta@f, e1@beta@g, 
                           kappa.transpose, order, super, e1@alpha@continuous&e1@beta@continuous)
              attr(res$g,"dim") <- c(N,K,order,p-1)
              inv.beta <- matFct(mesh=e1@alpha@mesh,f=res$f,g=aperm(res$g,c(1,2,4,3)),continuous=e1@alpha@continuous&e1@beta@continuous)
              # approximate tilde(gamma)
              res <- .Call('operatorCalc_inv_gamma', PACKAGE = 'operatorCalc', e1@alpha@f, e1@alpha@g, e1@gamma@f, e1@gamma@g, 
                           kappa.transpose, order, super, e1@alpha@continuous&e1@gamma@continuous)
              attr(res$g,"dim") <- c(N,K,order,p-1)
              inv.gamma <- matFct(mesh=e1@alpha@mesh,f=res$f,g=aperm(res$g,c(1,2,4,3)),continuous=e1@alpha@continuous&e1@gamma@continuous)
            }
            
            # solve projection part
            if (L==0) {
              inv.epsilon <- e1@epsilon
              inv.delta   <- e1@delta
            } else {
              if (K==0) {
                inv.epsilon <- inverse(transpose(e1@alpha),e1@epsilon)
                inv.delta   <- -inverse(e1@alpha,e1@delta)%*%base::solve(diag(L)+integral(transpose(inv.epsilon)*e1@delta))
              } else {
                inv.epsilon <- inverse(transpose(e1@alpha),e1@epsilon)+inv.gamma*backward(transpose(inv.beta)*e1@epsilon)
                inv.delta   <- -(inverse(e1@alpha,e1@delta)+inv.beta*forward(transpose(inv.gamma)*e1@delta))*
                  base::solve(diag(L)+integral(transpose(inv.epsilon)*e1@delta))                
              }
            }
            
            # return result
            operator(alpha=inv.alpha,beta=inv.beta,gamma=inv.gamma,delta=inv.delta,epsilon=inv.epsilon)
          })
