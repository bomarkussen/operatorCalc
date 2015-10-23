# ------------------
# Class definition
# ------------------

symmOperator <- setClass("symmOperator", slots=c(alpha="matFct",delta="matFct",epsilon="matFct"))

setValidity("symmOperator", function(object) {
  msg <- NULL
  valid <- TRUE
  if (length(unique(c(object@alpha@mesh[1],
                      object@delta@mesh[1],
                      object@epsilon@mesh[1])))!=1) {
    valid <- FALSE
    msg <- c(msg,"Left endpoints of the domains of the coefficient functions must coinside")
  }
  if (length(unique(c(object@alpha@mesh[length(object@alpha@mesh)],
                      object@delta@mesh[length(object@delta@mesh)],
                      object@epsilon@mesh[length(object@epsilon@mesh)])))!=1) {
    valid <- FALSE
    msg <- c(msg,"Right endpoints of the domains of the coefficient functions must coinside")
  }
  if (dim(object@alpha@f)[1]!=dim(object@alpha@f)[2]) {
    valid <- FALSE
    msg <- c(msg,"Value and domain dimension of the multiplication operator must coinside")
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
          signature="symmOperator",
          definition=function(object) {
            cat("An object of class ",class(object),"\n",sep="")
            cat(" Symmetric lattice operator on the space of functions: [",object@alpha@mesh[1],",",
                object@alpha@mesh[length(object@alpha@mesh)],"] -> ",
                dim(object@alpha@f)[1],"-dimensional space\n",sep="")
            cat(" Number of components = ",dim(object@delta@f)[2],"\n",sep="")
            invisible(NULL)
          })

setMethod("superSample",
          signature=c(x="numeric",e1="symmOperator"),
          definition=function(x,e1) {
            symmOperator(alpha=superSample(x,e1@alpha),
                         delta=superSample(x,e1@delta),
                         epsilon=superSample(x,e1@epsilon))
          })

setMethod("[",
          signature=c(x="symmOperator",i="numeric"),
          definition=function(x,i) {
            l <- dim(x@delta@f)[2]; l <- min(1,l):l
            symmOperator(alpha=x@alpha[i,i],
                         delta=x@delta[i,l],
                         epsilon=x@epsilon[i,l])
          })

setMethod("-",
          signature=c(e1="symmOperator"),
          definition=function(e1) {
            symmOperator(alpha=-e1@alpha,
                         delta=-e1@delta,
                         epsilon=e1@epsilon)
          })

setMethod("+",
          signature=c(e1="symmOperator",e2="symmOperator"),
          definition=function(e1,e2) {
            symmOperator(alpha=e1@alpha+e2@alpha,
                         delta=e1@delta|e2@delta,
                         epsilon=e1@epsilon|e2@epsilon)
          })

setMethod("-",
          signature=c(e1="symmOperator",e2="symmOperator"),
          definition=function(e1,e2) {
            symmOperator(alpha=e1@alpha-e2@alpha,
                         delta=e1@delta|-e2@delta,
                         epsilon=e1@epsilon|e2@epsilon)
          })

setMethod("*",
          signature=c(e1="matrix",e2="symmOperator"),
          definition=function(e1,e2) {
            symmOperator(alpha=e1*e2@alpha,delta=e1*e2@delta,epsilon=e2@epsilon)
          })

setMethod("*",
          signature=c(e1="numeric",e2="symmOperator"),
          definition=function(e1,e2) {
            symmOperator(alpha=e1*e2@alpha,delta=(sign(e1)*sqrt(abs(e1)))*e2@delta,epsilon=sqrt(abs(e1))*e2@epsilon)
          })

# setMethod("inverse",
#           signature=c(e1="symmOperator"),
#           # Using symmetric structure it should be possible to improve
#           # speed and numerical precision
#           definition=function(e1) {
#             tmp <- inverse(as.operator(e1))
#             symmOperator(alpha=tmp@alpha,
#                          delta=tmp@delta,
#                          epsilon=tmp@epsilon)
#           })

setMethod("inverse",
          signature=c(e1="symmOperator"),
          # Use numerical minimization of misfit energy
          definition=function(e1) {
            # Brownian motion
            e1 <- symmOperator(
              alpha=superSample(seq(0,1,0.01),matFct(mesh=mesh(c(0,1)),f=array(c(1,1),dim=c(1,1,2)),g=array(0,dim=c(1,1,1,0)),continuous=TRUE)),
              delta=superSample(seq(0,1,0.01),matFct(mesh=mesh(c(0,1)),f=array(c(0,1),dim=c(1,1,2)),g=array(0,dim=c(1,1,1,0)),continuous=TRUE)),
              epsilon=superSample(seq(0,1,0.01),matFct(mesh=mesh(c(0,1)),f=array(c(1,1),dim=c(1,1,2)),g=array(0,dim=c(1,1,1,0)),continuous=TRUE))
            )

            # make e1 piecewise linear
            e1@alpha@g   <- array(0,dim=c(dim(e1@alpha@g)[1:3],0))
            e1@delta@g   <- array(0,dim=c(dim(e1@delta@g)[1:3],0))
            e1@epsilon@g <- array(0,dim=c(dim(e1@epsilon@g)[1:3],0))
            
            # make initial guess at inverse operator using only projection part
            tmp <- solve(diag(dim(e1@delta@f)[2])+integral(transpose(e1@epsilon)*inverse(e1@alpha,e1@delta)))
            alpha1    <- inverse(e1@alpha)
            delta1    <- -inverse(e1@alpha,e1@delta)*tmp
            epsilon1  <- inverse(transpose(e1@alpha),e1@epsilon)
            delta12   <- inverse(e1@alpha,e1@delta) + epsilon1*forward(transpose(delta1)*e1@delta) -
              delta1*forward(transpose(epsilon1)*e1@delta)
            epsilon12 <- transpose(e1@alpha)*epsilon1 + e1@delta*backward(transpose(e1@epsilon)*epsilon1) + 
              e1@epsilon*forward(transpose(e1@delta)*epsilon1)
            alpha1@g    <- array(0,dim=c(dim(alpha1@g)[1:3],0))
            delta1@g    <- array(0,dim=c(dim(delta1@g)[1:3],0))
            epsilon1@g  <- array(0,dim=c(dim(epsilon1@g)[1:3],0))
            delta12@g   <- array(0,dim=c(dim(delta12@g)[1:3],0))
            epsilon12@g <- array(0,dim=c(dim(epsilon12@g)[1:3],0))
            
            S <- list(alpha1=alpha1,delta1=delta1,epsilon1=epsilon1,
                      alpha2=e1@alpha,delta2=e1@delta,epsilon2=e1@epsilon,
                      delta12=delta12,epsilon12=epsilon12)

            # make inverse
            # S <- inverse(as.operator(e1))
            # S <- symmOperator(alpha=S@alpha,
            #                   delta=S@delta,
            #                   epsilon=S@epsilon)

            # initial energy
          for (iter2 in 1:2) {
            E <- E.symm.id(S)
            dS <- dE.symm.id(S)
            h <- dS
            newS <- S
#            cat("initiation: energy=",E,"\n")
            # decrease energy
            for (iter in 1:1000) {
              # find derivative
              #dS <- dE.symm.id(S,e1)
              # decrease energy along gradient
#              sq.dS <- sum(diag(integral(transpose(dS$alpha)*dS$alpha)+integral(transpose(dS$delta)*dS$delta)+integral(transpose(dS$epsilon)*dS$epsilon)))
#              sq.dS <- sum(diag(integral(transpose(dS$delta)*dS$delta)+integral(transpose(dS$epsilon)*dS$epsilon)))
#              sq.dS <- sum(diag(integral(transpose(h$delta1)*h$delta1)+integral(transpose(h$epsilon1)*h$epsilon1)+
#                                integral(transpose(h$delta12)*h$delta12)+integral(transpose(h$epsilon12)*h$epsilon12)))
#              lambda0 <- 4*2*E/sq.dS
              lambda0 <- 128
              
              for (subs in 1:20) {
                #newS@alpha <- S@alpha-lambda0*dS$alpha
                newS$delta1@f    <- S$delta1@f-lambda0*h$delta1
                newS$epsilon1@f  <- S$epsilon1@f-lambda0*h$epsilon1
                newS$delta12@f   <- S$delta12@f-lambda0*h$delta12
                newS$epsilon12@f <- S$epsilon12@f-lambda0*h$epsilon12
                newE <- E.symm.id(newS)
                if (newE<E) break
                lambda0 <- lambda0/2
              }
              if (iter %% 1 == 0) cat("iteration",iter,": subdivisions=",subs,", lambda=",lambda0,": energy=",newE)
              #cat("\n")
              # update
              if (subs < 20) {
                # conjugate gradient
                dSnew <- dE.symm.id(newS)
                gamma <- (sum((dSnew$delta1-0*dS$delta1)*dSnew$delta1) +
                          sum((dSnew$epsilon1-0*dS$epsilon1)*dSnew$epsilon1) +
                          sum((dSnew$delta12-0*dS$delta12)*dSnew$delta12) +
                          sum((dSnew$epsilon12-0*dS$epsilon12)*dSnew$epsilon12)) /
                  (sum(dS$delta1*dS$delta1)+sum(dS$epsilon1*dS$epsilon1)+
                   sum(dS$delta12*dS$delta12)+sum(dS$epsilon12*dS$epsilon12))
                gamma <- min(1,max(0,gamma))
                cat(", gamma=",gamma,"\n")
                h$delta1    <- dSnew$delta1    + gamma*h$delta1
                h$epsilon1  <- dSnew$epsilon1  + gamma*h$epsilon1
                h$delta12   <- dSnew$delta12   + gamma*h$delta12
                h$epsilon12 <- dSnew$epsilon12 + gamma*h$epsilon12
                dS <- dSnew
                E <- newE
                S <- newS
              } else break
            }        
          }
            
            tmp <- as.operator(symmOperator(alpha=S$alpha1,delta=S$delta1,epsilon=S$epsilon1))*as.operator(e1)
            #plot(tmp@alpha)
            plot(tmp@delta*transpose(tmp@epsilon))
            plot(tmp@delta)
            plot(tmp@epsilon)
            
            dS <- dE.symm.id(S)
            profE <- rep(0,21)
            lambda <- seq(-2*lambda0,2*lambda0,length.out=21)
            for (i in 1:21) {
            #  newS@alpha <- S@alpha+lambda[i]*dS@alpha
              newS$delta1@f    <- S$delta1@f+lambda[i]*dS$delta1
              newS$epsilon1@f  <- S$epsilon1@f+lambda[i]*dS$epsilon1
              newS$delta12@f   <- S$delta12@f+lambda[i]*dS$delta12
              newS$epsilon12@f <- S$epsilon12@f+lambda[i]*dS$epsilon12
              profE[i] <- E.symm.id(newS)
            }
            par(mfrow=c(1,1))
            plot(lambda,profE)
            
          })
