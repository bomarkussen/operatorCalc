# methods for standard generic functions: show,+,-,*,plot,[
# methods for added generic functions:
setGeneric("evaluate", function(x, e1, ...) standardGeneric("evaluate"))
setGeneric("superSample", function(x, e1, ...) standardGeneric("superSample"))
setGeneric("forward", function(e1, ...) standardGeneric("forward"))
setGeneric("backward", function(e1, ...) standardGeneric("backward"))
setGeneric("integral", function(e1, ...) standardGeneric("integral"))
setGeneric("triangular", function(e1, ...) standardGeneric("triangular"))
setGeneric("transpose", function(e1, ...) standardGeneric("transpose"))
setGeneric("inverse", function(e1, e2, ...) standardGeneric("inverse"))

# ------------------
# Global options
# ------------------

.onLoad <- function(libname, pkgname) {
  options("minimal.polynomial.order"=2)
  options("maximal.polynomial.order"=8)
}

.onAttach <- function(libname, pkgname) {
  options("minimal.polynomial.order"=2)
  options("maximal.polynomial.order"=8)
}

# ------------------
# Class definition
# ------------------

matFct <- setClass("matFct", slots=c(mesh="mesh",f="array",g="array",continuous="logical"))

setValidity("matFct", function(object) {
  msg <- NULL
  valid <- TRUE
  if (length(dim(object@f))!=3) {
    valid <- FALSE
    msg <- c(msg,"Slot f must be a 3-way array")
  }
  if (length(dim(object@g))!=4) {
    valid <- FALSE
    msg <- c(msg,"Slot g must be a 4-way array")
  }
  if (dim(object@f)[3]!=length(object@mesh)) {
    valid <- FALSE
    msg <- c(msg,"Number of sampling points in f must equal mesh size")
  }
  if (dim(object@g)[3]!=length(object@mesh)-1) {
    valid <- FALSE
    msg <- c(msg,"Number of sampling points in g must equal mesh size - 1")
  }
  if (any(dim(object@f)[1:2]!=dim(object@g)[1:2])) {
    valid <- FALSE
    msg <- c(msg,"First two dimensions of f and g must be equal")
  }
  if (length(object@continuous)!=dim(object@g)[3]) {
    valid <- FALSE
    msg <- c(msg,"Length of continuity attribute must equal mesh size - 1")
  }
  
  if (valid) TRUE else msg
})

# ---------
# Methods
# ---------

setMethod("show",
          signature="matFct",
          definition=function(object) {
            cat("An object of class ",class(object),"\n",sep="")
            cat(" Piecewise polynominal (order=",max(1,dim(object@g)[4]),") function: [",
                object@mesh[1],",",object@mesh[length(object@mesh)],"] -> ",
                paste(dim(object@f)[1:2],collapse="x")," matrices\n",sep="")
            if (all(object@continuous)) {
              cat(" Supposed to be continuous in approximate computations\n") 
            } else {
              cat(" Supposed to be discontinuous at ",sum(!object@continuous)," mesh points in approximate computations\n",sep="")
            }
            invisible(NULL)
          })

setMethod("[",
          signature=c(x="matFct",i="numeric",j="numeric"),
          definition=function(x,i,j) {
            matFct(mesh=x@mesh,f=x@f[i,j,,drop=FALSE],g=x@g[i,j,,,drop=FALSE],continuous=x@continuous)
          })

setMethod("plot",
          signature=c(x="numeric",y="matFct"),
          definition=function(x,y) {
            # evaluate function
            x <- sort(unique(c(pmax(min(x),pmin(max(x),y@mesh)),x)))
            yy <- evaluate(x,y)
            # make plot
            par(mfrow=dim(y@f)[1:2])
            for (n in 1:dim(y@f)[1]) for (m in 1:dim(y@f)[2]) {
              plot(x,yy[n,m,],type="n",xlab="time",ylab=paste(c("f[",n,",",m,"]"),collapse=""))
              ii <- 1
              for (jj in grep(FALSE,c(y@continuous[-length(y@continuous)],FALSE))) {
                kk <- match(y@mesh[jj+1],x)
                if (is.na(kk)) kk <- length(x)
                points(x[ii],yy[n,m,ii],pch=16)
                lines(x[ii:(kk-1)],yy[n,m,ii:(kk-1)])
                ii <- kk
                if (kk==length(x)) break
              }
            }
          })

setMethod("plot",
          signature=c(x="matFct",y="missing"),
          definition=function(x) {
            plot(seq(x@mesh[1],x@mesh[length(x@mesh)],length.out=200),x)            
          })

setMethod("evaluate",
          signature=c(x="numeric",e1="matFct"),
          definition=function(x,e1) {
            if (prod(dim(e1@f)[1:2])==0) {
              # trivial result
              return(array(0,dim=c(dim(e1@f)[1:2],length(x))))
            } else {
              # non-trivial result
              x <- pmin(pmax(x,e1@mesh[1]),e1@mesh[length(e1@mesh)])
              j <- findInterval(x,e1@mesh,rightmost.closed=TRUE)
              s <- (x-e1@mesh[j])/(e1@mesh[j+1]-e1@mesh[j])
              # call Rcpp
              return(.Call('operatorCalc_evaluate', PACKAGE = 'operatorCalc', j-1, s, e1@f, e1@g))
            }
          })

setMethod("superSample",
          signature=c(x="numeric",e1="matFct"),
          definition=function(x,e1) {
            # is the function extended to the left?
            left <- NULL
            if (min(x) < min(e1@mesh)) {
              left <- min(x)
              x <- x[x >= min(e1@mesh)]
            }
            # is the function extended to the right?
            right <- NULL
            if (max(x) > max(e1@mesh)) {
              right <- max(x)
              x <- x[x <= max(e1@mesh)]
            }
            # take dimensions
            N <- dim(e1@f)[1]
            M <- dim(e1@f)[2]
            r <- dim(e1@g)[4]
            x <- sort(union(x,e1@mesh))
            p <- length(x)
            if ((N==0) || (M==0)) {
              # trivial result
              x <- c(left,x,right)
              p <- length(x)
              return(matFct(mesh=mesh(x),f=array(0,dim=c(N,M,p)),
                            g=array(0,dim=c(N,M,p-1,0)),continuous=rep(TRUE,p-1)))
            } else {
              # non-trivial result
              # construct the super sampled function inside the old range
              j <- findInterval(x[-p],e1@mesh,rightmost.closed=TRUE)
              a <- (x[-p]-e1@mesh[j])/(e1@mesh[j+1]-e1@mesh[j])
              b <- (x[-1]-x[-p])/(e1@mesh[j+1]-e1@mesh[j])
              res <- .Call('operatorCalc_superSample', PACKAGE = 'operatorCalc', outer(1:r,1:r,choose), j-1, a, b, e1@f, e1@g)
              f <- res$f
              g <- array(res$g,dim=c(N,M,p-1,r))
              # construct continuity property
              continuous <- rep(TRUE,p-1)
              continuous[is.element(x[-1],e1@mesh)] <- e1@continuous
              # extend the function to the left if necessary
              if (!is.null(left)) {
                x <- c(left,x)
                f <- array(c(f[,,1],f),dim=c(N,M,p+1))
                g <- aperm(array(c(rep(0,N*M*r),aperm(g,c(1,2,4,3))),dim=c(N,M,r,p)),c(1,2,4,3))
                continuous <- c(TRUE,continuous)
                p <- p+1
              }
              # extend the function to the right if necessary
              if (!is.null(right)) {
                x <- c(x,right)
                f <- array(c(f,f[,,p]),dim=c(N,M,p+1))
                g <- aperm(array(c(aperm(g,c(1,2,4,3)),rep(0,N*M*r)),dim=c(N,M,r,p)),c(1,2,4,3))
                continuous <- c(continuous,TRUE)
                p <- p+1
              }
              # return result
              return(matFct(mesh=mesh(x),f=f,g=g,continuous=continuous))
            }
          })

setMethod("-",
          signature=c(e1="matFct"),
          definition=function(e1) {
            matFct(mesh=e1@mesh,f=-e1@f,g=-e1@g,continuous=e1@continuous)
          })

setMethod("+",
          signature=c(e1="matFct",e2="matFct"),
          definition=function(e1,e2) {
            # sanity check
            if (any(dim(e1@f)[1:2]!=dim(e2@f)[1:2])) 
              stop(paste("Mismatch of dimensions (",dim(e1@f)[1],",",dim(e1@f)[2],") and (",dim(e2@f)[1],",",dim(e2@f)[2],") in matrix sum",sep=""))
            # if e1 and e2 have different meshes then super sample both functions
            if (!setequal(e1@mesh,e2@mesh)) {
              e1 <- superSample(e2@mesh,e1)
              e2 <- superSample(e1@mesh,e2)
            }
            # find g function
            if (dim(e1@g)[4]==dim(e2@g)[4]) g <- e1@g+e2@g else {
              if (dim(e1@g)[4]<dim(e2@g)[4]) {
                q <- dim(e1@g)[4]
                g <- e2@g
                if (q>0) g[,,,1:q] <- g[,,,1:q,drop=FALSE]+e1@g
              } else {
                q <- dim(e2@g)[4]
                g <- e1@g
                if (q>0) g[,,,1:q] <- g[,,,1:q,drop=FALSE]+e2@g
              }
            }
            # return result
            matFct(mesh=e1@mesh,f=e1@f+e2@f,g=g,continuous=e1@continuous&e2@continuous)
          })

setMethod("-",
          signature=c(e1="matFct",e2="matFct"),
          definition=function(e1,e2) {
            # sanity check
            if (any(dim(e1@f)[1:2]!=dim(e2@f)[1:2])) 
              stop(paste("Mismatch of dimensions (",dim(e1@f)[1],",",dim(e1@f)[2],") and (",dim(e2@f)[1],",",dim(e2@f)[2],") in matrix differece",sep=""))
            # if e1 and e2 have different meshes then super sample both functions
            if (!setequal(e1@mesh,e2@mesh)) {
              e1 <- superSample(e2@mesh,e1)
              e2 <- superSample(e1@mesh,e2)
            }
            # find g function
            if (dim(e1@g)[4]==dim(e2@g)[4]) g <- e1@g-e2@g else {
              if (dim(e1@g)[4]<dim(e2@g)[4]) {
                q <- dim(e1@g)[4]
                g <- -e2@g
                if (q>0) g[,,,1:q] <- g[,,,1:q,drop=FALSE]+e1@g
              } else {
                q <- dim(e2@g)[4]
                g <- e1@g
                if (q>0) g[,,,1:q] <- g[,,,1:q,drop=FALSE]-e2@g
              }
            }
            # return result
            matFct(mesh=e1@mesh,f=e1@f-e2@f,g=g,continuous=e1@continuous&e2@continuous)
          })

setMethod("|",
          signature=c(e1="matFct",e2="matFct"),
          definition=function(e1,e2) {
            # sanity check
            if (dim(e1@f)[1]!=dim(e2@f)[1])
              stop(paste("Mismatch of dimensions (",dim(e1@f)[1],",",dim(e1@f)[2],") and (",dim(e2@f)[1],",",dim(e2@f)[2],") in column concatenation",sep=""))
            # if e1 and e2 have different meshes then super sample both functions
            if (!setequal(e1@mesh,e2@mesh)) {
              e1 <- superSample(e2@mesh,e1)
              e2 <- superSample(e1@mesh,e2)
            }
            # find f function
            N  <- dim(e1@f)[1]
            M1 <- dim(e1@f)[2]
            M2 <- dim(e2@f)[2]
            p  <- dim(e1@f)[3]
            f  <- aperm(array(c(aperm(e1@f,c(1,3,2)),
                                aperm(e2@f,c(1,3,2))),
                              dim=c(N,p,M1+M2)),c(1,3,2))
            # find g function
            r1 <- dim(e1@g)[4]
            r2 <- dim(e2@g)[4]
            if (r1==r2) {
              # e1 and e2 have same r
              g <- aperm(array(c(aperm(e1@g,c(1,3,4,2)),
                                 aperm(e2@g,c(1,3,4,2))),
                               dim=c(N,p-1,r1,M1+M2)),c(1,4,2,3))
            } else {
              if (r1>r2) {
                # zero pad forth dimension of e2@g to match that of e1@g
                g <- aperm(array(c(aperm(e1@g,c(1,3,4,2)),
                                   aperm(array(c(e2@g,rep(0,N*M2*(p-1)*(r1-r2))),dim=c(N,M2,p-1,r1)),c(1,3,4,2))),
                                   dim=c(N,p-1,r1,M1+M2)),c(1,4,2,3))
              } else {
                # zero pad forth dimension of e1@g to match that of e2@g
                g <- aperm(array(c(aperm(array(c(e1@g,rep(0,N*M1*(p-1)*(r2-r1))),dim=c(N,M1,p-1,r2)),c(1,3,4,2)),
                                   aperm(e2@g,c(1,3,4,2))),
                                 dim=c(N,p-1,r2,M1+M2)),c(1,4,2,3))
              }
            }
            # return result
            matFct(mesh=e1@mesh,f=f,g=g,continuous=e1@continuous&e2@continuous)
          })

setMethod("*",
          signature=c(e1="matFct",e2="matFct"),
          definition=function(e1,e2) {
            # sanity check
            if (dim(e1@f)[2]!=dim(e2@f)[1]) 
              stop(paste("Mismatch of dimensions (",dim(e1@f)[1],",",dim(e1@f)[2],") and (",dim(e2@f)[1],",",dim(e2@f)[2],") in matrix product",sep=""))
            # if e1 and e2 have different meshes then super sample both functions
            if (!setequal(e1@mesh,e2@mesh)) {
              e1 <- superSample(e2@mesh,e1)
              e2 <- superSample(e1@mesh,e2)
            }
            # take dimensions
            N <- dim(e1@f)[1]
            M <- dim(e2@f)[2]
            p <- length(e1@mesh)
            if ((dim(e1@f)[2]==0) || (N==0) || (M==0)) {
              # trivial answer due to empty inner sum or empty marginals
              return(matFct(mesh=e1@mesh,f=array(0,dim=c(N,M,p)),g=array(0,dim=c(N,M,p-1,0)),
                            continuous=rep(TRUE,p-1)))

            } else {
              # non-trivial answer
              # Call C function to find (f,g)
              res <- .Call('operatorCalc_multFct', PACKAGE = 'operatorCalc', e1@f, e1@g, e2@f, e2@g)
              attr(res$g,"dim") <- c(N,M,p-1,max(1,dim(e1@g)[4])+max(1,dim(e2@g)[4]))
              return(matFct(mesh=e1@mesh,f=res$f,g=res$g,continuous=e1@continuous&e2@continuous))
            }
          })

setMethod("*",
          signature=c(e1="matrix",e2="matFct"),
          definition=function(e1,e2) {
            # sanity check
            if (dim(e1)[2]!=dim(e2@f)[1])
              stop(paste("Mismatch of dimensions (",dim(e1)[1],",",dim(e1)[2],") and (",dim(e2@f)[1],",",dim(e2@f)[2],") in matrix sum",sep=""))
            # take dimensions
            N <- dim(e1)[1]
            M <- dim(e2@f)[2]
            p <- length(e2@mesh)
            if ((dim(e1)[2]==0) || (N==0) || (M==0)) {
              # trivial answer due to empty inner sum or empty marginals
              return(matFct(mesh=e2@mesh,f=array(0,dim=c(N,M,p)),g=array(0,dim=c(N,M,p-1,0)),
                            continuous=rep(TRUE,p-1)))
              
            } else {
              # non-trivial answer
              # Call C function to find (f,g) and return result
              res <- .Call('operatorCalc_multFct_left', PACKAGE = 'operatorCalc', e1, e2@f, e2@g)
              attr(res$g,"dim") <- c(N,M,p-1,dim(e2@g)[4])
              return(matFct(mesh=e2@mesh,f=res$f,g=res$g,continuous=e2@continuous))
            }
          })

setMethod("*",
          signature=c(e1="matFct",e2="matrix"),
          definition=function(e1,e2) {
            # sanity check
            if (dim(e1@f)[2]!=dim(e2)[1])
              stop(paste("Mismatch of dimensions (",dim(e1@f)[1],",",dim(e1@f)[2],") and (",dim(e2)[1],",",dim(e2)[2],") in matrix sum",sep=""))
            # take dimensions
            N <- dim(e1@f)[1]
            M <- dim(e2)[2]
            p <- length(e1@mesh)
            if ((dim(e2)[1]==0) || (N==0) || (M==0)) {
              # trivial answer due to empty inner sum or empty marginals
              return(matFct(mesh=e1@mesh,f=array(0,dim=c(N,M,p)),g=array(0,dim=c(N,M,p-1,0)),
                            continuous=rep(TRUE,p-1)))
            } else {
              # non-trivial answer
              # Call C function to find (f,g) and return result
              res <- .Call('operatorCalc_multFct_right', PACKAGE = 'operatorCalc', e1@f, e1@g, e2)
              attr(res$g,"dim") <- c(N,M,p-1,dim(e1@g)[4])
              return(matFct(mesh=e1@mesh,f=res$f,g=res$g,continuous=e1@continuous))
            }
          })

setMethod("*",
          signature=c(e1="numeric",e2="matFct"),
          definition=function(e1,e2) {
            matFct(mesh=e2@mesh,f=e1*(e2@f),g=e1*(e2@g),continuous=e2@continuous)
          })

setMethod("*",
          signature=c(e1="matFct",e2="numeric"),
          definition=function(e1,e2) {
            matFct(mesh=e1@mesh,f=(e1@f)*e2,g=(e1@g)*e2,continuous=e1@continuous)
          })

setMethod("forward",
          signature=c(e1="matFct"),
          definition=function(e1) {
            if (prod(dim(e1@f)[1:2])==0) {
              # trivial result
              return(matFct(mesh=e1@mesh,
                            f=array(0,dim=dim(e1@f)),
                            g=array(0,dim=c(dim(e1@g)[1:3],0)),
                            continuous=rep(TRUE,dim(e1@g)[3])))
            } else {
              # non-trivial result
              # call Rcpp and return result
              res <- .Call('operatorCalc_forward', PACKAGE = 'operatorCalc', diff(unclass(e1@mesh)), e1@f, e1@g)
              attr(res$g,"dim") <- c(dim(e1@g)[1:3],max(1,dim(e1@g)[4])+1)
              return(matFct(mesh=e1@mesh,f=res$f,g=res$g,continuous=rep(TRUE,dim(e1@g)[3])))
            }
          })

setMethod("backward",
          signature=c(e1="matFct"),
          definition=function(e1) {
            if (prod(dim(e1@f)[1:2])==0) {
              # trivial result
              return(matFct(mesh=e1@mesh,
                            f=array(0,dim=dim(e1@f)),
                            g=array(0,dim=c(dim(e1@g)[1:3],0)),
                            continuous=rep(TRUE,dim(e1@g)[3])))
            } else {
              # non-trivial result
              # call Rcpp and return result
              res <- .Call('operatorCalc_backward', PACKAGE = 'operatorCalc', diff(unclass(e1@mesh)), e1@f, e1@g)
              attr(res$g,"dim") <- c(dim(e1@g)[1:3],max(1,dim(e1@g)[4])+1)
              return(matFct(mesh=e1@mesh,f=res$f,g=res$g,continuous=rep(TRUE,dim(e1@g)[3])))
            }
          })

setMethod("integral",
          signature=c(e1="matFct"),
          definition=function(e1) {
            # call Rcpp and return result
            .Call('operatorCalc_integral', PACKAGE = 'operatorCalc', diff(unclass(e1@mesh)), e1@f, e1@g)
          })

setMethod("triangular",
          signature=c(e1="matFct"),
          definition=function(e1) {
            if (prod(dim(e1@f)[1:2])==0) {
              # trivial result
              return(f=array(0,dim=dim(e1@f)[1:3]))
            } else {
              # non-trivial result
              # call Rcpp and return result
              return(.Call('operatorCalc_triangular', PACKAGE = 'operatorCalc', diff(unclass(e1@mesh)), e1@f, e1@g))
            }
          })

setMethod("transpose",
          signature=c(e1="matFct"),
          definition=function(e1) {
            matFct(mesh=e1@mesh,f=aperm(e1@f,c(2,1,3)),g=aperm(e1@g,c(2,1,3,4)),continuous=e1@continuous)
          })

setMethod("inverse",
          signature=c(e1="matFct",e2="missing"),
          definition=function(e1) {
            # sanity check
            if (dim(e1@f)[1]!=dim(e1@f)[2])
              stop(paste("First argument of dimension (",dim(e1@f)[1],",",dim(e1@f)[2],") is non-square",sep=""))
            # Take properties from matFct object and global options
            order <- max(options()$"minimal.polynomial.order",min(options()$"maximal.polynomial.order",dim(e1@g)[4]))
            ifelse(order>1, super <- 2*order+1, super <- 0)
            continuous <- e1@continuous
            # Call Cpp
            res <- .Call('operatorCalc_inverse', PACKAGE = 'operatorCalc', e1@f, e1@g, order, super, continuous)
            # fix Rcpp result: piecewise constant function?
            if (order==0) {
              res$g <- array(res$g[,1:dim(e1@g)[3]],dim=c(dim(e1@g)[1:2],1,dim(e1@g)[3]))
              continuous <- c(rep(FALSE,dim(e1@g)[3]-1),TRUE)
            }
            # fix Rcpp result: piecewise linear function?
            if (order==1) {
              res$g <- array(0,dim=c(dim(e1@g)[1:2],0,dim(e1@g)[3]))
              continuous <- rep(TRUE,dim(e1@g)[3])
            }
            # fix Rcpp result: piecewise polynomial function?
            if (order>1) attr(res$g,"dim") <- c(dim(e1@g)[1:2],order,dim(e1@g)[3])
            # return result
            return(matFct(mesh=e1@mesh,f=res$f,g=aperm(res$g,c(1,2,4,3)),continuous=continuous))
          })

setMethod("inverse",
          signature=c(e1="matFct", e2="matFct"),
          definition=function(e1, e2) {
            # sanity check
            if (dim(e1@f)[1]!=dim(e1@f)[2])
              stop(paste("First argument of dimension (",dim(e1@f)[1],",",dim(e1@f)[2],") is non-square",sep=""))
            if (dim(e1@f)[2]!=dim(e2@f)[1])
              stop(paste("Mismatch of dimensions (",dim(e1@f)[1],",",dim(e1@f)[2],") and (",dim(e2@f)[1],",",dim(e2@f)[2],") in linear equation",sep=""))
            # Take properties from matFct object and global options
            order <- max(options()$"minimal.polynomial.order",min(options()$"maximal.polynomial.order",dim(e1@g)[4]+dim(e2@g)[4]))
            ifelse(order>1, super <- 2*order+1, super <- 0)
            # if e1 and e2 have different meshes then super sample both functions
            if (!setequal(e1@mesh,e2@mesh)) {
              e1 <- superSample(e2@mesh,e1)
              e2 <- superSample(e1@mesh,e2)
            }
            # pointwise solution of linear equation
            res <- .Call('operatorCalc_solveFct', PACKAGE = 'operatorCalc', e1@f, e1@g, e2@f, e2@g, order, super, e1@continuous&e2@continuous)
            attr(res$g,"dim") <- c(dim(e1@g)[1],dim(e2@g)[2],order,dim(e1@g)[3])
            # return result
            return(matFct(mesh=e1@mesh,f=res$f,g=aperm(res$g,c(1,2,4,3)),continuous=e1@continuous&e2@continuous))
          })
