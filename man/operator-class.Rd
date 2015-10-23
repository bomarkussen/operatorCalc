\name{operator-class}
\Rdversion{1.1}
\docType{class}
\alias{operator}
\alias{-,operator,ANY-method}
\alias{-,operator,operator-method}
\alias{*,operator,matFct-method}
\alias{*,operator,operator-method}
\alias{+,operator,operator-method}
\alias{"[",operator,numeric-method}
\alias{inverse,operator-method}
\alias{superSample,numeric,operator-method}

\title{Class \code{"operator"}}
\description{
So-called lattice operators, which are an algebra of integral operators on functions on a compact interval and with possibly multivariate values.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("operator", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{alpha}:}{Object of class \code{"matFct"}}
    \item{\code{beta}:}{Object of class \code{"matFct"}}
    \item{\code{gamma}:}{Object of class \code{"matFct"}}
    \item{\code{delta}:}{Object of class \code{"matFct"}}
    \item{\code{epsilon}:}{Object of class \code{"matFct"}}
  }
}
\section{Methods}{
  \describe{
    \item{+}{\code{signature(e1 = "operator", e2 = "operator")}: Sum of two lattice operators.}
    \item{-}{\code{signature(e1 = "operator", e2 = "operator")}: Difference of two lattice operators.}
    \item{-}{\code{signature(e1 = "operator")}: Change sign of a lattice operator.}
    \item{*}{\code{signature(e1 = "operator", e2 = "operator")}: Product of two lattice operators.}
    \item{*}{\code{signature(e1 = "operator", e2 = "matFct")}: Application of a lattice operator on a function.}
    \item{inverse}{\code{signature(e1 = "operator")}: Inverse of lattice operator.}
    \item{"["}{\code{signature(e1 = "operator", i = "numeric")}: Subsetting of value space.}
    \item{superSample}{\code{signature(x = "numeric", e1 = "operator")}: Resamples \code{operator} such that the new mesh equals the union of \code{x} and the old mesh. If the new mesh extends the function domain, then the operator coefficients are extended to be constant in the extended domain.}    
	 }
}
\references{See vignette.}
\author{Bo Markussen <bomar@math.ku.dk>}
\note{Dimensions must match when lattice operators are used.}
\seealso{
See Also as \code{\link{matFct}},
\code{\linkS4class{matFct}}.
}
\examples{
# -----------------------------------------------------
# Functions to generate random functions and operators
# -----------------------------------------------------

random.function <- function(a,b,p,r,N,M) {
  pp <- 2*(r+1)*p
  approximate(seq(a,b,length.out=pp+1),
              array(apply(array(rnorm((pp+1)*N*M)/sqrt(pp),dim=c(pp+1,N,M)),2:3,cumsum),dim=c(pp+1,N,M)),
              seq(1,pp+1,2*(r+1)),
              r,TRUE)}
              
random.operator <- function(a,b,p,r,N,T.order,P.order) {
  tmp <- random.function(a,b,p,r,N,N)
  tmp <- tmp*transpose(tmp)
  tmp@f <- tmp@f + 10*c(diag(N))
  operator(alpha=tmp,
           beta=random.function(a,b,p,r,N,T.order),
           gamma=random.function(a,b,p,r,N,T.order),
           delta=random.function(a,b,p,r,N,P.order),
           epsilon=random.function(a,b,p,r,N,P.order))}
           
# ------------------------
# Check of the algebra
# ------------------------

a <- -1
b <- 2
p <- 20
r <- 2
N <- 2

fct <- random.function(a,b,p,r,N,3)
op1 <- random.operator(a,b,p,r,N,2,1)
op2 <- random.operator(a,b,p,r,N,2,3)

plot(op1*op2*fct)

# The associative rule is satisfied up to numerical precision
x <- seq(0,1,0.01)
summary(evaluate(x,(op1*op2)*fct)-evaluate(x,op1*(op2*fct)))
summary(evaluate(x,(op1+op2)*fct)-evaluate(x,op1*fct+op2*fct))

# Operator inversion: Not yet sufficiently stable
system.time(tmp <- op1*inverse(op1))
plot(tmp@alpha)
plot(tmp@beta*transpose(tmp@gamma))
plot(tmp@delta*transpose(tmp@epsilon))
summary(evaluate(x,tmp*fct)-evaluate(x,fct))
}
\keyword{classes}
\keyword{math}
