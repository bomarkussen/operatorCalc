\name{matFct-class}
\Rdversion{1.1}
\docType{class}
\alias{matFct}
\alias{+,matFct,matFct-method}
\alias{-,matFct,matFct-method}
\alias{-,matFct,ANY-method}
\alias{*,matFct,matFct-method}
\alias{*,matrix,matFct-method}
\alias{*,matFct,matrix-method}
\alias{*,numeric,matFct-method}
\alias{*,matFct,numeric-method}
\alias{|,matFct,matFct-method}
\alias{"[",matFct,numeric,numeric-method}
\alias{transpose,matFct-method}
\alias{inverse,matFct,numeric,logical-method}
\alias{solve,matFct,matFct-method}
\alias{forward,matFct-method}
\alias{backward,matFct-method}
\alias{integral,matFct-method}
\alias{triangular,matFct-method}
\alias{evaluate,numeric,matFct-method}
\alias{superSample,numeric,matFct-method}
\alias{plot,matFct,missing-method}
\alias{plot,numeric,matFct-method}

\title{Class \code{"matFct"}}
\description{
Piecewise polynomial continuous matrix value function defined on an interval.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("matFct", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{mesh}:}{Object of class \code{"mesh"}}
    \item{\code{f}:}{Object of class \code{"array"}}
    \item{\code{g}:}{Object of class \code{"array"}}
    \item{\code{continuous}:}{Object of class \code{"logical"}}
  }
}
\section{Methods}{
  \describe{
    \item{+}{\code{signature(e1 = "matFct", e2 = "matFct")}: Pointwise sum of two matrix valued functions.}
    \item{-}{\code{signature(e1 = "matFct", e2 = "matFct")}: Pointwise difference of two matrix valued functions.}    
    \item{*}{\code{signature(x = "matFct", y = "matFct")}: Pointwise matrix multiplication of two matrix valued functions.}
    \item{*}{\code{signature(x = "matrix", y = "matFct")}: Pointwise matrix multiplication of a constant matrix and a matrix valued function.}
    \item{*}{\code{signature(x = "matFct", y = "matrix")}: Pointwise matrix multiplication of a matrix valued function and a constant matrix.}
    \item{*}{\code{signature(x = "numeric", y = "matFct")}: Pointwise multiplication of a constant scalar and a matrix valued function.}
    \item{*}{\code{signature(x = "matFct", y = "numeric")}: Pointwise multiplication of a matrix valued function and a constant scalar.}
    \item{inverse}{\code{signature(e1 = "matFct", order = "numeric", continuous = "logical")}: Pointwise inverse of a matrix valued function. The inverse is exact at the mesh points, but in general only approximate between mesh points. The \code{order} argument controls the order of the polynomial interpolation between mesh points, and the \code{continuous} argument controls whether the result is enforced to be continuous.}
    \item{inverse}{\code{signature(a = "matFct", b = "matFct")}: Pointwise solution of linear equation. The solution is exact at the mesh points, but in general only approximate between mesh points.}
    \item{|}{\code{signature(e1 = "matFct", e2 = "matFct")}: Pointwise concatenation of the columns of two matrix valued functions.}
    \item{"["}{\code{signature(e1 = "matFct", i = "numeric", j = "numeric")}: Subsetting of matrix value.}
    \item{transpose}{\code{signature(e1 = "matFct")}: Pointwise transpose of a matrix valued function.}
    \item{forward}{\code{signature(e1 = "matFct")}: Forward integral of a matrix valued function, i.e. the definite integral from the left end point of the domain to the argument.}
    \item{backward}{\code{signature(e1 = "matFct")}: Backward integral of a matrix valued function, i.e. the definite integral from the argument to the right end point of the domain.}
    \item{integral}{\code{signature(e1 = "matFct")}: Integral of a matrix valued function, i.e. the definite integral from the left to the right end point of the domain.}
    \item{triangular}{\code{signature(e1 = "matFct")}: Returns vector of so-called triangular integrals. To be elaborated later.}    
    \item{evaluate}{\code{signature(x = "numeric", e1 = "matFct")}: Evaluation of \code{matFct} at the values \code{x}. Returns an \code{(N,M,p)-array}, where (N,M) is the dimension of the value space and p is the length of \code{x}.}
    \item{superSample}{\code{signature(x = "numeric", e1 = "matFct")}: Resamples \code{matFct} such that the new mesh equals the union of \code{x} and the old mesh. If the new mesh extends the function domain, then the function is extended to be constant in the extended domain.}
    \item{plot}{\code{signature(x = "matFct")}: Plot a matrix valued function on its domain.}
    \item{plot}{\code{signature(x = "numeric", y = "matFct")}: Plot a matrix valued function at the points specified in \code{x}.}
	 }
}
\note{
If the matrix dimensions of the arguments do not match the mathematical operation, e.g. the second and the first dimensions of the factors in a matrix product are not equal, then the operation is restricted to the left most columns or upper most rows, respectively, and an associated warning is issued.

Formally the '\code{|}' method belongs to the S4 group \code{Logic}. In our context this implies that the concatenation is done after matrix multiplication and sum, unless specified otherwise using parentheses.
}
%%\references{
%%  ~~put references to the literature/web site here~~
%%}
\author{
Bo Markussen <bomar@math.ku.dk>
}
\seealso{
\code{\link{mesh}}
%%\code{\linkS4class{mesh-class}}
}
\examples{
showClass("matFct")
}
\keyword{classes}
