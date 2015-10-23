\name{mesh-class}
\Rdversion{1.1}
\docType{class}
\alias{mesh}

\title{Class \code{"mesh"}}
\description{
Sampling mesh on an interval [a,b]. A mesh is a strictly increasing sequence of numbers. The squence must have length at least two, where the first and last entries define the left and right endpoints a and b, respectively.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("mesh",x)} or \code{mesh(x)}.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"numeric"}}
  }
}
\section{Extends}{
Class \code{"\linkS4class{numeric}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "numeric", distance 2.
}
\section{Methods}{
No methods defined with class "mesh" in the signature.
}
\author{
Bo Markussen <bomar@math.ku.dk>
}
\examples{
mesh(seq(0,1,0.1))
## An object of class mesh
##  Partition of [0,1] in 10 segments. 5-number statistics of segment lengths:
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##     0.1     0.1     0.1     0.1     0.1     0.1 
}
\keyword{classes}
