\name{approximate}
\alias{approximate}
\title{
Find \code{matFct} object that approximates a data sample
}
\description{
Find \code{matFct} object that approximates a data sample.
}
\usage{
approximate(t.obs, matData, index = 1:length(t.obs), order = 1, continuous = TRUE)
}
\arguments{
  \item{t.obs}{Numeric vector with sampling times. Must have length at least two and be strictly increasing.}
  \item{matData}{A 3-way array with sampling of the function, where the first dimension corresponds to sampling points in \code{t.obs}.}
  \item{index}{Numeric vector with indices of the sampling points that will be used as knot points. Must have length at least two, be strictly increasing, and with values in the range from 1 to \code{length(t.obs)}.}
  \item{order}{Order of the piecewise polynomial interpolation.}
  \item{continuous}{Logical encoding where the resulting function is enforced to be continuous. May also the a vector specifying this separately for each interval. This option only have effect when \code{order} is at least two.}
}
\details{The resulting function is exact at the mesh points (except possibly at the right endpoint of the sampling interval for piecewise constant functions), but in general only approximate between the mesh points. The \code{order} option specifies the order of the polynomial interpolation between mesh points; zero'th order gives piecewise constant and possibly discontinuous functions, first order gives piecewise linear and continuous functions, and second and higher order polynomials are found by least squares. The \code{continuous} option specifies whether the resulting function is enforced to be continuous. This may also be a vector with one entry less than \code{index} specifying the continuity separately for each sampling interval. The \code{continuous} option only has effect when the polynomial order is at least two.}
\value{An object of class \code{matFct}.}
\author{Bo Markussen (bomar@math.ku.dk)}
\note{Functions of class \code{matFct} are always continuous at the right endpoint of the sampling interval. Thus, if \code{order=0} a piecewise constant function is returned, where the value at the right endpoint is taken to be identical to the value at the second last endpoint.}
\seealso{\code{\link{matFct}}}
\examples{
# Make a sample of a random function.
N <- 2
M <- 3
t.obs <- seq(0,1,0.001)
index <- seq(1,length(t.obs),length.out=11)
matData <- apply(array(rnorm(length(t.obs)*N*M)/sqrt(1000),dim=c(length(t.obs),N,M)),2:3,cumsum)
show(approximate(t.obs,matData,index,order=4))

# plot the approximation at different polynomial orders
plot(approximate(t.obs,matData,index,order=0))
plot(approximate(t.obs,matData,index,order=1))
plot(approximate(t.obs,matData,index,order=2))
plot(approximate(t.obs,matData,index,order=8))

}
\keyword{ math }
