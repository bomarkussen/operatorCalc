approximate <- function(t.obs,matData,index=1:length(t.obs),order=1,continuous=TRUE) {
  # some sanity checking
  if (!is.array(matData) || (length(dim(matData))!=3)) stop("Data entry must be a 3-way array")
  if (length(t.obs)!=dim(matData)[1]) stop("Number of sampling points and samples must be equal")
  if (length(t.obs)<length(index)) stop("There can not be more knot than sampling points")
  if (length(index)<2) stop("There must be a least two knot points")
  if (order<0) stop("The polynomial order must be a non-negative integer")
  
  # piecewise constant function?
  if (order==0) {
    matData <- aperm(matData,c(2,3,1))    
    return(matFct(mesh=mesh(t.obs[index]),
                  f=matData[,,index,drop=FALSE],
                  g=array(aperm(-apply(matData[,,index,drop=FALSE],c(1,2),diff),c(2,3,1)),
                          dim=c(dim(matData)[1:2],length(index)-1,1)),
                  continuous=c(rep(FALSE,length(index)-2),TRUE)))
  }
  # piecewise linear function?
  if (order==1) {
    return(matFct(mesh=mesh(t.obs[index]),
                  f=aperm(matData[index,,,drop=FALSE],c(2,3,1)),
                  g=array(0,dim=c(dim(matData)[2:3],length(index)-1,0)),
                  continuous=rep(TRUE,length(index)-1)))
  }
  # piecewise polynomial function?
  if (order>1) {
    if (length(continuous)!=length(index)-1) continuous <- rep(continuous[1],length(index)-1)
    continuous[length(index)-1] <- TRUE
    g <- .Call('operatorCalc_approximate', PACKAGE = 'operatorCalc', t.obs, matData, index-1, diff(index)-1, continuous, order)
    return(matFct(mesh=mesh(t.obs[index]),
                  f=aperm(matData[index,,,drop=FALSE],c(2,3,1)),
                  g=aperm(array(g,dim=c(dim(matData)[2:3],order,length(index)-1)),
                          c(1,2,4,3)),
                  continuous=continuous))
  }
}
