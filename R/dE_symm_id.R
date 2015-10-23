dE.symm.id <- function(S) {
  # S: list: alpha1, delta1, epsilon1, alpha2, delta2, epsilon2, delta12, epsilon12

  # compute representation of S1*S2
  alpha12 <- S$alpha1*S$alpha2; alpha12@f <- alpha12@f - c(diag(dim(alpha12@f)[1]))
  delta12 <- inverse(S$alpha2,S$delta2) + S$epsilon1*forward(transpose(S$delta1)*S$delta2) -
    S$delta1*forward(transpose(S$epsilon1)*S$delta2)
  epsilon12 <- transpose(S$alpha2)*S$epsilon1 + S$delta2*backward(transpose(S$epsilon2)*S$epsilon1) + 
    S$epsilon2*forward(transpose(S$delta2)*S$epsilon1)

  # derivative wrt alpha1
  Dalpha1 <- alpha12*transpose(S$alpha2) +
    (delta12-S$delta12)*transpose(S$delta2)
  
  # derivative wrt delta1
  Ddelta1 <- S$delta1*backward(transpose(S$epsilon12)*S$epsilon12) + 
    S$delta12*backward(transpose(S$epsilon2)*S$epsilon12) +
    S$delta2*backward(transpose(delta12-S$delta12)*S$epsilon1) -
    (delta12-S$delta12)*forward(transpose(S$delta2)*S$epsilon1)
  
  # derivative wrt epsilon1
  Depsilon1 <- (delta12-S$delta12)*forward(transpose(S$delta2)*S$delta1) -
    S$delta2*backward(transpose(delta12-S$delta12)*S$delta1) +
    S$alpha2*(epsilon12-S$epsilon12) +
    S$epsilon2*forward(transpose(S$delta2)*(epsilon12-S$epsilon12)) +
    S$delta2*backward(transpose(S$epsilon2)*(epsilon12-S$epsilon12))
  
  # derivative wrt alpha2
  Dalpha2 <- transpose(S$alpha1)*alpha12 +
    S$epsilon1*(epsilon12-S$epsilon12)
  
  # derivative wrt delta2
  Ddelta2 <- transpose(S$alpha1)*(delta12-S$delta12) +
    S$delta1*backward(transpose(S$epsilon1)*(delta12-S$delta12)) -
    S$epsilon1*backward(transpose(S$delta1)*(delta12-S$delta12)) +
    (epsilon12-S$epsilon12)*backward(transpose(S$epsilon1)*S$epsilon2) +
    S$epsilon1*backward(transpose(epsilon12-S$epsilon12)*S$epsilon2)
  
  # derivative wrt epsilon2
  Depsilon2 <- S$epsilon12*forward(transpose(S$delta1)*S$delta12) +
    S$epsilon2*forward(transpose(S$delta12)*S$delta12) + 
    S$epsilon1*forward(transpose(epsilon12-S$epsilon12)*S$delta2) +
    S$delta1*backward(transpose(S$epsilon1)*(epsilon12-S$epsilon12))
  
  # derivative wrt delta12
  Ddelta12 <- S$delta1*backward(transpose(S$epsilon12)*S$epsilon2) +
    S$delta12*backward(transpose(S$epsilon2)*S$epsilon2) +
    (S$delta12 - delta12)
  
  # derivative wrt epsilon12
  Depsilon12 <- S$epsilon12*forward(transpose(S$delta1)*S$delta1) +
    S$epsilon2*forward(transpose(S$delta12)*S$delta1) +
    (S$epsilon12 - epsilon12)
 
  # make results piecewise linear
  #Dalpha1@g    <- array(0,dim=c(dim(Dalpha1@g)[1:3],0))
  #Ddelta1@g    <- array(0,dim=c(dim(Ddelta1@g)[1:3],0))
  #Depsilon1@g  <- array(0,dim=c(dim(Depsilon1@g)[1:3],0))
  #Dalpha2@g    <- array(0,dim=c(dim(Dalpha2@g)[1:3],0))
  #Ddelta2@g    <- array(0,dim=c(dim(Ddelta2@g)[1:3],0))
  #Depsilon2@g  <- array(0,dim=c(dim(Depsilon2@g)[1:3],0))
  #Ddelta12@g   <- array(0,dim=c(dim(Ddelta12@g)[1:3],0))
  #Depsilon12@g <- array(0,dim=c(dim(Depsilon12@g)[1:3],0))
  
  # return result
  return(list(alpha1=triangular(Dalpha1),delta1=triangular(Ddelta1),epsilon1=triangular(Depsilon1),
              alpha2=triangular(Dalpha2),delta2=triangular(Ddelta2),epsilon2=triangular(Depsilon2),
              delta12=triangular(Ddelta12),epsilon12=triangular(Depsilon12)))
}
