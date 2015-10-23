E.symm.id <- function(S) {
  # S: list: alpha1, delta1, epsilon1, alpha2, delta2, epsilon2, delta12, epsilon12
  
  # compute representation of S1*S2
  alpha12 <- S$alpha1*S$alpha2; alpha12@f <- alpha12@f - c(diag(dim(alpha12@f)[1]))
  delta12 <- S$alpha1*S$delta2 + 
    S$epsilon1*forward(transpose(S$delta1)*S$delta2) -
    S$delta1*forward(transpose(S$epsilon1)*S$delta2)
  epsilon12 <- transpose(S$alpha2)*S$epsilon1 + 
    S$delta2*backward(transpose(S$epsilon2)*S$epsilon1) + 
    S$epsilon2*forward(transpose(S$delta2)*S$epsilon1)

  # return energy
  return(0.5*sum(diag(integral(transpose(alpha12)*alpha12)))+
         0.5*sum(diag(integral(S$epsilon12*forward(transpose(S$delta1)*S$delta1)*transpose(S$epsilon12))))+
         0.5*sum(diag(integral(S$epsilon12*forward(transpose(S$delta1)*S$delta12)*transpose(S$epsilon2))))+
         0.5*sum(diag(integral(S$epsilon2*forward(transpose(S$delta12)*S$delta1)*transpose(S$epsilon12))))+
         0.5*sum(diag(integral(S$epsilon2*forward(transpose(S$delta12)*S$delta12)*transpose(S$epsilon2))))+
         0.5*sum(diag(integral(transpose(S$delta12-delta12)*(S$delta12-delta12))))+
         0.5*sum(diag(integral(transpose(S$epsilon12-epsilon12)*(S$epsilon12-epsilon12))))
  )
}
