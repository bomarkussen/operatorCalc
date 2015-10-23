# Test of operatorCalc package
library(operatorCalc)

# mesh class
showMethods(classes="mesh")
(me <- mesh(seq(0,1,length.out=11)))

# matFct class
showMethods(classes="matFct")
(mF  <- matFct(mesh=me,f=array(rbind(sin(3.1415*me),cos(3.1415*me)),dim=c(1,2,11)),g=array(0,dim=c(1,2,10,0)),continuous=rep(TRUE,10)))
(mF2 <- matFct(mesh=me,f=array(rbind(me,-me),dim=c(1,2,11)),g=array(0,dim=c(1,2,10,0)),continuous=rep(TRUE,10)))
plot(mF)
plot(mF+mF2)
plot((mF+mF2)-mF)

# matFct: superSample
(mF3 <- matFct(mesh=mesh(c(0,1)),f=array(1:4,dim=c(2,2,2)),g=array(1:8,dim=c(2,2,1,2)),continuous=FALSE))
x <- seq(0.1,0.9,0.01)
matplot(x,t(matrix(evaluate(x,mF3),4,length(x))))
(mF4 <- superSample(x=mesh(sort(c(-1,2,sample(seq(0,1,0.005),20)))),mF3))
matlines(x,t(matrix(evaluate(x,mF3),4,length(x))))
summary(evaluate(x,mF3)-evaluate(x,mF4))
mF3@mesh
mF4@mesh
plot(mF4@mesh)
# speed
system.time(for (i in 1:1000) res <- evaluate(seq(0,1,0.01),mF3))
system.time(for (i in 1:1000) res <- superSample(seq(0,1,0.01),mF3))

# matFct: sum
x <- seq(-0.5,1.5,0.01)
matplot(x,t(matrix(evaluate(x,mF2+mF4[1,1:2]),2,length(x))))
matlines(x,t(matrix(evaluate(x,mF2),2,length(x))))
matlines(x,t(matrix(evaluate(x,mF4)[1,,],2,length(x))))
(mF2+mF4[1,1:2])@mesh
summary(evaluate(x,mF2)+evaluate(x,mF4)[1,,,drop=FALSE]-evaluate(x,mF2+mF4[1,1:2]))

# matFct: product
x <- seq(-0.1,1.1,0.01)
matplot(x,t(matrix(evaluate(x,mF2*mF4),2,length(x))))
matlines(x,t(matrix(evaluate(x,mF2),2,length(x))))
matlines(x,t(matrix(evaluate(x,mF4),4,length(x))))
(mF2*mF4)@mesh
for (i in sort(sample(1:length(x),40))) {
  cat("t=",x[i],":",evaluate(x,mF2)[,,i]%*%evaluate(x,mF4)[,,i]-evaluate(x,mF2*mF4)[,,i],"\n")
}
# speed
mF4 <- superSample(mF2@mesh,mF4)
mF2 <- superSample(mF4@mesh,mF2)
system.time(for (i in 1:1000) res1 <- mF2*mF4)

# matFct: forward
x <- seq(0,2,0.01)
mF <- matFct(mesh=mesh(c(0,1,2)),f=array(c(0,1,-2),dim=c(1,1,3)),g=array(rnorm(2),dim=c(1,1,2,1)),continuous=rep(FALSE,2))
plot(x,evaluate(x,mF))
plot(x,evaluate(x,forward(mF)))
plot(x,evaluate(x,forward(forward(mF))))
plot(x,evaluate(x,forward(forward(forward(mF)))))
# speed
system.time(for (i in 1:1000) res <- forward(mF4))

# matFct: forward
x <- seq(0,2,0.01)
mF <- matFct(mesh=mesh(c(0,1,2)),f=array(c(0,1,2),dim=c(1,1,3)),g=array(0,dim=c(1,1,2,0)),continuous=FALSE)
plot(x,evaluate(x,mF))
plot(x,evaluate(x,forward(mF)))
plot(x,evaluate(x,forward(forward(mF))))
plot(x,evaluate(x,forward(forward(forward(mF)))))
# speed
system.time(for (i in 1:1000) res <- backward(mF4))

# matFct: backward
x <- seq(0,1,0.01)
mF <- matFct(mesh=mesh(c(0,0.5,1)),f=array(c(0,0.5,1),dim=c(1,1,3)),g=array(0,dim=c(1,1,2,0)),continuous=rep(FALSE,2))
plot(x,evaluate(x,mF))
plot(x,evaluate(x,backward(mF)))
plot(x,evaluate(x,backward(backward(mF))))
plot(x,evaluate(x,backward(backward(backward(mF)))))
# speed
system.time(for (i in 1:1000) res <- backward(mF4))

# matFct: triangular
mF <- matFct(mesh=mesh(c(0,1)),f=array(c(1,1),dim=c(1,1,2)),g=array(0,dim=c(1,1,1,0)),continuous=FALSE)
x <- seq(0,1,0.01)
triangular(superSample(seq(0,1,0.5),mF))

# matFct: concatenation
mF <- matFct(mesh=mesh(c(0,1)),f=array(c(1,1),dim=c(1,1,2)),g=array(1,dim=c(1,1,1,1)),continuous=FALSE)
mF@f
forward(mF)@f
(forward(mF)|mF)@f
mF@g
forward(mF)@g
(forward(mF)|mF)@g
# speed
system.time(for (i in 1:1000) res <- mF4|mF4)


# --------------------------
# Approximation of sample
# --------------------------

tmp <- apply(0.01*array(rnorm(4*101),dim=c(101,2,2)),c(2,3),cumsum)
plot(approximate(seq(0,1,0.01),tmp,seq(1,101,10),order=0))
plot(approximate(seq(0,1,0.01),tmp,seq(1,101,10),order=1))
plot(approximate(seq(0,1,0.01),tmp,seq(1,101,10),order=2))

mF1 <- approximate(seq(0,1,0.01),tmp,seq(1,101,10),order=2)
mF2 <- approximate(c(0,0.5,1),array(rnorm(12),dim=c(3,2,2)),order=0)
plot(mF2)
