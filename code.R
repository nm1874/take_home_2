#1. 
A <- matrix(c(3,1,-1,-2),2); A
v0<- c(1,4)
P <- eigen(A)$vectors; P
PInv <- solve(P);
eigen(A)$values

Aexp <- function(t) P%*%diag(c(exp(eigen(A)$values[1]*t),exp(eigen(A)$values[2]*t)))%*%PInv

t<-1
Aexp(t)%*%v0

t<-2
Aexp(t)%*%v0

F <- function(v) c(3*v[1]-v[2], v[1]-2*v[2])
jacobian(F, Aexp(1)%*%v0)
jacobian(F, Aexp(2)%*%v0)

#2. 
#xz=y^2
#2(x+y)=3z
#Define a function from R^3 to R^3 to use in f(v) = 0
f <-function(v) c(v[1]^2+v[2]^2+v[3]^2-49,v[1]*v[3]-v[2]^2, 2*(v[1]+v[2])-3*v[3])
#Make a good initial guess
v0 <- c(3,4,5); f(v0)
#Calculate the 3 x 3 Jacobian matrix
A <- jacobian(f, v0); A
#Invert the Jacobian
AInv =solve(A); AInv
#Use the same update formula as in the single-variable case
v1 <- v0 - AInv%*%f(v0); v1; f(v1)
#Repeat to improve the approximation
A <- jacobian(f, v1)
v2 <- v1 - solve(A)%*%f(v1); v2; f(v2)

A <- jacobian(f, v2)
v3 <- v2 - solve(A)%*%f(v2); f(v3)
v3     #this is a very good approximate solution to our original equations


#3.
f<-function(v) c(sqrt((v[1]-v[3])^2 + (v[2]-v[4])^2)-5, sqrt((v[1])^2 + (v[6]-v[2])^2)-sqrt(20), sqrt((v[3])^2 + (v[6]-v[4])^2)-sqrt(85), sqrt((v[3]-v[5])^2 + (v[4])^2)-sqrt(26))

v0<- c(4,9,7,5,8,11)
f(v0)

Df<-jacobian(f, v0); Df

#Write DF = [A | B]
A <- Df[,1:4]; B<- Df[,5:6]; A; B
AInv <- solve(A)
#The derivative of g is -A^{-1)B
#Suppose we change z to 3.14. What x and y satisfy the constraints?
h <- c(0.3, -0.2)   #the increment to z
v <- v0 + c(-AInv%*%B%*%h,h); v

#A's new coordinates are approximately (4.048276,  8.896552)
#B's new coordinates are approximately (7.144828,  4.968966)


#4. 

library(numDeriv)

f <- function(x,y) x^3 - 8*x^2+y^2-x*y+20*x-3*y  #needed for contour lines
fVec <- function(v) v[1]^3-8*v[1]^2+v[2]^2-v[1]*v[2]+20*v[1]-3*v[2] #needed for numDeriv
#We can search for maxima and minima by plotting contour lines.
y <- x <- seq(1, 5, 0.2)
z <- outer(x,y,f)    #matrix of function values
filled.contour(x,y,z) #colors may help identify maximum/minimum.
contour(x,y,z)

#The closed contour for value 4 is suggestive. It is centered near (4,4).
x <- seq(2, 5, .2)
y <- seq(2, 5, 0.2)
z <- outer(x,y,f)    #matrix of function values
contour(x,y,z)

#It looks as though we have found a minimum near (4, 3.5).
x <- seq(3.5, 4.5, 0.02)
y <- seq(3, 4, 0.02)
z <- outer(x,y,f)    #matrix of function values
contour(x,y,z)
contour(x,y,z, levels = c(4.0,4.02, 4.04, 4.06, 4.08, 4.09, 4.11))

#Now we have a good approximation.
#Use the partial derivatives to get two equations in two unknowns
Df <- function(v) c(3*v[1]^2-16*v[1]-v[2]+20, 2*v[2]-v[1]-3)
v0 <- c(4, 3.5)
Df(v0)       #pretty close to zero
A <- jacobian(Df, v0); A
v1 <- v0 - solve(A)%*%Df(v0); v1; Df(v1)
grad(fVec, v1)  

H <- hessian(fVec, v0); H; det(H) #positive det confirms extremum
sum(H*diag(c(1,1)))   #positive trace confirms a minimum


#We have found a local max. Go back and look for another extremum.
#Replot the original contour lines
y <- x <- seq(1, 5, 0.2)
z <- outer(x,y,f)    #matrix of function values
filled.contour(x,y,z) #colors may help identify maximum/minimum.
contour(x,y,z)

#The 10 contour gets close to itself near (1.5,2)
abline(h=2, col = "green"); abline(v=1.5, col = "green")

#Zoom in for a closer look
x <- seq(1, 2, 0.2)
y <- seq(1, 3, 0.2)
z <- outer(x,y,f)    #matrix of function values
filled.contour(x,y,z)
contour(x,y,z)

#Something is happening for values between 1.6 and 2.3
contour(x,y,z, levels = seq(10, 10.6, 0.01))
#we have found a saddle point near (1.6, 2.3)
#Zoom in for a closer look.
x <- seq(1.5, 1.7, 0.01)
y <- seq(2.2, 2.4, 0.01)
z <- outer(x,y,f)    #matrix of function values
filled.contour(x,y,z)
contour(x,y,z)
#Something is happening for values between 4.56 and 4.58.
contour(x,y,z, levels = seq(10.3, 10.4, 0.0005))

#It looks as though we have found a saddle point near (1.56, 2.27)
v0 <- c(1.56, 2.27)
Df(v0)       #pretty close to zero.
A <- jacobian(Df, v0); A
v1 <- v0 - solve(A)%*%Df(v0); v1; Df(v1) #one iteration of Newton's method works wonders
f(v1[1], v1[2])      #this is the value of the function at the saddle point.
#Notice that there are larger and smaller function values nearby.
#This is not an extremum!

#Have a look at the Hessian matrix
grad(fVec,v1)       #confirms our critical point
H <- hessian(fVec, v1); H; det(H) #negative determinant confirms
