# Demographic stochasticity in discrete time
# Goals:
# 1. Learn a bit about Markov chains and
# create a simulator illustrated with a Davis
# rain model
# 2. Select a Parus major model via time series assuming process error (demographic stochasticity)
# 3. Simulate the model, find quasi-stationary distributions, and solve for extinction times. 
rm(list = ls())

# Goal 1. 

P=rbind(c(0.9055235 ,0.01363398, 0.08084251),
        c(0.5691057, 0.15718157, 0.27371274),
        c(0.4038997, 0.07195915, 0.52414113))
# create a simulator
# inputs: X1 initial state, P, Tf length of run
# output a vector X of all the states during run
MC=function(X1=1,P,Tf=100){
X=numeric(Tf)
X[1]=X1
k=dim(P)[1]
for(n in 2:Tf){
X[n]=sample(1:k,size = 1,prob=P[X[n-1],])
}
plot(X,type="l")
return(X)
}

out=MC(P=P,Tf=36500)
sum(out==1)/36500

stuff=eigen(t(P))
v=stuff$vectors[,1]

# Goal 2
N.data=c(21,30,31,32,20,20,21,31,27,24,49,27,41,51,86,43,39,54,46,45,32,33,38,30,46,36,39,26,33,14,30,26,19,27,39,28,46,39,30,33,48,26,35,18)

plot(N.data,type="b")

# after fitting some models we got

fit=function(N)N*exp(0.05)/(1+0.1*N)

# simulate this model

N1=N.data[1]
Tf=10000
N=numeric(Tf)
N[1]=21
for(t in 2:Tf){
N[t]=rpois(n = 1,lambda = fit(N[t-1]))
}
N
# create the Q matrix
chop=500
Q=matrix(NA,chop,chop)
for(i in 1:chop){
  Q[i,]=dpois(x=1:chop,lambda=fit(i))
}
stuff=eigen(t(Q))
v=Re(stuff$vectors[,1])
v=v/sum(v)
hist(N,freq=FALSE)
lines(v,lwd=2)


# find extinction times

# we already found one which is given 
# by 1/(1-lambda) where lambda is the leading
# eigenvalue of Q

1/(1-Re(stuff$values[1])) # "intrinsic" mean time to extinction

# how about the normal PVA estimates
# mean time to extinction given initially i individuals
# this is given by (Id-Q)^{-1}%*%ones 
# where ones is a column vectors of ones

extinct=solve(diag(chop)-Q)%*%matrix(1,chop,1)

plot(1:50,extinct[1:50],type="b")
