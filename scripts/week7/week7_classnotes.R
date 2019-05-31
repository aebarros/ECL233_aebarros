#week 7 class notes
rm(list=ls()) #clear env
library(bbmle)


#Demographic stocasticity in deiscrete time
#Goals:
#1) Learn a bit about Markov chains and create a simulator issustrated with a Davis rain model
#2) Select a Parus major model via time series assuming process error (Demographic stocasticity)
#3) Simulate the model, find quasi-stationary distributions, and solve for extinction times.

#-----------------
#Goal 1 Markov Chains in a finite state space
#state space S={1,2,....,k} in the model for rain, 3 states (didn't rain, very little rain, more than .01 in of rain)
#markov chain will be a random walk along the state space
#X[n] = state at time n
#P[ij]= P(X[n+1] = j|Xn=i)  this provides a k x k matrix (a probability of going from each state to another)
P=rbind(c(0.9055235 ,0.01363398, 0.08084251),
        c(0.5691057, 0.15718157, 0.27371274),
        c(0.4038997, 0.07195915, 0.52414113))
#three rain states 1,2,3 (0, <0.01, >0.01 inches of rain), matrix P is a probability of moving from one state to another (3x3 matrix)
#check to make sure rows of P add up to 1
rowSums(P) #all rows add to 1

#create a simulator
#inputs: X1 initial state, matrix P, tf=length of run
#output will be a vector x of all the states during the run to tf

MC=function(X1=1,P,tf=100){
  X=numeric(tf) #make vector to plug in
  X[1]=X1
  k=dim(P)[1] #gives us k as the number of dimensions in the matrix
  for(n in 2:tf){
    X[n]=sample(1:k,size=1,prob=P[X[n-1],])
    }
  plot(X,type="l") #gives us a plot of our rain states
  return(X)}
MC(P=P)
out=MC(P=P,tf=3650)
#interested in what fraction of days is dry?
sum(out==1)/3650 #~84%, this will be more stable at a longer time rur, here we are determining the limiting behavior of the Markov Chain. So if the MC, if P is irriducible (can get from any state to anyother state), then teh fraction of i's after n steps converges (as n goes to inf) to a vector which adds up to 1 and satisfies the vP=v.... [not sure if this is useful]

stuff=eigen(t(P))
v=stuff$vectors[,1]
v/sum(v) #this is normalizing the vector and tells us what fraction of each state we spend time at

#----------------
#Goal 2
N.data=c(21,30,31,32,20,20,21,31,27,24,49,27,41,51,86,43,39,54,46,45,32,33,38,30,46,36,39,26,33,14,30,26,19,27,39,28,46,39,30,33,48,26,35,18)

plot(N.data,type="both")
#try to fit a model to this data assuming there is not observation error (even though of course there is) and instead there is only demographic stocasticity, and not accounting for env stocasticity.
#Model of form: N[t+1] ~ Poisson(N[t]exp(a+b*N[t])) 
#the above: exp(a+b*N[t]) is the fitness of an individual, or the mean # of offspring
#so each individual then draws a mean from a Poisson distribution (which are positive whole integers)

#after fitting some models offscreen we got 

#fit=function(N)N*exp(1.15-.15*N) #can change a and b values to effect final plot
fit=function(N)N*exp(1.15-.15*N)
#simulate the model

N1=N.data[1]
tf=1000
N=numeric(tf)
N[1]=N1
for(t in 2:tf){
N[t]=rpois(n=1,lambda= fit(N[t-1]))
}
plot(N,type="l") #with these numbers, this population doesn't often go extinct

#with this model we are seeing LONG transients
#we will create a Q probability matrix (matrix without the 0 state?)
chop=200
Q=matrix(NA,chop,chop)
for(i in 1:chop){
  Q[i,]=dpois(x=1:chop,lambda=fit(i))
}

stuff=eigen(t(Q))
v=Re(stuff$vectors[,1])
v=v/sum(v)
v
sum(v) # =1 so it's right

hist(N,freq=FALSE,col=rgb(1,0,0,0.5))
lines(v,lwd=2)

#-----------
#Goal 3
#find extinction times
#we already found 1, where 1/(1-lambda) where lambda is the leading eigenvalue of Q

1/(1-Re(stuff$values[1])) #"intrinsic mean time to extinction" gives us expected time to extinction 700 million years
#here lambda is the probability of not going extinct which is VERY close to 1

#now want to find mean time to extinction given i initial individuals
#this is given by (Id-Q)^{-1}%*%ones where ones is a column vector of ones

extinct=solve(diag(chop)-Q)%*%matrix(1,chop,1) #gives inverse of ID matrix?
plot(1:100,extinct[1:100],type="b") #gives us a plot of mean time to extinction based on initial starting population

